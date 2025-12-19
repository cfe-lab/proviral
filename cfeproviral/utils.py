import csv
import logging
import os
import re
import typing
from typing import TextIO, Mapping, Dict, Set, List, Iterable, Tuple, Optional, Literal

import yaml
import json
import shutil
from cfeproviral.version import get_version, get_cfeintact_version
import pandas as pd
import glob
from pathlib import Path
import mappy
from aligntools import Cigar, CigarHit
from csv import DictWriter, DictReader
from itertools import groupby
from operator import itemgetter

logger = logging.getLogger(__name__)

Backend = Optional[Literal["CFEIntact", "HIVSeqinR"]]

LEFT_PRIMER_END = 666
RIGHT_PRIMER_START = 9604


def load_yaml(afile):
    with open(afile) as f:
        return yaml.load(f, Loader=yaml.SafeLoader)


def dump_yaml(data, afile):
    with open(afile, 'w') as o:
        yaml.dump(data, o)


def reverse_and_complement(seq):
    return ''.join(complement_dict[nuc] for nuc in reversed(seq))


def write_annot(adict, filepath):
    """
        adict: {gene (str): [start (int), stop (int)], ...}
    """
    with open(filepath, 'w', newline='') as o:
        fieldnames = ['gene', 'start', 'stop']
        writer = DictWriter(o, fieldnames)
        writer.writeheader()
        for gene, (start, stop) in adict.items():
            writer.writerow({
                'gene': gene.upper(),
                'start': start,
                'stop': stop
            })


def write_fasta(adict, filepath_or_fileobject):
    """Writes a fasta file from a dictionary

    Args:
        adict (dict): A dictionary where the keys are header/sequence names and the values are the sequences
        filepath_or_fileobject (Path/File): A path or file-like object

    Returns:
        Path/File: Path or file-like object
    """
    try:
        with open(filepath_or_fileobject, 'w') as o:
            for key, value in adict.items():
                o.write(f'>{key}\n{value}\n')
    except TypeError:
        for key, value in adict.items():
            filepath_or_fileobject.write(f'>{key}\n{value}\n')
    return filepath_or_fileobject


def read_csv(csvfile):
    with open(csvfile) as f:
        reader = DictReader(f)
        for row in reader:
            yield row


def csv_to_bed(csvfile, target_name='HXB2', offset_start=0, offset_stop=0):
    with open(csvfile) as f:
        reader = DictReader(f)
        with open(f'{csvfile}.bed', 'w', newline='') as o:
            fieldnames = [
                'chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand',
                'thickStart', 'thickEnd'
            ]
            writer = DictWriter(o, fieldnames, delimiter='\t')
            for row in reader:
                writer.writerow({
                    'chrom': target_name,
                    'chromStart': int(row['start']) + offset_start,
                    'chromEnd': int(row['stop']) + offset_stop,
                    'name': row['gene'].upper(),
                    'score': 0,
                    'strand': '+',
                    'thickStart': int(row['start']) + offset_start,
                    'thickEnd': int(row['stop']) + offset_stop
                })


class MyContextManager(object):
    def __init__(self, file_name, method):
        self.file_obj = open(file_name, method)

    def __enter__(self):
        return self.file_obj

    def __exit__(self, type, value, traceback):
        self.file_obj.close()


def handle_file_or_fileobject(file_or_fileobject, method):
    try:
        file_or_fileobject.read()
        file_or_fileobject.seek(0)
        return file_or_fileobject
    except AttributeError:
        return MyContextManager(file_or_fileobject, method)


def read_fasta(fasta_file):
    with handle_file_or_fileobject(fasta_file, 'r') as fasta:
        name, seq = None, []
        for line in fasta:
            line = line.rstrip()
            if line.startswith(">"):
                if name:
                    yield (name, ''.join(seq))
                name, seq = line, []
            else:
                seq.append(line)
        if name:
            yield (name, ''.join(seq))


def get_genes(annot, pos, offset_start=0, offset_end=0, offset=None):
    genes = []
    if offset is not None:
        offset_start = offset_end = offset
    """
        annot (dict) => {gene: [start, end], ...}
        ref (str),
        pos (int)
    """
    for gene in annot:
        if ((annot[gene][0] - offset_start) <= pos <=
            (annot[gene][1] - offset_end)):
            genes.append(gene)
    return genes


def modify_reference(refseq):
    newseq = refseq

    # Second round primers (regex for mixture positions)
    fwd = 'GCGCCCGAACAGGGAC.TGAAA.CGAAAG'
    rev = 'TAAGCCTCAATAAAGCTTGCCTTGAGTGC'
    fwd_pos = re.search(fwd, newseq).start()
    rev_pos = [x for x in re.finditer(rev, newseq)][1].start()

    fwd_offset = fwd_pos + len(fwd)

    # Trim primers
    newseq = newseq[fwd_offset:rev_pos]

    assert newseq[:10] == 'GGAAACCAGA'

    # Reading frame correction (remove 1 nucleotide)
    newseq = newseq[:5107] + newseq[5108:]
    assert newseq[5102:5112] == 'CATTTCAGAA'

    # Fix premature stop codon
    pos = 8499
    assert newseq[pos] == 'A'
    newseq = newseq[:pos] + 'G' + newseq[pos + 1:]
    assert newseq[pos] == 'G'

    return newseq


def modify_annot(annot):
    newannot = {}
    for gene, (start, stop) in annot.items():
        if gene in genes_of_interest:
            newannot[gene] = [start, stop]
    # newannot = dict(annot)
    # Offset by second round fwd primer trim (666 trimmed from start)
    offset = -666
    for gene, (start, stop) in newannot.items():
        start += offset
        stop += offset
        if start <= 0:
            continue
        # Account for 1 nucleotide removal at index 5107 (pos 5108)
        if start >= 5108:
            start -= 1
            stop -= 1
        elif stop >= 5108:
            stop -= 1
        # Convert to 0 base
        newannot[gene] = [start - 1, stop - 1]

    return newannot


def splice_genes(query, target, cigar_hits, annotation):
    """
    Extract gene coordinates from aligned sequences using CigarHit objects.
    
    Uses aligntools' coordinate mapping for direct global coordinate translation.
    
    Args:
        query: query sequence string
        target: reference sequence string
        cigar_hits: List of CigarHit objects from alignment
        annotation: annotation dictionary mapping positions to genes
    
    Returns:
        dict: gene name -> [start_pos, end_pos] in query coordinates
    """
    results = {}
    
    for hit in cigar_hits:
        # CigarHit.coordinate_mapping works with GLOBAL coordinates
        mapping = hit.coordinate_mapping
        
        # For each gene, map its reference coordinates to query coordinates
        for gene, (gene_start, gene_end) in annotation.items():
            # Check if gene overlaps with the aligned region
            if gene_end < hit.r_st or gene_start > hit.r_ei:
                continue
            
            # Clip gene to aligned region
            clipped_start = max(gene_start, hit.r_st)
            clipped_end = min(gene_end, hit.r_ei)
            
            # Map reference positions to query positions
            # We must iterate through all reference positions because some may be
            # deleted (not present in the alignment). We collect all mapped query
            # positions and take the min/max to get the gene boundaries.
            query_positions = []
            for ref_pos in range(clipped_start, clipped_end + 1):
                query_pos = mapping.ref_to_query.get(ref_pos)
                if query_pos is not None:
                    query_positions.append(query_pos)
            
            if query_positions:
                query_start = min(query_positions)
                query_end = max(query_positions)
                if gene not in results:
                    results[gene] = [query_start, query_end]
                else:
                    results[gene][0] = min(results[gene][0], query_start)
                    results[gene][1] = max(results[gene][1], query_end)
                
    return results


def coords_to_genes(coords, query):
    genes = {
        gene: query[coords[0]:coords[1] + 1]
        for gene, coords in coords.items()
    }
    return genes


def splice_aligned_genes(query, target, cigar_hits, annotation):
    """
    Extract gene coordinates and aligned sequences using CigarHit objects.
    
    Args:
        query: query sequence string
        target: reference sequence string
        cigar_hits: List of CigarHit objects from alignment
        annotation: annotation dictionary mapping positions to genes
    
    Returns:
        tuple: (results dict, sequences dict) 
            - results: gene name -> [start_pos, end_pos]
            - sequences: gene name -> list of nucleotides
    """
    results = {}
    sequences = {}
    
    for cigar_hit in cigar_hits:
        # Use the coordinate mapping from CigarHit
        mapping = cigar_hit.cigar.coordinate_mapping
        
        # Iterate through reference positions in the aligned region
        # Note: CigarHit uses inclusive ranges
        for ref_pos in range(cigar_hit.r_st, cigar_hit.r_ei + 1):
            # Skip if reference position is outside target sequence
            if ref_pos >= len(target):
                logger.warning(f'{ref_pos} not in range of target')
                break
            
            target_nuc = target[ref_pos]
            
            # Reference position relative to alignment start
            ref_pos_in_alignment = ref_pos - cigar_hit.r_st
            
            # Get genes at this reference position
            genes = get_genes(annotation, ref_pos)
            
            # Check if this reference position maps to a query position
            if ref_pos_in_alignment in mapping.ref_to_query:
                # Position maps - there's a nucleotide here
                query_pos_in_alignment = mapping.ref_to_query[ref_pos_in_alignment]
                query_pos = cigar_hit.q_st + query_pos_in_alignment
                query_nuc = query[query_pos]
                
                print(ref_pos, query_pos, target_nuc, query_nuc, 
                      target_nuc == query_nuc, genes)
                
                for gene in genes:
                    if gene not in results:
                        results[gene] = [query_pos, query_pos]
                        sequences[gene] = [query_nuc]
                    else:
                        results[gene][0] = min(results[gene][0], query_pos)
                        results[gene][1] = max(results[gene][1], query_pos)
                        sequences[gene].append(query_nuc)
            else:
                # Position doesn't map - it's a deletion in the query
                print(f'Deletion at ref pos {ref_pos}')
                for gene in genes:
                    if gene in sequences:
                        sequences[gene].append('-')
            
            print('=' * 50)
        print('new alignment row'.center(50, '~'))
    
    return results, sequences


# Removes all files in a directory
def clean_dir(directory):
    for filename in os.listdir(directory):
        file_path = os.path.join(directory, filename)
        try:
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)
        except Exception as e:
            print('Failed to delete %s. Reason: %s' % (file_path, e))


def mappy_hit_to_cigar_hit(hit, query_seq):
    """
    Convert a mappy alignment hit to aligntools CigarHit object.
    
    Args:
        hit: mappy alignment hit object
        query_seq: the query sequence string (unused, kept for compatibility)
    
    Returns:
        CigarHit: aligntools CigarHit object
    """
    cigar = Cigar(hit.cigar)
    
    # Convert from half-open [r_st, r_en) to inclusive [r_st, r_ei]
    return CigarHit(
        cigar=cigar,
        r_st=hit.r_st,
        r_ei=hit.r_en - 1,
        q_st=hit.q_st,
        q_ei=hit.q_en - 1
    )


def align(target_seq,
          query_seq,
          query_name,
          outdir=Path(os.getcwd()).resolve()):
    """
    Align query sequence to target sequence using mappy.
    
    Uses mappy (Python binding for minimap2) and returns aligntools CigarHit objects.
    
    Args:
        target_seq: reference sequence to align to
        query_seq: query sequence to align
        query_name: name for the query sequence
        outdir: output directory for alignment files (for logging)
    
    Returns:
        List[CigarHit]: List of CigarHit objects, or empty list if alignment fails
    """
    outdir = outdir / query_name
    if os.path.isdir(outdir):
        shutil.rmtree(outdir)
    os.makedirs(outdir)
    
    # Write the query and target fasta (for debugging/reference)
    write_fasta({query_name: query_seq}, outdir / 'query.fasta')
    write_fasta({'MOD_HXB2': target_seq}, outdir / 'target.fasta')
    
    log_path = outdir / 'minimap2.log'
    
    try:
        # Create aligner with no preset (matches minimap2 -a default behavior)
        aligner = mappy.Aligner(seq=target_seq)
        
        if not aligner:
            raise RuntimeError("Failed to create mappy aligner")
        
        # Perform alignment
        hits = list(aligner.map(query_seq))
        
        # Write log file (for debugging)
        with log_path.open('w') as log_file:
            log_file.write("mappy alignment completed\n")
            log_file.write(f"Query: {query_name} ({len(query_seq)} bp)\n")
            log_file.write(f"Reference: MOD_HXB2 ({len(target_seq)} bp)\n")
            log_file.write(f"Hits: {len(hits)}\n")
        
        # Convert mappy hits to CigarHit objects
        cigar_hits = [mappy_hit_to_cigar_hit(hit, query_seq) for hit in hits]
        
        return cigar_hits
        
    except Exception as e:
        # Write error to log
        with log_path.open('w') as log_file:
            log_file.write(f"mappy alignment failed: {str(e)}\n")
        logger.error('Alignment with mappy failed! Details in %s.', log_path)
        return []

CFEINTACT_ERRORS_TABLE = [
    'UnknownNucleotide',
    'NonHIV',
    'LongDeletion',
    'InternalInversion',
    'Scramble',
    'APOBECHypermutation',
    'MajorSpliceDonorSiteMutated',
    'PackagingSignalDeletion',
    'PackagingSignalNotComplete',
    'RevResponseElementDeletion',
    'MisplacedORF',
    'WrongORFNumber',
    'Deletion',
    'Insertion',
    'InternalStop',
    'Frameshift',
    'MutatedStopCodon',
    'MutatedStartCodon',
    'SequenceDivergence',
    ]

def iterate_cfeintact_verdicts_1(directory: Path, intact: Set[str] = set()) -> Iterable[Tuple[str, str]]:
    intact = set()

    def get_verdict(SEQID: str, all_defects) -> Tuple[str, str]:
        if all_defects:
            ordered = sorted(all_defects, key=CFEINTACT_ERRORS_TABLE.index)
            verdict = ordered[0]
        else:
            verdict = "Intact"

        return (SEQID, verdict)

    with open(os.path.join(directory, 'holistic.csv'), 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row["intact"] == "True":
                intact.add(row["qseqid"])
                SEQID = row["qseqid"]
                yield get_verdict(SEQID, all_defects=[])

    with open(os.path.join(directory, 'defects.csv'), 'r') as f:
        reader = csv.DictReader(f)
        grouped = groupby(reader, key=itemgetter('qseqid'))
        for sequence_name, defects in grouped:
            if sequence_name not in intact:
                all_defects = [defect['code'] for defect in defects]
                yield get_verdict(sequence_name, all_defects)


def iterate_cfeintact_verdicts(outpath: Path) -> Iterable[Tuple[str, str]]:
    intact: Set[str] = set()

    for directory in outpath.glob('cfeintact*'):
        yield from iterate_cfeintact_verdicts_1(directory, intact)


def get_cfeintact_verdicts(name, outpath):
    column_names = ['SEQID', 'MyVerdict']
    data = iterate_cfeintact_verdicts(outpath)
    return pd.DataFrame(data, columns=column_names)


def iterate_hivseqinr_verdicts_1(directory: Path) -> Iterable[Tuple[str, str]]:
    path = directory / "Results_Final" / "Output_MyBigSummary_DF_FINAL.csv"
    if not path.is_file():
        logger.error("Missing HIVSeqinR result in %r.", str(path))
        return

    with path.open() as fd:
        reader = csv.DictReader(fd)
        for row in reader:
            yield (row["SEQID"], row["MyVerdict"])


def iterate_hivseqinr_verdicts(outpath: Path) -> Iterable[Tuple[str, str]]:
    for directory in outpath.glob('hivseqinr*'):
        yield from iterate_hivseqinr_verdicts_1(directory)


def get_hivseqinr_verdicts(name, outpath):
    column_names = ['SEQID', 'MyVerdict']
    data = iterate_hivseqinr_verdicts(outpath)
    return pd.DataFrame(data, columns=column_names)


def generate_table_precursor(name, outpath, add_columns=None):
    # Output csv
    precursor_path: Path = outpath / 'table_precursor.csv'

    # Load filtered sequences
    filtered_path = outpath / (name + '_filtered.csv')
    filtered = pd.read_csv(filtered_path)
    # Load hivseqinr data or Cfeintact results

    if any(outpath.glob('cfeintact*')):
        results = get_cfeintact_verdicts(name, outpath)
    elif any(outpath.glob('hivseqinr*')):
        results = get_hivseqinr_verdicts(name, outpath)
    else:
        raise RuntimeError("Neither CFEIntact nor HIVSeqinR directory exists.")

    try:
        # Assign new columns based on split
        results[['name', 'sample', 'reference',
                'seqtype']] = results['SEQID'].str.split('::', expand=True)
        # Merge
        merged = results.merge(filtered, on='sample')
    except ValueError:
        with precursor_path.open('w') as output_file:
            writer = DictWriter(output_file,
                                ['sample', 'sequence', 'MyVerdict'] +
                                genes_of_interest)
            writer.writeheader()
        return precursor_path

    for gene in genes_of_interest:
        merged[gene] = None

    data = {}
    for index, row in merged.iterrows():
        folder = outpath / row['sample']
        genes_fasta_path = folder / 'genes.fasta'
        if not os.path.isfile(genes_fasta_path):
            print(f'No genes for {row["sample"]}')
            continue
        genes_fasta = read_fasta(genes_fasta_path)
        genes = dict([x for x in genes_fasta])
        for gene in genes_of_interest:
            try:
                seq = genes[f'>{gene}']
            except KeyError:
                seq = None
            data.setdefault(gene, []).append(seq)
    for gene, seqs in data.items():
        merged[gene] = seqs

    if add_columns:
        for key, val in add_columns.items():
            merged[key] = val
    # Version columns should come from the filtered dataframe
    # If they're not there (shouldn't happen with updated code), add them
    if 'cfeproviral_version' not in merged.columns:
        merged['cfeproviral_version'] = get_version()
    if 'cfeintact_version' not in merged.columns:
        merged['cfeintact_version'] = get_cfeintact_version()
    if 'micall_version' not in merged.columns:
        merged['micall_version'] = None
    
    if not results.empty:
        merged[['sample', 'sequence', 'MyVerdict'] + genes_of_interest + ['cfeproviral_version', 'cfeintact_version', 'micall_version']].to_csv(
            precursor_path, index=False)
    else:
        merged[['sample', 'sequence'] + genes_of_interest + ['cfeproviral_version', 'cfeintact_version', 'micall_version']].to_csv(precursor_path,
                                                                  index=False)
    return precursor_path


def generate_table_precursor_2(hivseqinr_resultsfile, filtered_file,
                               genes_file, table_precursorfile):
    try:
        seqinr = pd.read_csv(hivseqinr_resultsfile)
    except FileNotFoundError:
        logger.error('hivseqinr could not produce results!')
        # Create an empty results file
        table_precursorfile.touch()
        return False
    # Assign new columns based on split
    # Make sure this matches the join in primer_finder run()
    seqinr[['reference', 'seqtype']] = seqinr['SEQID'].str.split('::',
                                                                 expand=True)

    # Load filtered sequences
    filtered = pd.read_csv(filtered_file)

    # Merge
    merged = seqinr.merge(filtered,
                          left_index=True,
                          right_index=True,
                          how='outer')
    for gene in genes_of_interest:
        merged[gene] = None

    genes_fasta = read_fasta(genes_file)
    genes = dict([x for x in genes_fasta])
    for gene in genes_of_interest:
        try:
            seq = genes[f'>{gene}']
        except KeyError:
            seq = None
        merged[gene] = seq

    # Output csv
    merged['cfeproviral_version'] = get_version()
    merged['cfeintact_version'] = get_cfeintact_version()
    # micall_version should come from filtered dataframe
    # If not present, set to None
    if 'micall_version' not in merged.columns:
        merged['micall_version'] = None
    merged[['sequence', 'MyVerdict'] + genes_of_interest + ['cfeproviral_version', 'cfeintact_version', 'micall_version']].to_csv(
        table_precursorfile, index=False)
    return table_precursorfile


def get_softclipped_region(query, cigar_hits):
    """
    Extract the soft-clipped region from alignment if present.
    
    Mappy doesn't include soft-clipping in the CIGAR string - it's implicit.
    If the alignment starts at q_st > 0, those bases are soft-clipped.
    
    Args:
        query: query sequence string
        cigar_hits: List of CigarHit objects
    
    Returns:
        str or None: soft-clipped sequence if present, None otherwise
    """
    if not cigar_hits:
        logger.warning('No alignment hits!')
        return None
    
    first_hit = cigar_hits[0]
    
    # Soft-clipped bases are before q_st
    if first_hit.q_st > 0:
        return query[:first_hit.q_st]
    
    return None


def sequence_to_coords(query, target, cigar_hits, annotations):
    """
    Convert sequence to coordinates using soft-clipped region.
    
    Args:
        query: query sequence string
        target: reference sequence string
        cigar_hits: List of CigarHit objects
        annotations: annotation dictionary
    
    Returns:
        dict or None: gene coordinates
    """
    softclip = get_softclipped_region(query, cigar_hits)
    if softclip is None:
        return None
    import cfeproviral.probe_finder
    finder = cfeproviral.probe_finder.ProbeFinder(softclip, target)
    if not finder.valid:
        return None
    # query_match = target[finder.start:len(finder.contig_match)]
    target_match = finder.contig_match
    matchlen = len(target_match)
    coords = {}
    for pos in range(finder.start, matchlen - 1):
        genes = get_genes(annotations, pos)
        for gene in genes:
            coords.setdefault(gene, [pos, pos])[1] += 1
    return coords


def merge_coords(coords1, coords2):
    """Return new coordinates based on merging coords1 and coords2

    Args:
        coords1 (dict): A dictionary of coords with key=gene_name, value=[start, end]
        coords2 (dict): Same as coords1
    """
    new_coords = {}
    for gene in coords2:
        if gene not in coords1:
            new_coords[gene] = coords2[gene][:]
            continue
        new_coords[gene] = [
            min(coords2[gene][0], coords1[gene][0]),
            max(coords1[gene][1], coords2[gene][1])
        ]
    for gene in coords1:
        if gene not in coords2:
            new_coords[gene] = coords1[gene][:]
            continue
        new_coords[gene] = [
            min(coords1[gene][0], coords2[gene][0]),
            max(coords2[gene][1], coords1[gene][1])
        ]
    return new_coords


def filter_valid(df):
    # Remove any row that has no errors
    filtered = df[(~df['error'].isna())
                  | (~df['fwd_error'].isna())
                  | (~df['rev_error'].isna())].copy()
    # Remove contig not max errors and V3 errors
    filtered = filtered[(filtered['error'] != 'contig not MAX')
                        & (filtered['error'] != 'is V3 sequence')]
    # Set error field for duplicates
    #filtered.loc[filtered.duplicated(subset='sample', keep=False),
    #             'error'] = 'duplicate'
    filtered = filtered[(~filtered['reference'].str.contains('reverse'))
                        & (~filtered['reference'].str.contains('unknown'))]
    return filtered


def genFailureSummary(contigs_df, conseqs_df, outpath):
    filtered_contigs = filter_valid(contigs_df)
    filtered_conseqs = filter_valid(conseqs_df)
    contigs_simple = filtered_contigs[[
        'sample', 'run_name', 'reference', 'error', 'fwd_error', 'rev_error',
        'cfeproviral_version', 'cfeintact_version', 'micall_version'
    ]]
    conseqs_simple = filtered_conseqs[[
        'sample', 'run_name', 'reference', 'error', 'fwd_error', 'rev_error',
        'cfeproviral_version', 'cfeintact_version', 'micall_version'
    ]]
    concat = pd.concat([contigs_simple, conseqs_simple])
    outfile = outpath / 'failure_summary.csv'
    concat.to_csv(outfile, index=False)
    return outfile


# Checks if a value is nan, handles strings
def isNan(num):
    return num != num


def get_samples_from_cascade(cascade_csv: typing.IO,
                             default_sample_name: Optional[str] = None):
    all_samples = {}
    reader = DictReader(cascade_csv)
    if reader.fieldnames and 'sample' in reader.fieldnames:
        for row in reader:
            all_samples[row['sample']] = {
                'remap': int(row['remap'])
            }
        return all_samples
    rows = list(reader)
    assert len(rows) == 1, len(rows)
    remap_count = int(rows[0]['remap'])
    return {default_sample_name: {'remap': remap_count}}


## Define some variables
cwd = Path(os.path.realpath(__file__)).parent

# Define genes of interest
genes_of_interest = load_yaml(cwd / 'genes_of_interest.yaml')

# Define some global variables
mixture_dict = {
    'W': 'AT',
    'R': 'AG',
    'K': 'GT',
    'Y': 'CT',
    'S': 'CG',
    'M': 'AC',
    'V': 'AGC',
    'H': 'ATC',
    'D': 'ATG',
    'B': 'TGC',
    'N': 'ATGC',
    '-': 'ATGC'
}

complement_dict = {
    'A': 'T',
    'C': 'G',
    'G': 'C',
    'T': 'A',
    'W': 'S',
    'R': 'Y',
    'K': 'M',
    'Y': 'R',
    'S': 'W',
    'M': 'K',
    'B': 'V',
    'D': 'H',
    'H': 'D',
    'V': 'B',
    '*': '*',
    'N': 'N',
    '-': '-'
}

hxb2 = next(read_fasta(cwd / 'hxb2.fasta'))[1]
mod_hxb2 = modify_reference(hxb2)

annot = {
    x['gene']: [int(x['start']), int(x['stop'])]
    for x in read_csv(cwd / 'annot.csv')
}

mod_annot = modify_annot(annot)

# Note these are 1-based indicies
primers = {
    'fwd': {
        'seq': 'GCGCCCGAACAGGGACYTGAAARCGAAAG',
        'nomix': 'GCGCCCGAACAGGGACCTGAAAGCGAAAG',
        # convert to 0-base index
        'hxb2_start': 638 - 1,
        'hxb2_end': 666,
        'direction': 'fwd'
    },
    'rev': {
        'seq': 'TAAGCCTCAATAAAGCTTGCCTTGAGTGC',
        'nomix': 'TAAGCCTCAATAAAGCTTGCCTTGAGTGC',
        # convert to 0-base index
        'hxb2_start': 9604 - 1,
        'hxb2_end': 9632,
        'direction': 'rev'
    }
}
