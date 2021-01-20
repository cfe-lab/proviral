import math
import os
import re
from typing import Dict
import yaml
import shutil
import subprocess as sp
import pandas as pd
import glob
import numpy as np
from pathlib import Path
from csv import DictWriter, DictReader
from logger import logger


def load_yaml(afile):
    with open(afile) as f:
        return yaml.load(f, Loader=yaml.SafeLoader)


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


def write_fasta(adict, filename):
    with open(filename, 'w') as o:
        for key, value in adict.items():
            o.write(f'>{key}\n{value}\n')
    return filename


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


def split_cigar(row):
    pattern = re.compile('(\d+)([A-Z])')
    cigar = re.findall(pattern, row[5])
    return cigar


def read_fasta(fasta_file):
    with open(fasta_file) as fasta:
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


#     newannot = dict(annot)
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


def splice_genes(query, target, samfile, annotation):
    results = {}
    for i, row in samfile.iterrows():
        # Subtract 1 to convert target position to zero-base
        target_pos = int(row[3]) - 1
        query_pos = None
        for size, op in row['cigar']:
            size = int(size)
            # logger.debug(f'size: {size}, op: {op}')
            # logger.debug(f'target_pos: {target_pos}, query_pos: {query_pos}')
            # If the first section is hard-clipped the query should start at
            # the first non-hard-clipped-based. The target should also be offset
            # by the size of the hard-clip
            if op == 'H' and query_pos is None:
                query_pos = size
                # logger.debug('='*50)
                continue
            elif query_pos is None:
                query_pos = 0
            if op == 'S':
                query_pos += size
                # logger.debug('='*50)
                continue
            elif op in ('M', '=', 'X'):
                for j in range(size):
                    try:
                        target_nuc = target[target_pos]
                    except IndexError:
                        logger.warning(f'{target_pos} not in range of target')
                        break
                    query_nuc = query[query_pos]
                    match = (target_nuc == query_nuc)
                    genes = get_genes(annotation, target_pos)
                    # logger.debug(target_pos, query_pos, target_nuc, query_nuc, match, genes)
                    for gene in genes:
                        if gene not in results:
                            results[gene] = [query_pos, query_pos]
                        elif query_pos > results[gene][1]:
                            results[gene][1] = query_pos
                    query_pos += 1
                    target_pos += 1
            elif op == 'D':
                target_pos += size
                # logger.debug('='*50)
                continue
            elif op == 'I':
                query_pos += size
                # logger.debug('='*50)
                continue
            # logger.debug('='*50)
        # logger.debug('new alignment row'.center(50, '~'))
    return results


def coords_to_genes(results, query):
    genes = {
        gene: query[coords[0]:coords[1] + 1]
        for gene, coords in results.items()
    }
    return genes


# This function has not been tested yet
def get_sequences(query, target, samfile, annotation):
    results = {}
    sequences = {}
    for i, row in samfile.iterrows():
        # Subtract 1 to convert target position to zero-base
        target_pos = int(row[3]) - 1
        genes = get_genes(annotation, target_pos)
        query_pos = None
        for size, op in row['cigar']:
            print(f'size: {size}, op: {op}')
            print(f'target_pos: {target_pos}, query_pos: {query_pos}')
            size = int(size)
            # If the first section is hard-clipped the query should start at
            # the first non-hard-clipped-based. The target should also be offset
            # by the size of the hard-clip
            if op == 'H' and query_pos is None:
                query_pos = size
                print('=' * 50)
                continue
            elif query_pos is None:
                query_pos = 0
            if op == 'S':
                query_pos += size
                print('=' * 50)
                continue
            elif op in ('M', '=', 'X'):
                for i in range(size):
                    try:
                        target_nuc = target[target_pos]
                    except IndexError:
                        print(f'{target_pos} not in range of target')
                        break
                    query_nuc = query[query_pos]
                    match = (target_nuc == query_nuc)
                    genes = get_genes(annotation, target_pos)
                    print(target_pos, query_pos, target_nuc, query_nuc, match,
                          genes)
                    for gene in genes:
                        if gene not in results:
                            results[gene] = [query_pos, query_pos]
                            sequences[gene] = [query_nuc]
                        elif query_pos > results[gene][1]:
                            results[gene][1] = query_pos
                            sequences[gene].append(query_nuc)
                    query_pos += 1
                    if op != 'I':
                        target_pos += 1
            elif op == 'D':
                target_pos += size
                for i in range(size):
                    sequences[gene].append('-')
                print('=' * 50)
                continue
            elif op == 'I':
                query_pos += size
                print('=' * 50)
                continue
            print('=' * 50)
        print('new alignment row'.center(50, '~'))
    return results, sequences


def load_samfile(samfile_path):
    result = pd.read_table(samfile_path, skiprows=2, header=None)
    result['cigar'] = result.apply(split_cigar, axis=1)
    return result


def aligner_available(aligner_path='minimap2'):
    cmd = [aligner_path]
    try:
        process = sp.run(cmd)
        if process.returncode == 0:
            return True
    except Exception:
        return False


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


def align(target_seq,
          query_seq,
          query_name,
          outdir=Path(os.getcwd()).resolve(),
          aligner_path='minimap2'):
    if not aligner_available(aligner_path):
        return None
    outdir = outdir / query_name
    if os.path.isdir(outdir):
        shutil.rmtree(outdir)
    os.makedirs(outdir)
    # Write the query fasta
    query_fasta_path = write_fasta({query_name: query_seq},
                                   outdir / 'query.fasta')
    # Write the target fasta
    target_fasta_path = write_fasta({'MOD_HXB2': target_seq},
                                    outdir / 'target.fasta')
    cmd = [aligner_path, '-a', target_fasta_path, query_fasta_path]
    alignment_path = outdir / 'alignment.sam'
    with open(alignment_path, 'w') as alignment:
        process = sp.run(cmd, stdout=alignment, errors=True)
    if process.returncode != 0:
        raise UserWarning(f'Alignment failed! Error: {process.stderr}')
    else:
        return alignment_path


def generate_table_precursor(name, outpath):
    # Load filtered sequences
    filtered_path = outpath / f'{name}_filtered.csv'
    filtered = pd.read_csv(filtered_path)
    # Load hivseqinr data
    seqinr_paths = glob.glob(
        str(outpath / 'hivseqinr*' / 'Results_Final' /
            'Output_MyBigSummary_DF_FINAL.csv'))
    parts = []
    for path in seqinr_paths:
        if not os.path.isfile(path):
            continue
        part = pd.read_csv(path)
        parts.append(part)
    # seqinr = pd.read_csv(seqinr_path)
    try:
        seqinr = pd.concat(parts)
        # Assign new columns based on split
        seqinr[['name', 'sample', 'reference',
                'seqtype']] = seqinr['SEQID'].str.split('::', expand=True)
        # Merge
        merged = seqinr.merge(filtered, on='sample')
    except ValueError:
        print('No hivseqinr runs found')
        seqinr = None
        merged = filtered

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

    # Output csv
    outfile = outpath / 'table_precursor.csv'
    if parts:
        merged[['sample', 'sequence', 'MyVerdict'] + genes_of_interest].to_csv(
            outfile, index=False)
    else:
        merged[['sample', 'sequence'] + genes_of_interest].to_csv(outfile,
                                                                  index=False)
    return outfile


def get_softclipped_region(query, alignment):
    try:
        size, op = alignment.iloc[0]['cigar'][0]
    except IndexError:
        logger.warning('No alignment!')
        return
    if op != 'S':
        logger.warning('Alignment does not start with softclip')
        return
    size = int(size)
    return query[:size]


def sequence_to_coords(query, target, alignment_path, annot):
    aln = load_samfile(alignment_path)
    softclip = get_softclipped_region(query, aln)
    if softclip is None:
        return
    import probe_finder
    finder = probe_finder.ProbeFinder(softclip, target)
    if not finder.valid:
        return None
    # query_match = target[finder.start:len(finder.contig_match)]
    target_match = finder.contig_match
    matchlen = len(target_match)
    coords = {}
    for pos in range(finder.start, matchlen - 1):
        genes = get_genes(annot, pos)
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
    #    import pdb; pdb.set_trace()
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
        'sample', 'run_name', 'reference', 'error', 'fwd_error', 'rev_error'
    ]]
    conseqs_simple = filtered_conseqs[[
        'sample', 'run_name', 'reference', 'error', 'fwd_error', 'rev_error'
    ]]
    concat = pd.concat([contigs_simple, conseqs_simple])
    outfile = outpath / 'failure_summary.csv'
    concat.to_csv(outfile, index=False)
    return outfile


# Checks if a value is nan, handles strings
def isNan(num):
    return num != num


# Given a sample name, return whether or not the sample is proviral-related
def is_proviral(sample_name):
    if (('GAGGAG' in sample_name) or ('VIR' in sample_name) or
        ('NEF-HIV') in sample_name) or ('V3LOOP' in sample_name) or (
            'HLA' in sample_name) or ('HCV' in sample_name):
        return False
    return True


def convert_none(data, translation='none'):
    if data is None:
        return translation
    return data


def genOutcomeSummary(contigs_df, conseqs_df, outpath):
    data = {}

    contigs_df = contigs_df.where(pd.notnull(contigs_df), None)
    conseqs_df = conseqs_df.where(pd.notnull(conseqs_df), None)

    max_failed = 0

    # Go through all the conseqs
    for index, row in conseqs_df.iterrows():
        seqtype = 'conseq'
        sample = row['sample']
        # passed = isNan(row['error']) and isNan(row['fwd_error']) and isNan(
        #     row['rev_error'])
        passed = row['error'] is None and row['fwd_error'] is None and row[
            'rev_error'] is None
        # If sample not in data yet
        if row['sample'] not in data:
            data[sample] = {
                'sample': sample,
                'run': row['run_name'],
                'conseq_passed': False,
                'contig_passed': False,
                'reference': None,
                'seqlen': None,
                'seqtype': None,
                'sequence': None,
                'failed': [],
                'error': None
            }
            if passed:
                data[sample]['reference'] = row['reference']
                data[sample]['seqlen'] = row['seqlen']
                data[sample]['conseq_passed'] = True
                data[sample]['sequence'] = row['sequence']
                data[sample]['seqtype'] = row['seqtype']
        # Else if sample is already in data
        else:
            # If we have already seen a conseq that passed, set them both to fail since this means multiple passing conseqs
            if passed and data[sample]['conseq_passed']:
                logger.critical('Sample "%s" already has a passed sequence!' %
                                sample)
                data[sample]['conseq_passed'] = False
                data[sample]['contig_passed'] = False
                data[sample]['error'] = 'Sample has multiple passed sequences'
                data[sample]['sequence'] = None
                data[sample]['seqtype'] = None
            # If we have not seen a conseq that passed, set to pass
            elif passed:
                data[sample]['reference'] = row['reference']
                data[sample]['seqlen'] = row['seqlen']
                data[sample]['conseq_passed'] = True
                data[sample]['sequence'] = row['sequence']
                data[sample]['seqtype'] = row['seqtype']
            # If not passed
            else:
                ## Determine type of error for certain cases
                # If conseq is not max, do not record it
                if ((row['error'] == 'contig not MAX')
                        or (row['error'] == 'is V3 sequence')):
                    continue
                elif not is_proviral(sample):
                    row['error'] = 'Sample is non-proviral'
                elif row['reference'] is None:
                    pass
                elif any(
                    [x in row['reference'] for x in ('reverse', 'unknown')]):
                    row['error'] = 'Sample does not align to HIV'
                nfailed = len(data[sample]['failed'])
                if nfailed > max_failed:
                    max_failed = nfailed
                data[sample]['failed'].append({
                    f'fail_error_{nfailed}':
                    row['error'],
                    f'fail_fwd_err_{nfailed}':
                    row['fwd_error'],
                    f'fail_rev_err_{nfailed}':
                    row['rev_error'],
                    f'fail_seqtype_{nfailed}':
                    seqtype,
                    f'fail_seqlen_{nfailed}':
                    row['seqlen'],
                    f'fail_sequence_{nfailed}':
                    row['sequence'],
                    f'fail_ref_{nfailed}':
                    row['reference']
                })

    # Go through all the contigs
    for index, row in contigs_df.iterrows():
        seqtype = 'contig'
        sample = row['sample']
        # passed = isNan(row['error']) and isNan(row['fwd_error']) and isNan(
        #     row['rev_error'])
        passed = row['error'] is None and row['fwd_error'] is None and row[
            'rev_error'] is None
        # If sample not in data yet, print a warning because we already went through all of the conseqs and so we should have captured every sample
        if row['sample'] not in data:
            logger.warning(
                'Sample "%s" not found in conseqs but was in contigs?!' %
                sample)
        # Else if sample is already in data
        else:
            if passed and data[sample]['contig_passed']:
                logger.critical('Sample "%s" already has a passed sequence!' %
                                sample)
                data[sample]['contig_passed'] = False
                data[sample]['conseq_passed'] = False
                data[sample]['error'] = 'Sample has multiple passed sequences'
                data[sample]['sequence'] = None
                data[sample]['seqtype'] = None
            # If the contig passes
            elif passed:
                data[sample]['contig_passed'] = True
                data[sample]['reference'] = row['reference']
                data[sample]['seqlen'] = row['seqlen']
                data[sample]['sequence'] = row['sequence']
                data[sample]['seqtype'] = row['seqtype']
            else:
                if not is_proviral(sample):
                    row['error'] = 'Sample is non-proviral'
                elif row['reference'] is None:
                    pass
                elif any(
                    [x in row['reference'] for x in ('reverse', 'unknown')]):
                    row['error'] = 'Sample does not align to HIV'
                nfailed = len(data[sample]['failed'])
                if nfailed > max_failed:
                    max_failed = nfailed
                data[sample]['failed'].append({
                    f'fail_error_{nfailed}':
                    row['error'],
                    f'fail_fwd_err_{nfailed}':
                    row['fwd_error'],
                    f'fail_rev_err_{nfailed}':
                    row['rev_error'],
                    f'fail_seqtype_{nfailed}':
                    seqtype,
                    f'fail_seqlen_{nfailed}':
                    row['seqlen'],
                    f'fail_sequence_{nfailed}':
                    row['sequence'],
                    f'fail_ref_{nfailed}':
                    row['reference']
                })
    outfile = outpath / 'outcome_summary.csv'

    fieldnames = [
        'sample', 'run', 'passed', 'reference', 'seqtype', 'seqlen', 'sequence'
    ]
    for i in range(max_failed):
        fieldnames += [
            f'fail_error_{i}', f'fail_fwd_err_{i}', f'fail_rev_err_{i}',
            f'fail_seqtype_{i}', f'fail_seqlen_{i}', f'fail_sequence_{i}',
            f'fail_ref_{i}'
        ]

    # Write the rows
    with open(outfile, 'w', newline='') as csvfile:
        writer = DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for sample in data:
            data[sample]['passed'] = data[sample]['conseq_passed'] or data[
                sample]['contig_passed']

            # Zabrina's request to remove all failed seqs if sample passed
            if data[sample]['passed']:
                data[sample]['failed'] = []

            # Zabrina's request to simply display a single error if all failures are due to not aligning to HIV
            count_non_hiv = 0
            is_hiv_indicies = []
            for i, fail in enumerate(data[sample]['failed']):
                if fail[f'fail_error_{i}'] == 'Sample does not align to HIV':
                    count_non_hiv += 1
                else:
                    is_hiv_indicies.append(i)
            # If the number of failures is equal to the count of non-hiv failures then all failures are due to non-hiv
            if count_non_hiv == len(
                    data[sample]['failed']) and not data[sample]['passed']:
                data[sample]['error'] = 'Sample does not align to HIV'
                data[sample]['failed'] = []
            # Otherwise at least one sequence was not non-hiv so we should display only the error for that sequence
            else:
                new_failed = []
                for i in is_hiv_indicies:
                    new_failed.append(data[sample]['failed'][i])
                data[sample]['failed'] = new_failed
            for fail in data[sample]['failed']:
                for k, v in fail.items():
                    data[sample][k] = v
            data[sample] = {
                k: v
                for k, v in data[sample].items() if k in fieldnames
            }
            writer.writerow(data[sample])
    return outfile


def getSamplesFromCascade(cascade_csv):
    all_samples = {}
    reader = DictReader(cascade_csv)
    for row in reader:
        all_samples[row['sample']] = int(row['remap'])
    return all_samples


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