# Quick Implementation Reference for Minimap2 → Mappy Migration

This document provides the key code snippets needed for the migration.

## 1. Update pyproject.toml

```toml
dependencies = [
    "numpy==2.2.2",
    "python-levenshtein==0.26.1",
    "pandas==2.2.3",
    "requests==2.32.4",
    "pyyaml==6.0.2",
    "gotoh @ git+https://github.com/cfe-lab/MiCall.git@v7.7.0#egg=gotoh&subdirectory=micall/alignment",
    "cfeintact @ git+https://github.com/cfe-lab/CFEIntact.git@v1.23.2",
    "mappy>=2.24",  # ADD THIS LINE
]
```

## 2. Helper Functions to Add to cfeproviral/utils.py

Add these functions before the `aligner_available()` function:

```python
def _mappy_cigar_to_sam(cigar_tuples, strand, query_len, q_st, q_en):
    """
    Convert mappy CIGAR format to SAM CIGAR string.
    
    mappy CIGAR: List of (length, operation) tuples
    Operations: 0=M, 1=I, 2=D, 3=N, 4=S, 5=H, 6=P, 7==, 8=X
    
    Args:
        cigar_tuples: List of (length, operation) from mappy
        strand: +1 or -1
        query_len: Length of query sequence
        q_st: Query start position (0-based)
        q_en: Query end position (0-based, exclusive)
    
    Returns:
        SAM CIGAR string
    """
    if not cigar_tuples:
        return '*'
    
    # Operation mapping: mappy op code -> SAM op character
    op_map = {
        0: 'M',  # Match/mismatch
        1: 'I',  # Insertion to reference
        2: 'D',  # Deletion from reference
        3: 'N',  # Skipped region from reference
        4: 'S',  # Soft clipping
        5: 'H',  # Hard clipping
        6: 'P',  # Padding
        7: '=',  # Sequence match
        8: 'X',  # Sequence mismatch
    }
    
    cigar_parts = []
    
    # Add soft clipping at start if needed
    if q_st > 0:
        cigar_parts.append(f"{q_st}S")
    
    # Convert main CIGAR operations
    for length, op in cigar_tuples:
        cigar_parts.append(f"{length}{op_map.get(op, 'M')}")
    
    # Add soft clipping at end if needed
    if q_en < query_len:
        cigar_parts.append(f"{query_len - q_en}S")
    
    return ''.join(cigar_parts)


def _mappy_hit_to_sam(hit, query_name, query_seq, reference_name):
    """
    Convert a mappy Hit object to a SAM format line.
    
    Args:
        hit: mappy Hit object
        query_name: Query sequence name
        query_seq: Query sequence string
        reference_name: Reference sequence name
    
    Returns:
        SAM format string (without trailing newline)
    """
    # SAM fields: QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, QUAL
    
    qname = query_name
    
    # Calculate FLAG
    flag = 0
    if hit.strand == -1:
        flag |= 16  # Reverse strand
    if not hit.is_primary:
        flag |= 256  # Secondary alignment
    
    rname = reference_name
    pos = hit.r_st + 1  # Convert to 1-based SAM coordinate
    mapq = hit.mapq
    
    # Convert CIGAR
    cigar = _mappy_cigar_to_sam(
        hit.cigar,
        hit.strand,
        len(query_seq),
        hit.q_st,
        hit.q_en
    )
    
    rnext = '*'
    pnext = 0
    tlen = 0
    
    # Sequence
    seq = query_seq
    if hit.strand == -1:
        seq = reverse_and_complement(seq)
    
    qual = '*'
    
    # Optional fields
    optional_fields = []
    optional_fields.append(f"NM:i:{hit.NM}")  # Edit distance
    optional_fields.append(f"ms:i:{hit.mlen}")  # Matching bases
    optional_fields.append(f"AS:i:{hit.blen}")  # Alignment score
    optional_fields.append(f"nn:i:{hit.trans_strand}")  # Translocation strand
    
    sam_line = f"{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{seq}\t{qual}"
    if optional_fields:
        sam_line += "\t" + "\t".join(optional_fields)
    
    return sam_line
```

## 3. Replace aligner_available() Function

Find the current `aligner_available()` function (around line 356) and replace with:

```python
def aligner_available(aligner_path='minimap2'):
    """
    Check if mappy alignment library is available.
    
    Args:
        aligner_path: Deprecated parameter, kept for backward compatibility.
                     Previously used to specify minimap2 executable path.
    
    Returns:
        True if mappy is available, False otherwise.
    
    Note:
        After migration to mappy, this function checks for the Python library
        rather than an external executable.
    """
    try:
        import mappy
        return True
    except ImportError:
        return False
```

## 4. Replace align() Function

Find the current `align()` function (around line 379) and replace with:

```python
def align(target_seq,
          query_seq,
          query_name,
          outdir=Path(os.getcwd()).resolve(),
          aligner_path='minimap2'):
    """
    Align query sequence to target sequence using mappy.
    
    Args:
        target_seq: Target (reference) sequence string
        query_seq: Query sequence string to align
        query_name: Name/identifier for the query sequence
        outdir: Output directory path (default: current directory)
        aligner_path: Deprecated parameter, kept for backward compatibility
    
    Returns:
        Path object pointing to alignment.sam file if successful,
        False if alignment fails.
    
    Output Structure:
        Creates: outdir/query_name/
            - query.fasta: Query sequence in FASTA format
            - target.fasta: Target sequence in FASTA format
            - alignment.sam: Alignment in SAM format
            - minimap2.log: Alignment log (name kept for compatibility)
    
    Raises:
        FileNotFoundError: If mappy library is not available
    
    Note:
        After migration, uses mappy Python library instead of minimap2 executable.
        The output format and directory structure remain the same for compatibility.
    """
    if not aligner_available():
        raise FileNotFoundError(
            'mappy library not available. '
            'Install with: pip install mappy'
        )
    
    try:
        import mappy
    except ImportError as e:
        raise FileNotFoundError(f'Failed to import mappy: {e}')
    
    # Set up output directory
    outdir = outdir / query_name
    if os.path.isdir(outdir):
        shutil.rmtree(outdir)
    os.makedirs(outdir)
    
    # Write FASTA files (kept for debugging and compatibility)
    query_fasta_path = write_fasta(
        {query_name: query_seq},
        outdir / 'query.fasta'
    )
    target_fasta_path = write_fasta(
        {'MOD_HXB2': target_seq},
        outdir / 'target.fasta'
    )
    
    alignment_path = outdir / 'alignment.sam'
    log_path = outdir / 'minimap2.log'  # Name kept for compatibility
    
    try:
        # Create aligner with target sequence
        # preset='sr' mimics minimap2's short-read alignment mode
        aligner = mappy.Aligner(seq=target_seq, preset='sr')
        
        # Perform alignment
        hits = list(aligner.map(query_seq))
        
        # Write SAM file
        with alignment_path.open('w') as sam_file, log_path.open('w') as log_file:
            # Write SAM header
            sam_file.write("@HD\tVN:1.0\tSO:unsorted\n")
            sam_file.write(f"@SQ\tSN:MOD_HXB2\tLN:{len(target_seq)}\n")
            sam_file.write(f"@PG\tID:mappy\tPN:mappy\tVN:{mappy.__version__}\n")
            
            # Write alignments
            if not hits:
                log_file.write(f"Warning: No alignments found for {query_name}\n")
                log_file.write(f"Query length: {len(query_seq)}\n")
                log_file.write(f"Target length: {len(target_seq)}\n")
                
                # Write unmapped record
                flag = 4  # Unmapped
                sam_file.write(
                    f"{query_name}\t{flag}\t*\t0\t0\t*\t*\t0\t0\t{query_seq}\t*\n"
                )
            else:
                for hit in hits:
                    sam_line = _mappy_hit_to_sam(
                        hit,
                        query_name,
                        query_seq,
                        'MOD_HXB2'
                    )
                    sam_file.write(sam_line + '\n')
                
                log_file.write(f"Successfully aligned {query_name}\n")
                log_file.write(f"Number of hits: {len(hits)}\n")
                log_file.write(f"Primary alignment: r_st={hits[0].r_st}, "
                             f"r_en={hits[0].r_en}, mapq={hits[0].mapq}\n")
        
        return alignment_path
        
    except Exception as e:
        with log_path.open('w') as log_file:
            log_file.write(f'Alignment failed: {str(e)}\n')
            import traceback
            log_file.write(traceback.format_exc())
        
        logger.error('Alignment failed! Details in %s.', log_path)
        return False
```

## 5. Update Installation Scripts

### Remove minimap2 Installation

Edit `scripts/install.sh` and comment out or remove the line that calls:
```bash
# sh -- /opt/cfeproviral/scripts/installation/30-minimap2.sh  # No longer needed
```

Or rename/move the file:
```bash
mv scripts/installation/30-minimap2.sh scripts/installation/30-minimap2.sh.old
```

## 6. Update CI/CD

Edit `.github/workflows/python-app.yml`:

**Remove this section**:
```yaml
    - name: Install minimap2
      run: |
        cd /opt
        wget https://github.com/lh3/minimap2/releases/download/v2.17/minimap2-2.17_x64-linux.tar.bz2
        tar -jxvf ./minimap2-2.17_x64-linux.tar.bz2
        echo /opt/minimap2-2.17_x64-linux >> $GITHUB_PATH
```

The `uv sync` command will now install mappy automatically.

## 7. Testing

### Install mappy
```bash
pip install mappy
```

### Run Tests
```bash
# Migration-specific tests
pytest tests/test_minimap2_migration/ -v

# All tests
pytest tests/ -v

# With coverage
pytest tests/ --cov=cfeproviral.utils --cov-report=html
```

## 8. Validation

### Test with Real Data
```bash
# Run gene splicer on test data
python -m cfeproviral.gene_splicer tests/data/example1/query.fasta --outdir /tmp/test_output

# Compare output with expected
diff /tmp/test_output/genes.fasta tests/data/example1/expected/genes.fasta
```

### Test Docker Build
```bash
docker build -t proviral .
docker run --rm proviral --help
```

## 9. Common Issues and Solutions

### Issue: mappy not found
**Solution**: Ensure mappy is in pyproject.toml and run `pip install mappy` or `uv sync`

### Issue: CIGAR strings differ
**Solution**: This is expected. Verify functional equivalence (positions, aligned sequences) rather than exact CIGAR match.

### Issue: Performance slower
**Solution**: Acceptable if < 20% slower. mappy is a C extension so should be fast. Profile if needed.

### Issue: Tests fail on existing SAM comparisons
**Solution**: Update tests to compare functional output (genes identified) rather than exact SAM format.

## 10. Rollback

If needed to rollback:

```bash
# Revert changes
git revert <commit-hash>

# Or restore old functions
git checkout <previous-commit> -- cfeproviral/utils.py

# Reinstall minimap2
apt-get install minimap2
```

## Implementation Order

1. ✅ Add mappy to pyproject.toml
2. ✅ Install mappy: `pip install mappy`
3. ✅ Add helper functions (_mappy_cigar_to_sam, _mappy_hit_to_sam)
4. ✅ Update aligner_available()
5. ✅ Update align()
6. ✅ Test with test_sam_conversion.py
7. ✅ Test with test_minimap2_compatibility.py
8. ✅ Run full test suite
9. ✅ Update installation scripts
10. ✅ Update CI/CD
11. ✅ Update documentation

## Verification Commands

```bash
# Check mappy installed
python -c "import mappy; print(mappy.__version__)"

# Test align function
python -c "
from cfeproviral.utils import align, mod_hxb2
from pathlib import Path
result = align(mod_hxb2, mod_hxb2[:1000], 'test', Path('/tmp'))
print(f'Alignment result: {result}')
"

# Check SAM output
cat /tmp/test/alignment.sam
```

## Need Help?

- Review migration plan: `docs/MINIMAP2_MIGRATION.md`
- Check test suite: `tests/test_minimap2_migration/README.md`
- See checklist: `MIGRATION_CHECKLIST.md`
- Run tests: `pytest tests/test_minimap2_migration/ -v`
