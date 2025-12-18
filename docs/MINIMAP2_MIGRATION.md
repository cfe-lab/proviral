# Migration Plan: minimap2 Executable to mappy (Python Binding)

## Executive Summary

This document outlines the plan to migrate from using minimap2 as an external executable to using mappy, the official Python binding for minimap2. This migration will simplify deployment, eliminate external dependencies, and provide better error handling.

## Current State Analysis

### Current Implementation

**Primary Location**: `cfeproviral/utils.py`

**Key Functions**:

1. `aligner_available(aligner_path='minimap2')` (line 356)
   - Checks if minimap2 executable is available
   - Uses `subprocess.run()` to test with `-h` flag
   - Returns boolean

2. `align(target_seq, query_seq, query_name, outdir, aligner_path='minimap2')` (line 379)
   - Main alignment function
   - Writes FASTA files for target and query
   - Calls minimap2 via subprocess: `minimap2 -a target.fasta query.fasta`
   - Outputs SAM format to `alignment.sam`
   - Creates log file at `minimap2.log`
   - Returns Path to SAM file or False on failure

**Usage Pattern**:
- Called from `cfeproviral/gene_splicer.py` (line 29)
- Used in HIV genome alignment for gene splicing
- Results parsed by `utils.load_samfile()`

### Current Dependencies

1. **Installation**: `scripts/installation/30-minimap2.sh`
   - System package: `apt-get install -y minimap2`

2. **Docker**: `Dockerfile`
   - Installs via install.sh script

3. **CI/CD**: `.github/workflows/python-app.yml`
   - Downloads specific version (v2.17) from GitHub releases
   - Manual installation in CI environment

4. **Test Files**:
   - `tests/test_gene_splicer/test_gene_splicer.py` - skips tests on Windows
   - Pre-generated SAM files for testing

### Current Command

```bash
minimap2 -a <target.fasta> <query.fasta> > alignment.sam 2> minimap2.log
```

**Flags**:
- `-a`: Output in SAM format (instead of PAF)

## Proposed Changes

### 1. Replace subprocess calls with mappy

**New Dependencies**:
```toml
# Add to pyproject.toml
dependencies = [
    ...
    "mappy>=2.24",  # Python bindings for minimap2
]
```

### 2. Update `cfeproviral/utils.py`

#### New Implementation of `aligner_available()`

```python
def aligner_available(aligner_path='minimap2'):
    """Check if mappy is available (aligner_path parameter kept for backward compatibility)."""
    try:
        import mappy
        return True
    except ImportError:
        return False
```

#### New Implementation of `align()`

```python
def align(target_seq,
          query_seq,
          query_name,
          outdir=Path(os.getcwd()).resolve(),
          aligner_path='minimap2'):  # Kept for backward compatibility
    """
    Align query sequence to target sequence using mappy.
    
    Args:
        target_seq: Target sequence string
        query_seq: Query sequence string
        query_name: Name of the query sequence
        outdir: Output directory path
        aligner_path: Deprecated, kept for backward compatibility
    
    Returns:
        Path to alignment.sam file or False on failure
    """
    if not aligner_available():
        raise FileNotFoundError('mappy library not available. Install with: pip install mappy')
    
    import mappy
    
    outdir = outdir / query_name
    if os.path.isdir(outdir):
        shutil.rmtree(outdir)
    os.makedirs(outdir)
    
    # Still write FASTA files for debugging/logging purposes
    query_fasta_path = write_fasta({query_name: query_seq},
                                   outdir / 'query.fasta')
    target_fasta_path = write_fasta({'MOD_HXB2': target_seq},
                                    outdir / 'target.fasta')
    
    alignment_path = outdir / 'alignment.sam'
    log_path = outdir / 'minimap2.log'
    
    try:
        # Create aligner with target sequence
        aligner = mappy.Aligner(seq=target_seq, preset='sr')
        
        # Perform alignment
        hits = list(aligner.map(query_seq))
        
        # Convert to SAM format
        with alignment_path.open('w') as sam_file, log_path.open('w') as log_file:
            # Write SAM header
            sam_file.write("@HD\tVN:1.0\tSO:unsorted\n")
            sam_file.write(f"@SQ\tSN:MOD_HXB2\tLN:{len(target_seq)}\n")
            sam_file.write(f"@PG\tID:mappy\tPN:mappy\tVN:{mappy.__version__}\n")
            
            # Write alignments
            if not hits:
                log_file.write(f"Warning: No alignments found for {query_name}\n")
                # Write unmapped record
                flag = 4  # Unmapped
                sam_file.write(f"{query_name}\t{flag}\t*\t0\t0\t*\t*\t0\t0\t{query_seq}\t*\n")
            else:
                for hit in hits:
                    sam_line = _mappy_hit_to_sam(hit, query_name, query_seq, 'MOD_HXB2')
                    sam_file.write(sam_line + '\n')
                
                log_file.write(f"Successfully aligned {query_name}: {len(hits)} hit(s)\n")
        
        return alignment_path
        
    except Exception as e:
        with log_path.open('w') as log_file:
            log_file.write(f'Alignment failed: {str(e)}\n')
        logger.error('Alignment failed! Details in %s.', log_path)
        return False


def _mappy_hit_to_sam(hit, query_name, query_seq, reference_name):
    """
    Convert a mappy Hit object to a SAM format line.
    
    Args:
        hit: mappy Hit object
        query_name: Query sequence name
        query_seq: Query sequence
        reference_name: Reference sequence name
    
    Returns:
        SAM format string (without newline)
    """
    # SAM fields: QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, QUAL
    
    qname = query_name
    
    # FLAG: Calculate SAM flag
    flag = 0
    if hit.strand == -1:
        flag |= 16  # Reverse strand
    if not hit.is_primary:
        flag |= 256  # Secondary alignment
    
    rname = reference_name
    pos = hit.r_st + 1  # Convert to 1-based
    mapq = hit.mapq
    
    # Convert CIGAR from mappy format to SAM string
    cigar = _mappy_cigar_to_sam(hit.cigar, hit.strand, len(query_seq), hit.q_st, hit.q_en)
    
    rnext = '*'
    pnext = 0
    tlen = 0
    
    # Sequence and quality
    seq = query_seq
    if hit.strand == -1:
        # Reverse complement for reverse strand
        seq = reverse_and_complement(seq)
    
    qual = '*'  # Quality scores not available
    
    # Build optional fields
    optional = f"NM:i:{hit.NM}\tms:i:{hit.mlen}\tAS:i:{hit.blen}\tnn:i:{hit.trans_strand}"
    
    return f"{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{seq}\t{qual}\t{optional}"


def _mappy_cigar_to_sam(cigar_tuples, strand, query_len, q_st, q_en):
    """
    Convert mappy CIGAR format to SAM CIGAR string.
    
    mappy CIGAR: List of (length, operation) tuples
    Operations: 0=M, 1=I, 2=D, 3=N, 4=S, 5=H, 6=P, 7==, 8=X
    
    Args:
        cigar_tuples: List of (length, operation) from mappy
        strand: +1 or -1
        query_len: Length of query sequence
        q_st: Query start position
        q_en: Query end position
    
    Returns:
        SAM CIGAR string
    """
    if not cigar_tuples:
        return '*'
    
    # Operation mapping
    op_map = {
        0: 'M',  # Match/mismatch
        1: 'I',  # Insertion
        2: 'D',  # Deletion
        3: 'N',  # Skipped region
        4: 'S',  # Soft clipping
        5: 'H',  # Hard clipping
        6: 'P',  # Padding
        7: '=',  # Sequence match
        8: 'X',  # Sequence mismatch
    }
    
    cigar_str = []
    
    # Add soft clipping at start if needed
    if q_st > 0:
        cigar_str.append(f"{q_st}S")
    
    # Convert CIGAR operations
    for length, op in cigar_tuples:
        cigar_str.append(f"{length}{op_map.get(op, 'M')}")
    
    # Add soft clipping at end if needed
    if q_en < query_len:
        cigar_str.append(f"{query_len - q_en}S")
    
    return ''.join(cigar_str)
```

### 3. Update Installation Scripts

**Delete or comment out**: `scripts/installation/30-minimap2.sh`

**Update**: `scripts/install.sh`
- Remove call to `30-minimap2.sh`
- mappy will be installed via pip/uv as a Python dependency

### 4. Update CI/CD

**`.github/workflows/python-app.yml`**:

Remove minimap2 installation step:
```yaml
# DELETE THIS SECTION
# - name: Install minimap2
#   run: |
#     cd /opt
#     wget https://github.com/lh3/minimap2/releases/download/v2.17/minimap2-2.17_x64-linux.tar.bz2
#     tar -jxvf ./minimap2-2.17_x64-linux.tar.bz2
#     echo /opt/minimap2-2.17_x64-linux >> $GITHUB_PATH
```

mappy will be installed automatically via `uv sync`

### 5. Update Tests

**`tests/test_gene_splicer/test_gene_splicer.py`**:

Remove the Windows skip comment:
```python
# REMOVE: # I have to skip this section because minimap2 does not run on windows
```

Tests should now work on all platforms including Windows!

### 6. Update Documentation

**Files to update**:
- `docs/introduction.md` - Update mention of minimap2
- `docs/interpretation.md` - Log files will still exist
- `docs/installation.md` - Simplify installation instructions
- `docs/steps.md` - Update description

**Key message**: 
- Installation is now simpler (no external dependencies)
- Works on all platforms including Windows
- Functionally equivalent behavior

## Benefits

### 1. Simplified Installation
- No need for system package managers
- No need for manual downloads in CI
- Single `pip install` or `uv sync` handles everything

### 2. Cross-Platform Compatibility
- Works on Windows (current implementation doesn't)
- Consistent behavior across Linux, macOS, Windows

### 3. Better Error Handling
- Python exceptions instead of subprocess errors
- Better debugging with Python stack traces
- No need to parse stderr logs

### 4. Performance
- Potential performance improvement (no subprocess overhead)
- Direct memory access (no file I/O for alignment)
- Can keep aligner object in memory for multiple queries

### 5. Maintainability
- Fewer external dependencies to manage
- No version mismatches between systems
- Easier to debug

## Risks and Mitigations

### Risk 1: Output Format Differences
**Risk**: mappy might produce slightly different SAM output than minimap2 executable

**Mitigation**: 
- Comprehensive test suite to verify equivalence
- SAM conversion function thoroughly tested
- Pre-generated expected outputs to validate against

### Risk 2: Performance Differences
**Risk**: Python binding might be slower than native executable

**Mitigation**:
- Performance tests included in test suite
- mappy is a C extension, should be comparable speed
- Can revert if performance is unacceptable

### Risk 3: Breaking Existing Workflows
**Risk**: External tools/scripts might depend on minimap2 executable

**Mitigation**:
- Keep backward compatible function signatures
- Maintain same output file structure
- Document migration clearly

### Risk 4: CIGAR String Conversion
**Risk**: Converting mappy CIGAR format to SAM might introduce bugs

**Mitigation**:
- Detailed unit tests for CIGAR conversion
- Test with real HIV sequences
- Compare outputs with existing SAM files

## Migration Steps

### Phase 1: Preparation (Current)
- [x] Analyze current usage
- [x] Create test suite
- [ ] Run baseline tests with current implementation
- [ ] Document expected behavior

### Phase 2: Implementation
- [ ] Add mappy to dependencies
- [ ] Implement new align() function
- [ ] Implement SAM conversion helpers
- [ ] Update aligner_available()
- [ ] Keep old code commented for reference

### Phase 3: Testing
- [ ] Run all existing tests
- [ ] Run new compatibility tests
- [ ] Compare outputs with pre-generated SAM files
- [ ] Performance testing
- [ ] Test on multiple platforms

### Phase 4: Integration
- [ ] Update installation scripts
- [ ] Update CI/CD pipeline
- [ ] Update documentation
- [ ] Remove old code

### Phase 5: Validation
- [ ] Run full pipeline on test data
- [ ] Compare results with previous runs
- [ ] Verify gene splicing outputs
- [ ] User acceptance testing

## Testing Strategy

### Unit Tests
- Test SAM format conversion
- Test CIGAR string conversion
- Test flag calculation
- Test error handling

### Integration Tests
- Test full alignment workflow
- Test with real HIV sequences
- Test with existing test data
- Compare with pre-generated outputs

### Compatibility Tests
- Side-by-side comparison with minimap2 executable
- Verify SAM format compliance
- Test edge cases (empty, short, mismatched sequences)

### Performance Tests
- Benchmark alignment speed
- Memory usage comparison
- Large sequence handling

## Rollback Plan

If migration fails:

1. **Immediate Rollback**:
   - Revert to previous commit
   - Restore minimap2 executable dependency
   - Re-enable installation scripts

2. **Partial Rollback**:
   - Keep mappy as optional dependency
   - Add fallback to minimap2 executable
   - User can choose via config

```python
# Fallback implementation
def align(..., use_mappy=True):
    if use_mappy and aligner_available_mappy():
        return align_with_mappy(...)
    else:
        return align_with_executable(...)
```

## Success Criteria

Migration is successful if:

1. ✅ All existing tests pass
2. ✅ New compatibility tests pass
3. ✅ Output SAM files are functionally equivalent
4. ✅ Performance is acceptable (within 20% of current)
5. ✅ Works on Windows, Linux, and macOS
6. ✅ Installation is simpler
7. ✅ No regression in gene splicing results
8. ✅ Documentation is updated

## Timeline

- **Phase 1 (Preparation)**: 1-2 days
- **Phase 2 (Implementation)**: 2-3 days
- **Phase 3 (Testing)**: 3-5 days
- **Phase 4 (Integration)**: 1-2 days
- **Phase 5 (Validation)**: 2-3 days

**Total Estimated Time**: 9-15 days

## References

- mappy documentation: https://github.com/lh3/minimap2/tree/master/python
- SAM format specification: https://samtools.github.io/hts-specs/SAMv1.pdf
- minimap2 paper: https://doi.org/10.1093/bioinformatics/bty191

## Contact

For questions about this migration:
- Review migration test suite: `tests/test_minimap2_migration/`
- Check implementation: `cfeproviral/utils.py`
- See documentation updates: `docs/MINIMAP2_MIGRATION.md`
