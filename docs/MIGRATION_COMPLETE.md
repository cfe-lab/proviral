# Migration to mappy - Completion Report

## Status: ✅ COMPLETE

The migration from the minimap2 executable to the mappy Python binding has been successfully completed.

## Summary

- **Migration Date**: 2024
- **Total Tests**: 194 (all passing)
- **Migration Tests**: 59 (all passing)
- **Existing Tests**: 135 (all passing)

## Changes Made

### 1. New Functions in `cfeproviral/utils.py`

#### `mappy_available()` (Line ~361)
- Checks if the mappy module is available
- Returns `True` if mappy can be imported, `False` otherwise

#### `mappy_hit_to_sam_line()` (Line ~383)
- Converts a mappy Hit object to a SAM format line
- Handles:
  - FLAG calculation (reverse strand detection)
  - CIGAR string conversion from mappy format to SAM format
  - Sequence reverse complementation for reverse strand alignments
  - All 11 mandatory SAM fields

#### `align_with_mappy()` (Line ~423)
- Main alignment function using mappy
- Creates the same directory structure as the original function
- Writes:
  - `query.fasta` - query sequence
  - `target.fasta` - reference sequence
  - `alignment.sam` - SAM format alignment
  - `minimap2.log` - log file with alignment info
- Returns: Path to `alignment.sam` or `False` on failure

### 2. Modified Functions

#### `align()` (Line ~542)
- Now tries mappy first (if available)
- Falls back to minimap2 executable if mappy is unavailable
- Maintains full backward compatibility
- Updated docstring to reflect new behavior

## Key Technical Details

### Alignment Parameters

**Critical Discovery**: The original code uses `minimap2 -a` with **NO preset flag** (`-x`).

- Original: `minimap2 -a target.fasta query.fasta`
- Mappy equivalent: `mappy.Aligner(seq=target_seq)` (no preset parameter)

This was confirmed in test `test_exact_current_invocation_with_no_preset` which verified identical results:
- Position: 1001
- CIGAR: 5000M
- MAPQ: 60

### SAM Format Conversion

The conversion from mappy Hit objects to SAM format handles:

1. **CIGAR Operations**: Conversion from mappy tuples `[(length, op), ...]` to SAM strings
   - mappy op 0 → M (match/mismatch)
   - mappy op 1 → I (insertion)
   - mappy op 2 → D (deletion)
   - mappy op 4 → S (soft clipping)
   - mappy op 7 → = (match)
   - mappy op 8 → X (mismatch)

2. **FLAG Field**: Calculated from hit.strand
   - 0x10 for reverse strand
   - 0x4 for unmapped reads

3. **Coordinates**: 1-based SAM position from 0-based mappy position

4. **Sequence**: Reverse complemented for reverse strand alignments

### Directory Structure

The output structure is identical to the original implementation:

```
outdir/
└── query_name/
    ├── query.fasta
    ├── target.fasta
    ├── alignment.sam
    └── minimap2.log
```

## Benefits of Migration

1. **No External Dependencies**: No need for minimap2 executable to be installed
2. **Better Integration**: Native Python objects instead of subprocess calls
3. **Cleaner Error Handling**: Python exceptions instead of parsing stderr
4. **Performance**: Potentially faster for multiple alignments (no process spawning)
5. **Backward Compatibility**: Falls back to minimap2 if mappy unavailable

## Testing

### Migration-Specific Tests (59 tests)

1. **test_minimap2_compatibility.py** (22 tests)
   - Baseline minimap2 behavior
   - Basic mappy functionality
   - Equivalence validation (including critical no-preset test)
   - Edge cases
   - Performance validation
   - SAM format compatibility

2. **test_sam_conversion.py** (25 tests)
   - CIGAR string conversion
   - SAM flag calculation
   - SAM line construction
   - SAM header generation
   - Reverse complement utilities

3. **test_integration_real_data.py** (12 tests)
   - Real HIV sequence alignments
   - Current align() function validation
   - Gene splicer integration
   - Backward compatibility

### Existing Tests (135 tests)

All existing tests pass without modification, confirming that the migration maintains full backward compatibility.

## Validation Results

### Test Results
```
===========================
== 194 passed in 22.57s ===
===========================
```

### Key Test Validations

1. **Exact Behavior Replication**: Confirmed identical alignment results for current invocation
2. **Real Data**: Successfully aligned real HIV sequences (HXB2, GAG)
3. **Gene Splicer**: Integration with gene_splicer workflow validated
4. **Edge Cases**: Handled empty, short, long, mismatched sequences
5. **SAM Format**: Proper SAM format with valid headers and alignment lines

## Known Limitations

1. **Quality Scores**: mappy doesn't provide quality scores, so `*` is used in the QUAL field
2. **Secondary Alignments**: Currently only primary alignments are reported
3. **Paired-End**: Only single-end alignments supported (matching original behavior)

## Future Enhancements (Optional)

1. Add support for mappy preset parameters (e.g., 'sr', 'map-ont') via optional parameter
2. Add support for secondary alignments
3. Add progress reporting for large sequence alignments
4. Consider using aligntools library for more sophisticated CIGAR handling

## Migration Checklist Status

- ✅ Analyze codebase for minimap2 usage
- ✅ Create comprehensive test suite
- ✅ Implement mappy-based functions
- ✅ Update align() function
- ✅ Test with real data
- ✅ Validate all tests pass
- ✅ Document changes

## Deployment Recommendations

1. Ensure `mappy>=2.24` is in `pyproject.toml` dependencies (already present)
2. Update documentation to mention mappy as the primary alignment method
3. Keep minimap2 as optional fallback for systems where mappy cannot be installed
4. Consider updating CI/CD to test both mappy and minimap2 paths

## Conclusion

The migration has been successfully completed with:
- ✅ 100% test pass rate (194/194)
- ✅ Full backward compatibility
- ✅ Identical alignment results
- ✅ Clean code architecture
- ✅ Comprehensive documentation

The codebase now uses mappy as the primary alignment method while maintaining the ability to fall back to the minimap2 executable if needed.
