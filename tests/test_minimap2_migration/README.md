# Minimap2 Migration Test Suite

This directory contains a comprehensive test suite to ensure safe migration from minimap2 executable to mappy (Python binding).

## Test Files

### 1. `test_minimap2_compatibility.py`
**Purpose**: Validate that minimap2 executable and mappy produce equivalent results.

**Test Classes**:
- `TestMinimap2Baseline`: Capture current minimap2 behavior
  - Verify minimap2 executable is available
  - Test basic alignment
  - Validate SAM output format
  
- `TestMappyCompatibility`: Test mappy library functionality
  - Test basic import and API
  - Test alignment with sequences and files
  - Verify CIGAR format
  
- `TestMinimap2ToMappyEquivalence`: Compare outputs
  - CIGAR string comparison
  - Alignment position comparison
  - Ensure functional equivalence
  
- `TestCurrentCodebaseCompatibility`: Verify existing code
  - Test `utils.align()` function signature
  - Test `aligner_available()` function
  - Document expected behavior
  
- `TestEdgeCases`: Handle special cases
  - Empty sequences
  - Very short sequences
  - Mismatched sequences
  - Sequences with N (unknown nucleotides)
  - Long sequences (HIV genome ~9kb)
  
- `TestPerformance`: Ensure acceptable speed
- `TestSAMFormatCompatibility`: Verify SAM format conversion

### 2. `test_sam_conversion.py`
**Purpose**: Test SAM format conversion utilities.

**Test Classes**:
- `TestCigarConversion`: CIGAR string conversion
  - Match operations (M)
  - Insertions (I) and deletions (D)
  - Soft clipping (S)
  - Complex CIGAR strings
  - Edge cases
  
- `TestSamFlagCalculation`: SAM FLAG field
  - Forward/reverse strand
  - Primary/secondary alignments
  - Unmapped reads
  
- `TestSamLineConstruction`: Complete SAM line
  - All required fields
  - Optional fields
  - Reverse strand handling
  - Unmapped reads
  
- `TestSamHeaderGeneration`: SAM header
  - Minimal header
  - Program information
  - Format compliance
  
- `TestReverseComplement`: DNA sequence operations

### 3. `test_integration_real_data.py`
**Purpose**: Integration tests with actual HIV sequences from the test suite.

**Test Classes**:
- `TestWithRealHIVSequences`: Use real test data
  - Load HXB2 reference genome
  - Test with existing SAM files
  - Align GAG gene
  - Compare with existing outputs
  
- `TestAlignFunctionWithRealData`: Test `align()` function
  - Verify function signature
  - Test output compatibility
  - Ensure directory structure
  
- `TestGeneSplicerIntegration`: Full workflow
  - Align → Load SAM → Splice genes → Extract sequences
  - Test complete gene identification pipeline
  
- `TestBackwardCompatibility`: Ensure no breaking changes
  - Parameter compatibility
  - Return value compatibility
  - Directory structure compatibility

## Running Tests

### Prerequisites

```bash
# Current system (before migration)
apt-get install minimap2

# After migration
pip install mappy
```

### Run All Tests

```bash
# From project root
pytest tests/test_minimap2_migration/ -v

# Run specific test file
pytest tests/test_minimap2_migration/test_minimap2_compatibility.py -v

# Run specific test class
pytest tests/test_minimap2_migration/test_sam_conversion.py::TestCigarConversion -v

# Run with coverage
pytest tests/test_minimap2_migration/ --cov=cfeproviral.utils --cov-report=html
```

### Expected Behavior

**Before Migration** (with minimap2 executable):
- `TestMinimap2Baseline` tests should PASS
- `TestMappyCompatibility` tests will SKIP (mappy not installed)
- `TestCurrentCodebaseCompatibility` tests should PASS

**After Migration** (with mappy):
- All tests should PASS
- No tests should SKIP (except on missing test data)
- Performance tests should complete within acceptable time

## Test Strategy

### Phase 1: Baseline (Current)
1. Run `TestMinimap2Baseline` to capture current behavior
2. Save outputs for comparison
3. Document any quirks or edge cases

### Phase 2: Migration Development
1. Install mappy: `pip install mappy`
2. Implement new `align()` and `aligner_available()` functions
3. Run `TestMappyCompatibility` to verify mappy usage
4. Run `test_sam_conversion.py` to verify conversion logic

### Phase 3: Validation
1. Run `TestMinimap2ToMappyEquivalence` to compare outputs
2. Run `test_integration_real_data.py` with real sequences
3. Run full test suite: `pytest tests/`
4. Compare results with baseline

### Phase 4: Regression Testing
1. Run all existing tests
2. Verify no behavioral changes
3. Test on multiple platforms (Linux, macOS, Windows)

## Success Criteria

✅ All tests pass
✅ SAM output is functionally equivalent
✅ Gene splicing produces same results
✅ Performance is acceptable (< 20% slower)
✅ Works on Windows (currently doesn't)
✅ No external dependencies required

## Known Issues

### mappy vs minimap2 Differences

1. **CIGAR Format**: mappy returns tuples `(length, op)`, minimap2 outputs strings
   - **Solution**: Conversion function `_mappy_cigar_to_sam()`
   
2. **Header Format**: May differ slightly
   - **Solution**: Generate compliant SAM headers manually
   
3. **Alignment Parameters**: Need to match minimap2's default behavior
   - **Solution**: Use `preset='sr'` (short reads) as default

4. **Optional Fields**: mappy may include different tags
   - **Solution**: Include essential tags only

## Maintenance

### Adding New Tests

When adding tests, follow this pattern:

```python
class TestNewFeature:
    """Description of what is being tested."""
    
    @pytest.fixture
    def mappy_module(self):
        """Import mappy if available."""
        try:
            import mappy
            return mappy
        except ImportError:
            pytest.skip("mappy not installed")
    
    def test_specific_behavior(self, mappy_module):
        """Test description."""
        # Arrange
        target_seq = "ATCG" * 100
        query_seq = "ATCG" * 50
        
        # Act
        aligner = mappy_module.Aligner(seq=target_seq)
        hits = list(aligner.map(query_seq))
        
        # Assert
        assert len(hits) > 0, "Should produce alignments"
```

### Updating Tests

When the code changes:
1. Update test expectations
2. Add regression tests for bugs
3. Document behavior changes
4. Update this README

## References

- [mappy documentation](https://github.com/lh3/minimap2/tree/master/python)
- [SAM format specification](https://samtools.github.io/hts-specs/SAMv1.pdf)
- [minimap2 paper](https://doi.org/10.1093/bioinformatics/bty191)
- [Migration plan](../docs/MINIMAP2_MIGRATION.md)

## Contact

For questions about these tests:
- See migration plan: `docs/MINIMAP2_MIGRATION.md`
- Check implementation: `cfeproviral/utils.py`
- Review test output: `pytest -v --tb=long`
