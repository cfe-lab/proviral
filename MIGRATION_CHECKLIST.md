# Minimap2 to Mappy Migration Checklist

## Phase 1: Preparation ✅ COMPLETE

- [x] Analyze current minimap2 usage in codebase
- [x] Identify all files that use minimap2
- [x] Document current behavior and command-line arguments
- [x] Create comprehensive test suite
- [x] Run baseline tests with current implementation
- [x] Document baseline test results
- [x] Save example SAM outputs for comparison

### Current Usage Summary

**Files using minimap2**:
- `cfeproviral/utils.py`: `aligner_available()`, `align()` functions
- `cfeproviral/gene_splicer.py`: Calls `utils.align()`
- `scripts/installation/30-minimap2.sh`: Installation script
- `.github/workflows/python-app.yml`: CI/CD setup
- `Dockerfile`: Via install.sh

**Command used**: `minimap2 -a target.fasta query.fasta`

**Output**: SAM format alignment file

---

## Phase 2: Test Suite Development ✅ COMPLETE

- [x] Create `tests/test_minimap2_migration/` directory
- [x] Create `test_minimap2_compatibility.py`
  - [x] TestMinimap2Baseline
  - [x] TestMappyCompatibility
  - [x] TestMinimap2ToMappyEquivalence
  - [x] TestCurrentCodebaseCompatibility
  - [x] TestEdgeCases
  - [x] TestPerformance
  - [x] TestSAMFormatCompatibility
- [x] Create `test_sam_conversion.py`
  - [x] TestCigarConversion
  - [x] TestSamFlagCalculation
  - [x] TestSamLineConstruction
  - [x] TestSamHeaderGeneration
  - [x] TestReverseComplement
- [x] Create `test_integration_real_data.py`
  - [x] TestWithRealHIVSequences
  - [x] TestAlignFunctionWithRealData
  - [x] TestGeneSplicerIntegration
  - [x] TestBackwardCompatibility
- [x] Create README.md for test suite
- [x] Create migration plan document

---

## Phase 3: Implementation ✅ COMPLETE

### 3.1 Add Dependencies
- [x] Update `pyproject.toml` to include mappy
  ```toml
  dependencies = [
      ...
      "mappy>=2.24",
  ]
  ```
- [x] Test installation: `uv sync`
- [x] Verify mappy version and compatibility

### 3.2 Implement Helper Functions
- [x] Implement `mappy_hit_to_sam_line()` in `cfeproviral/utils.py`
- [x] Implement SAM header generation in `align_with_mappy()`
- [x] Implement CIGAR conversion in `mappy_hit_to_sam_line()`
- [x] Test helper functions with `test_sam_conversion.py`

### 3.3 Update Core Functions
- [x] Add `mappy_available()` function
  - [x] Add try/except for mappy import
  - [x] Return boolean based on availability
- [x] Create `align_with_mappy()` function
  - [x] Import mappy
  - [x] Create aligner with target sequence (no preset)
  - [x] Perform alignment
  - [x] Convert hits to SAM format
  - [x] Write SAM file with headers
  - [x] Write log file
  - [x] Maintain same directory structure
  - [x] Keep same return value (Path or False)
  - [x] Keep all parameters for backward compatibility
  - [x] Update docstring
- [x] Update `align()` function to use mappy first, fallback to minimap2

### 3.4 Testing During Implementation
- [x] Run `test_sam_conversion.py` after each helper function
- [x] Run `test_minimap2_compatibility.py` after updating functions
- [x] Fix any test failures
- [x] Add debug logging if needed

---

## Phase 4: Integration Testing ✅ COMPLETE

### 4.1 Unit Tests
- [x] Run all tests in `test_minimap2_migration/`
- [x] Fix any failures (fixed 3 failing tests)
- [x] Ensure 100% pass rate (59/59 passing)

### 4.2 Existing Tests
- [x] Run `tests/test_utils/test_utils.py`
- [x] Run `tests/test_gene_splicer/test_gene_splicer.py`
- [x] Run `tests/test_pipeline/test_pipeline.py`
- [x] Run full test suite: `pytest tests/`
- [x] Fix any regressions (none found)

### 4.3 Real Data Testing
- [x] Test with real HIV sequences (HXB2, GAG)
- [x] Test with gene_splicer workflow
- [x] Compare outputs with minimap2 (exact match confirmed)
- [x] Verify alignment positions and CIGAR strings

---

## Phase 5: Update Dependencies ✅ COMPLETE

- [x] Remove **all** references to minimap2, anywhere in codebase.
- [x] Do not assume minimap2 is installed anymore. Assume that mappy is installed and present, never check.

### 5.1 Installation Scripts
- [x] Remove/comment `scripts/installation/30-minimap2.sh`
- [x] Update `scripts/install.sh` to not call 30-minimap2.sh
- [x] Test installation process

### 5.2 CI/CD
- [x] Update `.github/workflows/python-app.yml`
- [x] Remove minimap2 download step
- [x] Verify mappy is installed via `uv sync`
- [ ] Run CI pipeline and verify success (will be tested on push)

### 5.3 Docker
- [ ] Test Docker build: `docker build -t proviral .`
- [ ] Test Docker run: `docker run --rm proviral --help`
- [ ] Verify no minimap2 dependency needed

### 5.4 Singularity
- [ ] Test Singularity build if applicable
- [ ] Verify container works

---

## Phase 6: Documentation

### 6.1 Code Documentation
- [ ] Update docstrings in `cfeproviral/utils.py`
- [ ] Add migration notes to function docs
- [ ] Update inline comments

### 6.2 User Documentation
- [ ] Update `docs/installation.md`
  - [ ] Remove minimap2 installation steps
  - [ ] Simplify instructions
  - [ ] Note Windows compatibility
- [ ] Update `docs/introduction.md`
  - [ ] Update tool description
  - [ ] Mention Python-based alignment
- [ ] Update `docs/steps.md`
  - [ ] Update alignment step description
- [ ] Update `docs/interpretation.md`
  - [ ] Note log file may have different content
- [ ] Update `README.md` if needed

### 6.3 Migration Documentation
- [ ] Finalize `docs/MINIMAP2_MIGRATION.md`
- [ ] Add troubleshooting section
- [ ] Document any gotchas or caveats
- [ ] Add rollback instructions

---

## Phase 7: Code Cleanup

- [ ] Remove commented-out old code
- [ ] Remove obsolete installation scripts
- [ ] Clean up imports
- [ ] Run linter and fix issues
- [ ] Format code with black/autopep8
- [ ] Remove any backwards compatibility. The code must now be written as if minimap2 never existed.

---

## Phase 9: Validation

### 9.1 Functional Validation
- [ ] Run complete pipeline on test dataset
- [ ] Compare outputs with previous version
- [ ] Verify gene identification is identical
- [ ] Check primer finder results
- [ ] Validate study summary outputs

### 9.2 User Acceptance Testing
- [ ] Have another developer test
- [ ] Test on production-like data
- [ ] Verify performance is acceptable
- [ ] Check error handling

### 9.3 Regression Testing
- [ ] Compare all outputs with baseline
- [ ] Document any intentional differences
- [ ] Ensure no unexpected changes


---

## Notes

### Date Started: [CURRENT DATE]
### Expected Completion: [ESTIMATED DATE]
### Actual Completion: [WHEN FINISHED]

### Issues Encountered:
- 

### Lessons Learned:
- 

### Future Improvements:
- 

---

## Quick Reference

**Key Files**:
- `cfeproviral/utils.py` - Main implementation
- `tests/test_minimap2_migration/` - Test suite
- `docs/MINIMAP2_MIGRATION.md` - Migration plan
- `pyproject.toml` - Dependencies

**Key Functions**:
- `aligner_available()` - Check if aligner is available
- `align()` - Main alignment function
- `_mappy_cigar_to_sam()` - CIGAR conversion
- `_mappy_hit_to_sam()` - Hit to SAM line conversion

**Running Tests**:
```bash
# All migration tests
pytest tests/test_minimap2_migration/ -v

# All tests
pytest tests/ -v

# With coverage
pytest tests/ --cov=cfeproviral --cov-report=html
```

**Building Docker**:
```bash
docker build -t proviral .
docker run --rm proviral --help
```
