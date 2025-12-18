# Minimap2 to Mappy Migration Checklist

## Phase 1: Preparation ✓

- [x] Analyze current minimap2 usage in codebase
- [x] Identify all files that use minimap2
- [x] Document current behavior and command-line arguments
- [x] Create comprehensive test suite
- [ ] Run baseline tests with current implementation
- [ ] Document baseline test results
- [ ] Save example SAM outputs for comparison

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

## Phase 2: Test Suite Development ✓

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

## Phase 3: Implementation

### 3.1 Add Dependencies
- [ ] Update `pyproject.toml` to include mappy
  ```toml
  dependencies = [
      ...
      "mappy>=2.24",
  ]
  ```
- [ ] Test installation: `pip install mappy`
- [ ] Verify mappy version and compatibility

### 3.2 Implement Helper Functions
- [ ] Implement `_mappy_cigar_to_sam()` in `cfeproviral/utils.py`
- [ ] Implement `_mappy_hit_to_sam()` in `cfeproviral/utils.py`
- [ ] Implement `_generate_sam_header()` in `cfeproviral/utils.py`
- [ ] Test helper functions with `test_sam_conversion.py`

### 3.3 Update Core Functions
- [ ] Update `aligner_available()` function
  - [ ] Add try/except for mappy import
  - [ ] Keep parameter for backward compatibility
  - [ ] Update docstring
- [ ] Update `align()` function
  - [ ] Import mappy
  - [ ] Create aligner with target sequence
  - [ ] Perform alignment
  - [ ] Convert hits to SAM format
  - [ ] Write SAM file with headers
  - [ ] Write log file
  - [ ] Maintain same directory structure
  - [ ] Keep same return value (Path or False)
  - [ ] Keep all parameters for backward compatibility
  - [ ] Update docstring
- [ ] Keep old code commented for reference

### 3.4 Testing During Implementation
- [ ] Run `test_sam_conversion.py` after each helper function
- [ ] Run `test_minimap2_compatibility.py` after updating functions
- [ ] Fix any test failures
- [ ] Add debug logging if needed

---

## Phase 4: Integration Testing

### 4.1 Unit Tests
- [ ] Run all tests in `test_minimap2_migration/`
- [ ] Fix any failures
- [ ] Ensure 100% pass rate

### 4.2 Existing Tests
- [ ] Run `tests/test_utils/test_utils.py`
- [ ] Run `tests/test_gene_splicer/test_gene_splicer.py`
- [ ] Run `tests/test_pipeline/test_pipeline.py`
- [ ] Run full test suite: `pytest tests/`
- [ ] Fix any regressions

### 4.3 Real Data Testing
- [ ] Test with example1 data
- [ ] Test with example2 data
- [ ] Test with example3 data
- [ ] Test with example4 data
- [ ] Compare outputs with previous runs
- [ ] Verify gene splicing results match

### 4.4 Performance Testing
- [ ] Benchmark alignment speed
- [ ] Compare with previous performance
- [ ] Ensure < 20% slowdown
- [ ] Profile memory usage

---

## Phase 5: Update Dependencies

### 5.1 Installation Scripts
- [ ] Remove/comment `scripts/installation/30-minimap2.sh`
- [ ] Update `scripts/install.sh` to not call 30-minimap2.sh
- [ ] Test installation process

### 5.2 CI/CD
- [ ] Update `.github/workflows/python-app.yml`
- [ ] Remove minimap2 download step
- [ ] Verify mappy is installed via `uv sync`
- [ ] Run CI pipeline and verify success

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

---

## Phase 8: Platform Testing

### 8.1 Linux Testing
- [ ] Test on Ubuntu 22.04
- [ ] Test on Ubuntu 20.04
- [ ] Test on Debian if used
- [ ] Verify installation
- [ ] Run full test suite
- [ ] Test with real data

### 8.2 macOS Testing
- [ ] Test on macOS (if available)
- [ ] Verify installation
- [ ] Run test suite
- [ ] Test with real data

### 8.3 Windows Testing
- [ ] Test on Windows 10/11
- [ ] Verify installation (should work now!)
- [ ] Run test suite
- [ ] Test with real data
- [ ] Verify this is a NEW capability

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

## Phase 10: Release

### 10.1 Version Control
- [ ] Create feature branch for migration
- [ ] Commit changes with clear messages
- [ ] Create pull request
- [ ] Code review
- [ ] Merge to main branch

### 10.2 Release Notes
- [ ] Document changes in CHANGELOG
- [ ] Note breaking changes (none expected)
- [ ] List improvements (Windows support!)
- [ ] Update version number

### 10.3 Deployment
- [ ] Tag release
- [ ] Build and publish package
- [ ] Update documentation site
- [ ] Notify users of changes

---

## Phase 11: Post-Release

### 11.1 Monitoring
- [ ] Monitor for issues in first week
- [ ] Check CI/CD pipeline
- [ ] Review user feedback
- [ ] Address any bugs quickly

### 11.2 Documentation
- [ ] Add FAQ entries if needed
- [ ] Update troubleshooting guide
- [ ] Create migration guide for users

### 11.3 Cleanup
- [ ] Remove old minimap2 references
- [ ] Archive old test outputs
- [ ] Clean up migration-specific code

---

## Rollback Plan

If migration fails at any point:

### Immediate Rollback
- [ ] Revert commits
- [ ] Restore minimap2 executable dependency
- [ ] Re-enable installation scripts
- [ ] Rebuild Docker/Singularity containers
- [ ] Notify team

### Partial Rollback (Fallback Mode)
- [ ] Implement fallback to minimap2 executable
- [ ] Add configuration option for tool choice
- [ ] Document when to use each option

---

## Success Metrics

Migration is successful when:

- ✅ All tests pass (existing + new)
- ✅ Performance is acceptable
- ✅ Works on Windows, Linux, macOS
- ✅ No breaking changes for users
- ✅ Documentation is updated
- ✅ CI/CD works without manual minimap2 installation
- ✅ Docker builds successfully
- ✅ Gene splicing results match baseline

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
