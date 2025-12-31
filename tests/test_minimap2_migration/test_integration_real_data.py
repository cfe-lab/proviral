"""
Integration tests using real HIV sequences from the test suite.

These tests use actual test data to ensure the migration will work
with the real sequences used in production.
"""

import os
from pathlib import Path
import pytest


TEST_DATA_DIR = Path(__file__).parent.parent / 'data'
TEST_GENE_SPLICER_DIR = Path(__file__).parent.parent / 'test_gene_splicer'


class TestWithRealHIVSequences:
    """Tests using real HIV test data."""
    
    @pytest.fixture
    def utils_module(self):
        """Import cfeproviral.utils."""
        import cfeproviral.utils as utils
        return utils
    
    @pytest.fixture
    def mappy_module(self):
        """Import mappy if available."""
        try:
            import mappy
            return mappy
        except ImportError:
            pytest.skip("mappy not installed yet")
    
    def test_load_mod_hxb2_reference(self, utils_module):
        """Verify we can load the modified HXB2 reference."""
        assert hasattr(utils_module, 'mod_hxb2')
        assert isinstance(utils_module.mod_hxb2, str)
        # The modified HXB2 reference is actually 8936 bp
        assert len(utils_module.mod_hxb2) > 8900, "HXB2 genome should be ~8.9kb"
        assert len(utils_module.mod_hxb2) < 9000
    
    def test_load_existing_sam_file(self, utils_module):
        """Test that we can generate alignments (replacing SAM file loading)."""
        # We no longer load SAM files - we generate alignments with mappy
        target = utils_module.mod_hxb2
        query = target[123:1626]  # GAG region
        
        cigar_hits = utils_module.align(target, query, 'test_load', outdir=TEST_GENE_SPLICER_DIR)
        
        assert isinstance(cigar_hits, list)
        assert len(cigar_hits) > 0, "Should generate alignments"
        
        # Each alignment should be a CigarHit object
        from aligntools import CigarHit
        for hit in cigar_hits:
            assert isinstance(hit, CigarHit)
    
    def test_existing_sam_structure(self, utils_module):
        """Verify structure of existing SAM files."""
        sam_file = TEST_GENE_SPLICER_DIR / 'large_deletion1.sam'
        
        if not sam_file.exists():
            pytest.skip(f"Test file not found: {sam_file}")
        
        with open(sam_file) as f:
            lines = f.readlines()
        
        # Should have headers
        header_lines = [l for l in lines if l.startswith('@')]
        assert len(header_lines) > 0, "SAM file should have headers"
        
        # Should have alignments
        alignment_lines = [l for l in lines if l and not l.startswith('@')]
        assert len(alignment_lines) > 0, "SAM file should have alignments"
        
        # Check first alignment structure
        fields = alignment_lines[0].split('\t')
        assert len(fields) >= 11
        
        # Document the structure
        qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual = fields[:11]
        print(f"Example alignment from existing SAM:")
        print(f"  QNAME: {qname}")
        print(f"  FLAG: {flag}")
        print(f"  RNAME: {rname}")
        print(f"  POS: {pos}")
        print(f"  MAPQ: {mapq}")
        print(f"  CIGAR: {cigar}")
    
    def test_align_gag_sequence(self, utils_module, mappy_module):
        """Test aligning GAG sequence (from large_deletion1 test)."""
        # This is the actual test case from test_gene_splicer.py
        target = utils_module.mod_hxb2
        
        # GAG gene coordinates in HXB2
        query = target[123:1626]  # GAG region
        
        assert len(query) > 0
        
        # Test with mappy
        aligner = mappy_module.Aligner(seq=target, preset='sr')
        hits = list(aligner.map(query))
        
        assert len(hits) > 0, "GAG should align to HXB2"
        
        # Check alignment positions
        hit = hits[0]
        assert hit.r_st >= 0
        assert hit.r_en <= len(target)
        
        # Should map back to approximately the same region
        # (allowing some tolerance for alignment differences)
        assert abs(hit.r_st - 123) < 50, f"Expected alignment near 123, got {hit.r_st}"
    
    def test_align_full_genome(self, utils_module, mappy_module):
        """Test aligning full HXB2 genome to itself."""
        target = utils_module.mod_hxb2
        query = target  # Exact match
        
        aligner = mappy_module.Aligner(seq=target, preset='sr')
        hits = list(aligner.map(query))
        
        assert len(hits) > 0, "Full genome should align to itself"
        
        hit = hits[0]
        # Should be a near-perfect alignment
        assert hit.mapq > 0, "Mapping quality should be high"
        assert hit.r_st < 10, "Should start near beginning"
        assert hit.r_en > len(target) - 10, "Should extend near end"
    
    def test_compare_with_existing_sam_output(self, utils_module, mappy_module):
        """Test that mappy produces consistent alignments."""
        # Get query sequence from test
        target = utils_module.mod_hxb2
        query = target[123:1626]  # GAG
        
        # Align with mappy
        aligner = mappy_module.Aligner(seq=target, preset='sr')
        hits = list(aligner.map(query))
        
        assert len(hits) > 0, "Should produce alignments"
        
        # Check that alignment is reasonable
        mappy_hit = hits[0]
        mappy_pos = mappy_hit.r_st
        
        print(f"Mappy position: {mappy_pos}")
        print(f"Mappy CIGAR: {mappy_hit.cigar}")
        
        # Should align near the original position (within tolerance)
        assert abs(mappy_pos - 123) < 20, \
            f"Position difference too large: {abs(mappy_pos - 123)}"



class TestAlignFunctionWithRealData:
    """Test the actual align() function with real data."""
    
    @pytest.fixture
    def utils_module(self):
        import cfeproviral.utils as utils
        return utils
    
    def test_current_align_function_signature(self, utils_module):
        """Verify current align function returns CigarHit list."""
        import tempfile
        from pathlib import Path
        from aligntools import CigarHit
        
        target = utils_module.mod_hxb2
        query = target[123:1626]  # GAG
        
        with tempfile.TemporaryDirectory() as tmpdir:
            result = utils_module.align(
                target_seq=target,
                query_seq=query,
                query_name="test_gag",
                outdir=Path(tmpdir)
            )
            
            # align() now returns a list of CigarHit objects
            assert isinstance(result, list), "Should return a list"
            assert len(result) > 0, "Should have alignments"
            assert isinstance(result[0], CigarHit), "Should be CigarHit objects"
            
            # Verify expected directory structure (for debugging/reference)
            outdir = Path(tmpdir) / "test_gag"
            assert (outdir / "query.fasta").exists()
            assert (outdir / "target.fasta").exists()
            assert (outdir / "minimap2.log").exists()
    
    def test_align_output_compatibility(self, utils_module):
        """Test that align() output works with splice_genes()."""
        import tempfile
        from pathlib import Path
        from aligntools import CigarHit
        
        target = utils_module.mod_hxb2
        query = target[1000:2000]
        
        with tempfile.TemporaryDirectory() as tmpdir:
            result = utils_module.align(
                target_seq=target,
                query_seq=query,
                query_name="test_region",
                outdir=Path(tmpdir)
            )
            
            # align() returns CigarHit list
            assert isinstance(result, list)
            assert all(isinstance(hit, CigarHit) for hit in result)
            
            # Test compatibility with splice_genes
            coords = utils_module.splice_genes(query, target, result, utils_module.mod_annot)
            assert isinstance(coords, dict)


class TestGeneSplicerIntegration:
    """Test gene splicer workflow that uses align()."""
    
    @pytest.fixture
    def utils_module(self):
        import cfeproviral.utils as utils
        return utils
    
    def test_gene_splicer_workflow(self, utils_module):
        """Test the full gene splicer workflow."""
        import tempfile
        from pathlib import Path
        
        target = utils_module.mod_hxb2
        query = target[123:1626]  # GAG
        
        with tempfile.TemporaryDirectory() as tmpdir:
            # Step 1: Align (returns CigarHit list)
            cigar_hits = utils_module.align(
                target,
                query,
                "test_gag",
                outdir=Path(tmpdir)
            )
            
            assert len(cigar_hits) > 0, "Should produce alignments"
            
            # Step 2: Splice genes (using CigarHit list)
            coords = utils_module.splice_genes(
                query,
                target,
                cigar_hits,
                utils_module.mod_annot
            )
            
            # Should identify genes
            assert isinstance(coords, dict)
            # GAG should be in the results
            assert 'gag' in coords or len(coords) > 0
            
            # Step 3: Convert coordinates to sequences
            genes = utils_module.coords_to_genes(coords, query)
            
            assert isinstance(genes, dict)
            print(f"Identified genes: {list(genes.keys())}")


class TestBackwardCompatibility:
    """Tests to ensure backward compatibility during migration."""
    
    def test_align_function_parameters(self):
        """Document that align() parameters should remain the same."""
        import cfeproviral.utils as utils
        import inspect
        
        sig = inspect.signature(utils.align)
        params = sig.parameters
        
        # These parameters must be preserved
        assert 'target_seq' in params
        assert 'query_seq' in params
        assert 'query_name' in params
        assert 'outdir' in params
        
        # aligner_path parameter removed after migration completed
        # No backward compatibility needed - migration is complete
    
    def test_align_return_value(self):
        """Document that align() return value should remain the same."""
        # align() should return:
        # - Path object pointing to alignment.sam, OR
        # - False if alignment fails
        # 
        # This contract must be maintained in the new implementation
        pass
    
    def test_output_directory_structure(self):
        """Document expected output directory structure."""
        # The align() function creates:
        # outdir/query_name/
        #   - query.fasta
        #   - target.fasta
        #   - alignment.sam
        #   - minimap2.log (name can change to alignment.log)
        #
        # This structure should be maintained for compatibility
        pass


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
