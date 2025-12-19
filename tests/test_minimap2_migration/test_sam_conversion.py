"""
Unit tests for SAM format conversion utilities.

These tests focus on the conversion between mappy's internal format
and SAM format strings that the rest of the codebase expects.
"""

import pytest


class TestCigarConversion:
    """Tests for CIGAR string conversion from mappy format to SAM."""
    
    def test_cigar_match_only(self):
        """Test simple match CIGAR."""
        # mappy format: [(100, 0)] means 100 matches (operation 0 = M)
        # Expected SAM: "100M"
        cigar_tuples = [(100, 0)]
        result = _convert_cigar_to_sam(cigar_tuples, strand=1, query_len=100, q_st=0, q_en=100)
        assert result == "100M"
    
    def test_cigar_with_insertion(self):
        """Test CIGAR with insertion."""
        # 50M, 5I, 45M
        cigar_tuples = [(50, 0), (5, 1), (45, 0)]
        result = _convert_cigar_to_sam(cigar_tuples, strand=1, query_len=100, q_st=0, q_en=100)
        assert result == "50M5I45M"
    
    def test_cigar_with_deletion(self):
        """Test CIGAR with deletion."""
        # 50M, 5D, 50M
        cigar_tuples = [(50, 0), (5, 2), (50, 0)]
        result = _convert_cigar_to_sam(cigar_tuples, strand=1, query_len=100, q_st=0, q_en=100)
        assert result == "50M5D50M"
    
    def test_cigar_with_soft_clipping_start(self):
        """Test CIGAR with soft clipping at start."""
        # Query starts at position 10, so 10S at start
        cigar_tuples = [(90, 0)]
        result = _convert_cigar_to_sam(cigar_tuples, strand=1, query_len=100, q_st=10, q_en=100)
        assert result == "10S90M"
    
    def test_cigar_with_soft_clipping_end(self):
        """Test CIGAR with soft clipping at end."""
        # Query ends at position 90, so 10S at end
        cigar_tuples = [(90, 0)]
        result = _convert_cigar_to_sam(cigar_tuples, strand=1, query_len=100, q_st=0, q_en=90)
        assert result == "90M10S"
    
    def test_cigar_with_soft_clipping_both_ends(self):
        """Test CIGAR with soft clipping at both ends."""
        cigar_tuples = [(80, 0)]
        result = _convert_cigar_to_sam(cigar_tuples, strand=1, query_len=100, q_st=10, q_en=90)
        assert result == "10S80M10S"
    
    def test_cigar_complex(self):
        """Test complex CIGAR string."""
        # 5S 40M 3I 35M 2D 15M 5S
        cigar_tuples = [(40, 0), (3, 1), (35, 0), (2, 2), (15, 0)]
        result = _convert_cigar_to_sam(cigar_tuples, strand=1, query_len=100, q_st=5, q_en=95)
        assert result == "5S40M3I35M2D15M5S"
    
    def test_cigar_empty(self):
        """Test empty CIGAR (unmapped)."""
        cigar_tuples = []
        result = _convert_cigar_to_sam(cigar_tuples, strand=1, query_len=100, q_st=0, q_en=100)
        assert result == "*"
    
    def test_cigar_with_skipped_region(self):
        """Test CIGAR with skipped region (N operation)."""
        # For spliced alignments: 50M 1000N 50M
        cigar_tuples = [(50, 0), (1000, 3), (50, 0)]
        result = _convert_cigar_to_sam(cigar_tuples, strand=1, query_len=100, q_st=0, q_en=100)
        assert result == "50M1000N50M"


class TestSamFlagCalculation:
    """Tests for SAM FLAG field calculation."""
    
    def test_flag_forward_primary(self):
        """Test flag for forward strand primary alignment."""
        # Primary alignment on forward strand
        flag = _calculate_sam_flag(strand=1, is_primary=True, is_mapped=True)
        assert flag == 0
    
    def test_flag_reverse_primary(self):
        """Test flag for reverse strand primary alignment."""
        # Primary alignment on reverse strand
        flag = _calculate_sam_flag(strand=-1, is_primary=True, is_mapped=True)
        assert flag == 16  # Reverse strand flag
    
    def test_flag_forward_secondary(self):
        """Test flag for forward strand secondary alignment."""
        flag = _calculate_sam_flag(strand=1, is_primary=False, is_mapped=True)
        assert flag == 256  # Secondary alignment flag
    
    def test_flag_reverse_secondary(self):
        """Test flag for reverse strand secondary alignment."""
        flag = _calculate_sam_flag(strand=-1, is_primary=False, is_mapped=True)
        assert flag == 16 | 256  # Reverse + Secondary
    
    def test_flag_unmapped(self):
        """Test flag for unmapped read."""
        flag = _calculate_sam_flag(strand=1, is_primary=True, is_mapped=False)
        assert flag == 4  # Unmapped flag


class TestSamLineConstruction:
    """Tests for complete SAM line construction."""
    
    def test_sam_line_basic(self):
        """Test basic SAM line construction."""
        hit_data = {
            'q_st': 0,
            'q_en': 100,
            'r_st': 500,
            'r_en': 600,
            'mapq': 60,
            'cigar': [(100, 0)],
            'strand': 1,
            'is_primary': True,
        }
        
        sam_line = _construct_sam_line(
            hit_data,
            query_name="test_query",
            query_seq="A" * 100,
            reference_name="ref"
        )
        
        fields = sam_line.split('\t')
        assert len(fields) >= 11, "SAM line should have at least 11 fields"
        assert fields[0] == "test_query"  # QNAME
        assert fields[1] == "0"  # FLAG (forward primary)
        assert fields[2] == "ref"  # RNAME
        assert fields[3] == "501"  # POS (1-based)
        assert fields[4] == "60"  # MAPQ
        assert fields[5] == "100M"  # CIGAR
    
    def test_sam_line_reverse_strand(self):
        """Test SAM line for reverse strand alignment."""
        hit_data = {
            'q_st': 0,
            'q_en': 100,
            'r_st': 500,
            'r_en': 600,
            'mapq': 60,
            'cigar': [(100, 0)],
            'strand': -1,
            'is_primary': True,
        }
        
        sam_line = _construct_sam_line(
            hit_data,
            query_name="test_query",
            query_seq="A" * 100,
            reference_name="ref"
        )
        
        fields = sam_line.split('\t')
        assert fields[1] == "16", "FLAG should indicate reverse strand"
    
    def test_sam_line_with_soft_clipping(self):
        """Test SAM line with soft clipped regions."""
        hit_data = {
            'q_st': 10,
            'q_en': 90,
            'r_st': 500,
            'r_en': 580,
            'mapq': 60,
            'cigar': [(80, 0)],
            'strand': 1,
            'is_primary': True,
        }
        
        sam_line = _construct_sam_line(
            hit_data,
            query_name="test_query",
            query_seq="A" * 100,
            reference_name="ref"
        )
        
        fields = sam_line.split('\t')
        assert fields[5] == "10S80M10S", "CIGAR should include soft clipping"
    
    def test_sam_line_unmapped(self):
        """Test SAM line for unmapped read."""
        sam_line = _construct_unmapped_sam_line(
            query_name="test_query",
            query_seq="A" * 100
        )
        
        fields = sam_line.split('\t')
        assert len(fields) >= 11
        assert fields[1] == "4", "FLAG should indicate unmapped"
        assert fields[2] == "*", "RNAME should be *"
        assert fields[3] == "0", "POS should be 0"
        assert fields[4] == "0", "MAPQ should be 0"
        assert fields[5] == "*", "CIGAR should be *"


class TestSamHeaderGeneration:
    """Tests for SAM header generation."""
    
    def test_header_minimal(self):
        """Test minimal SAM header."""
        header = _generate_sam_header(
            reference_name="ref",
            reference_length=10000
        )
        
        lines = header.strip().split('\n')
        assert len(lines) >= 2
        assert lines[0].startswith("@HD")
        assert "@SQ" in header
        assert "ref" in header
        assert "10000" in header
    
    def test_header_with_program_info(self):
        """Test SAM header with program information."""
        header = _generate_sam_header(
            reference_name="ref",
            reference_length=10000,
            program_id="mappy",
            program_version="2.24"
        )
        
        assert "@PG" in header
        assert "mappy" in header
        assert "2.24" in header
    
    def test_header_format(self):
        """Test SAM header format compliance."""
        header = _generate_sam_header(
            reference_name="MOD_HXB2",
            reference_length=9719
        )
        
        lines = header.strip().split('\n')
        for line in lines:
            assert line.startswith('@'), "Header lines should start with @"
            assert '\t' in line, "Header fields should be tab-separated"


class TestReverseComplement:
    """Tests for reverse complement function (needed for reverse strand)."""
    
    def test_reverse_complement_simple(self):
        """Test reverse complement of simple sequence."""
        seq = "ATCG"
        result = _reverse_complement(seq)
        assert result == "CGAT"
    
    def test_reverse_complement_all_bases(self):
        """Test all DNA bases."""
        seq = "ATCGATCG"
        result = _reverse_complement(seq)
        assert result == "CGATCGAT"
    
    def test_reverse_complement_with_n(self):
        """Test reverse complement with N."""
        seq = "ATCGN"
        result = _reverse_complement(seq)
        assert result == "NCGAT"
    
    def test_reverse_complement_palindrome(self):
        """Test palindromic sequence."""
        seq = "GAATTC"  # EcoRI site
        result = _reverse_complement(seq)
        assert result == "GAATTC"


# Helper function stubs (these would be implemented in utils.py during migration)

def _convert_cigar_to_sam(cigar_tuples, strand, query_len, q_st, q_en):
    """
    Convert mappy CIGAR format to SAM CIGAR string.
    This is a stub for testing - actual implementation goes in utils.py
    """
    if not cigar_tuples:
        return '*'
    
    op_map = {0: 'M', 1: 'I', 2: 'D', 3: 'N', 4: 'S', 5: 'H', 6: 'P', 7: '=', 8: 'X'}
    
    cigar_str = []
    
    if q_st > 0:
        cigar_str.append(f"{q_st}S")
    
    for length, op in cigar_tuples:
        cigar_str.append(f"{length}{op_map.get(op, 'M')}")
    
    if q_en < query_len:
        cigar_str.append(f"{query_len - q_en}S")
    
    return ''.join(cigar_str)


def _calculate_sam_flag(strand, is_primary, is_mapped):
    """Calculate SAM FLAG field."""
    flag = 0
    if not is_mapped:
        flag |= 4
    if strand == -1:
        flag |= 16
    if not is_primary:
        flag |= 256
    return flag


def _construct_sam_line(hit_data, query_name, query_seq, reference_name):
    """Construct complete SAM line from hit data."""
    flag = _calculate_sam_flag(
        hit_data['strand'],
        hit_data['is_primary'],
        is_mapped=True
    )
    
    cigar = _convert_cigar_to_sam(
        hit_data['cigar'],
        hit_data['strand'],
        len(query_seq),
        hit_data['q_st'],
        hit_data['q_en']
    )
    
    pos = hit_data['r_st'] + 1  # 1-based
    
    seq = query_seq
    if hit_data['strand'] == -1:
        seq = _reverse_complement(seq)
    
    return f"{query_name}\t{flag}\t{reference_name}\t{pos}\t{hit_data['mapq']}\t{cigar}\t*\t0\t0\t{seq}\t*"


def _construct_unmapped_sam_line(query_name, query_seq):
    """Construct SAM line for unmapped read."""
    return f"{query_name}\t4\t*\t0\t0\t*\t*\t0\t0\t{query_seq}\t*"


def _generate_sam_header(reference_name, reference_length, program_id=None, program_version=None):
    """Generate SAM header."""
    header = "@HD\tVN:1.0\tSO:unsorted\n"
    header += f"@SQ\tSN:{reference_name}\tLN:{reference_length}\n"
    
    if program_id and program_version:
        header += f"@PG\tID:{program_id}\tPN:{program_id}\tVN:{program_version}\n"
    
    return header


def _reverse_complement(seq):
    """Get reverse complement of DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return ''.join(complement.get(base, base) for base in reversed(seq))


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
