class PrimerFinderErrors:
    def __init__(self) -> None:
        self.not_max = 'not MAX'
        self.is_v3 = 'is V3 sequence'
        self.low_internal_cov = 'low internal read coverage'
        self.non_tcga = 'sequence contained non-TCGA/gap'
        self.low_end_cov = 'low end read coverage'
        self.no_primer = 'primer was not found'
        self.failed_validation = 'primer failed validation'
        self.no_sequence = 'no contig/conseq constructed'
        self.non_hiv = 'sequence is non-hiv'
        self.multiple_passed = 'sample has multiple QC-passed sequences'
        self.primer_error = 'primer error'
        self.multiple_contigs = 'multiple contigs'
        self.hiv_but_failed = 'hiv but failed'
        self.non_proviral = 'sequence is non-proviral'