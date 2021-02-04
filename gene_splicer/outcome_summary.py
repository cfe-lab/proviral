import os
import csv
import pandas as pd
import numpy as np
from pathlib import Path
from gene_splicer.logger import logger
import gene_splicer.primer_finder as primer_finder

import gene_splicer.utils as utils


class OutcomeSummary:
    def __init__(self, conseqs_df, contigs_df, outpath):
        self.errors = primer_finder.PFErrors()
        self.data = {}
        self.path = outpath
        # To know how many columns to produce
        self.max_failed = 0
        self.create(conseqs_df, contigs_df)

    def addSample(self, row):
        self.data = {
            'conseq_passed': False,
            'contig_passed': False,
            'reference': '',
            'seqlen': '',
            'seqtype': '',
            'sequence': '',
            'failed': [],
            'error': ''
        }

    def setPassed(self, row):
        self.data['reference'] = row['reference']
        self.data['seqlen'] = row['seqlen']
        self.data['conseq_passed'] = True
        self.data['sequence'] = row['sequence']
        self.data['seqtype'] = row['seqtype']

    def setFailed(self, row, error):
        logger.critical('Sample "%s" already has a passed sequence!' %
                        row['sample'])
        self.data['conseq_passed'] = False
        self.data['contig_passed'] = False
        self.data['sequence'] = ''
        self.data['seqtype'] = ''
        self.data['error'] = error

    def addFailure(self, row, seqtype):
        nfailed = len(self.data['failed'])
        self.data['failed'].append({
            f'fail_error_{nfailed}': row['error'],
            f'fail_fwd_err_{nfailed}': row['fwd_error'],
            f'fail_rev_err_{nfailed}': row['rev_error'],
            f'fail_seqtype_{nfailed}': seqtype,
            f'fail_seqlen_{nfailed}': row['seqlen'],
            f'fail_sequence_{nfailed}': row['sequence'],
            f'fail_ref_{nfailed}': row['reference']
        })

    def handleEdgeCases(self, row):
        ## Determine type of error for certain cases
        ## Raises a specific exception to tell main loop to continue
        ## These error strings are set in primer_finder.find_primers() function
        # If conseq is not max, do not record it
        if ((row['error'] == self.errors.not_max)
                or (row['error'] == self.errors.is_v3)):
            raise IndexError
        # If the sample is not proviral (see function in utils for details)
        # Remove this since the dev version should already only process NFLHIVDNA
        # elif not utils.isProviral(row['sample']):
        #     row['error'] = 'Sample is non-proviral'
        # If the sample had no contigs/conseqs (i.e. remap == 0)
        elif row['error'] == self.errors.no_sequence:
            self.data['error'] = row['error']
            raise IndexError
        # If any reverse or unknown in reference, sample is non-HIV
        elif any([x in row['reference'] for x in ('reverse', 'unknown')]):
            row['error'] = self.errors.non_hiv

    def processConseqs(self, conseqs_df):
        for index, row in conseqs_df.iterrows():
            seqtype = 'conseq'
            passed = not row['error'] and not row['fwd_error'] and not row[
                'rev_error']
            self.addSample(row)
            if passed:
                self.setPassed(row)
            # If not passed
            else:
                try:
                    self.handleEdgeCases(row)
                except IndexError:
                    continue
                self.addFailure(row, seqtype=seqtype)

    def processContigs(self, contigs_df):
        # Go through all the contigs
        for index, row in contigs_df.iterrows():
            seqtype = 'contig'
            passed = not row['error'] and not row['fwd_error'] and not row[
                'rev_error']
            # If sample not in data yet, print a warning because we already went through all of the conseqs and so we should have captured every sample
            if passed and self.data['contig_passed']:
                self.setFailed(row, error=self.errors.multiple_passed)
            # If the contig passes
            elif passed:
                self.setPassed(row)
            else:
                try:
                    self.handleEdgeCases(row)
                except IndexError:
                    continue
                self.addFailure(row, seqtype=seqtype)

    def write(self):
        fieldnames = [
            'passed', 'error', 'reference', 'seqtype', 'seqlen', 'sequence'
        ]
        for i in range(self.max_failed):
            fieldnames += [
                f'fail_error_{i}', f'fail_fwd_err_{i}', f'fail_rev_err_{i}',
                f'fail_seqtype_{i}', f'fail_seqlen_{i}', f'fail_sequence_{i}',
                f'fail_ref_{i}'
            ]

        # Write the rows
        with open(self.path, 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            self.data = {k: v for k, v in self.data.items() if k in fieldnames}
            writer.writerow(self.data)

    # Perform some final filtering steps
    def reduce(self):
        # Set passed if either contig or conseq passes
        self.data['passed'] = self.data['conseq_passed'] or self.data[
            'contig_passed']

        # Zabrina's request to remove all failed seqs if sample passed
        if self.data['passed']:
            self.data['failed'] = []

        # Zabrina's request to simply display a single error if all failures are due to not aligning to HIV
        count_non_hiv = 0
        is_hiv_indicies = []
        for i, fail in enumerate(self.data['failed']):
            if fail[f'fail_error_{i}'] == self.errors.non_hiv:
                count_non_hiv += 1
            else:
                is_hiv_indicies.append(i)

        # If the number of failures is equal to the count of non-hiv failures then all failures are due to non-hiv
        if count_non_hiv == len(
                self.data['failed']) and not self.data['passed'] and len(
                    self.data['failed']) > 0:
            self.data['error'] = self.errors.non_hiv
            self.data['failed'] = []

        # Otherwise at least one sequence was not non-hiv so we should display only the error for that sequence
        else:
            new_failed = []
            for i in is_hiv_indicies:
                new_failed.append(self.data['failed'][i])
            self.data['failed'] = new_failed
        for i, fail in enumerate(self.data['failed']):
            newfail = {}
            for k, v in fail.items():
                # Reset the number
                k = k.rsplit('_', 1)[0] + f'_{i}'
                newfail[k] = v
                self.data[k] = v
            self.data['failed'][i] = newfail

        # Update the max number of failures
        nfailed = len(self.data['failed'])
        if nfailed > self.max_failed:
            self.max_failed = nfailed

        # Natalie's requests
        # 1. If sample has only one contig and it had primer failure -> primer error
        # 2. A sample that yielded only one contig, and it had low coverage -> LOW COVERAGE
        # 3. A sample that yielded multiple contigs, all of which SOLELY suffered primer failure -> MULTIPLE CONTIGS
        # 4. A sample that yielded multiple contigs, all of which SOLELY suffered low coverage ->  MULTIPLE CONTIGS
        # 5. A sample that yielded multiple contigs, some of which failed due to primer, some due to low coverage ->  MULTIPLE CONTIGS

        # Case 1 and 2
        if len(self.data['failed']) == 1 and not self.data['passed']:
            # Case 1
            if 'primer' in self.data['failed'][0][
                    'fail_fwd_err_0'] or 'primer' in self.data['failed'][0][
                        'fail_rev_err_0']:
                self.data['error'] = self.errors.primer_error
            # Case 2
            elif 'coverage' in self.data['failed'][0][
                    'fail_fwd_err_0'] or 'coverage' in self.data['failed'][0][
                        'fail_rev_error_0']:
                self.data['error'] = self.errors.low_end_cov
        # Case 3, 4, and 5
        elif not self.data['passed'] and len(self.data['failed']) > 1:
            # I can just set the error to multiple contigs?
            self.data['error'] = self.errors.multiple_contigs

    def create(self, conseqs_df, contigs_df):
        # Normalize nan to None
        contigs_df = contigs_df.where(pd.notnull(contigs_df), '')
        conseqs_df = conseqs_df.where(pd.notnull(conseqs_df), '')

        self.processConseqs(conseqs_df)
        self.processContigs(contigs_df)

        self.reduce()

        self.write()