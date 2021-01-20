import os
import csv
import pandas as pd
import numpy as np
from pathlib import Path
from logger import logger

import utils


class OutcomeSummary:
    def __init__(self, conseqs_df, contigs_df, outpath):
        self.data = {}
        self.path = outpath / 'outcome_summary.csv'
        # To know how many columns to produce
        self.max_failed = 0
        self.create(conseqs_df, contigs_df)

    def addSample(self, row):
        self.data[row['sample']] = {
            'sample': row['sample'],
            # Since I add suffix "_conseq", remove it
            'run': row['run_name'].rsplit('_', 1)[0],
            'conseq_passed': False,
            'contig_passed': False,
            'reference': None,
            'seqlen': None,
            'seqtype': None,
            'sequence': None,
            'failed': [],
            'error': None
        }

    def setPassed(self, row):
        self.data[row['sample']]['reference'] = row['reference']
        self.data[row['sample']]['seqlen'] = row['seqlen']
        self.data[row['sample']]['conseq_passed'] = True
        self.data[row['sample']]['sequence'] = row['sequence']
        self.data[row['sample']]['seqtype'] = row['seqtype']

    def setFailed(self, row, error):
        logger.critical('Sample "%s" already has a passed sequence!' %
                        row['sample'])
        self.data[row['sample']]['conseq_passed'] = False
        self.data[row['sample']]['contig_passed'] = False
        self.data[row['sample']]['sequence'] = None
        self.data[row['sample']]['seqtype'] = None
        self.data[row['sample']]['error'] = error

    def addFailure(self, row, seqtype):
        nfailed = len(self.data[row['sample']]['failed'])
        self.data[row['sample']]['failed'].append({
            f'fail_error_{nfailed}':
            row['error'],
            f'fail_fwd_err_{nfailed}':
            row['fwd_error'],
            f'fail_rev_err_{nfailed}':
            row['rev_error'],
            f'fail_seqtype_{nfailed}':
            seqtype,
            f'fail_seqlen_{nfailed}':
            row['seqlen'],
            f'fail_sequence_{nfailed}':
            row['sequence'],
            f'fail_ref_{nfailed}':
            row['reference']
        })

    def handleEdgeCases(self, row):
        ## Determine type of error for certain cases
        ## Raises a specific exception to tell main loop to continue
        ## These error strings are set in primer_finder.find_primers() function
        # If conseq is not max, do not record it
        if ((row['error'] == 'contig not MAX')
                or (row['error'] == 'is V3 sequence')):
            raise IndexError
        # If the sample is not proviral (see function in utils for details)
        elif not utils.isProviral(row['sample']):
            row['error'] = 'Sample is non-proviral'
        # If the sample had no contigs/conseqs (i.e. remap == 0)
        elif row['error'] == 'No contig/conseq constructed':
            self.data[row['sample']]['error'] = row['error']
            raise IndexError
        # If any reverse or unknown in reference, sample is non-HIV
        elif any([x in row['reference'] for x in ('reverse', 'unknown')]):
            row['error'] = 'Sample does not align to HIV'

    def processConseqs(self, conseqs_df):
        for index, row in conseqs_df.iterrows():
            seqtype = 'conseq'
            sample = row['sample']
            passed = row['error'] is None and row['fwd_error'] is None and row[
                'rev_error'] is None
            # If sample not in data yet
            if row['sample'] not in self.data:
                self.addSample(row)
                if passed:
                    self.setPassed(row)
            # Else if sample is already in self.data
            else:
                # If we have already seen a conseq that passed, set them both to fail since this means multiple passing conseqs
                if passed and self.data[sample]['conseq_passed']:
                    self.setFailed(
                        row, error='Sample has multiple passed sequences')
                # If we have not seen a conseq that passed, set to pass
                elif passed:
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
            sample = row['sample']
            passed = row['error'] is None and row['fwd_error'] is None and row[
                'rev_error'] is None
            # If sample not in data yet, print a warning because we already went through all of the conseqs and so we should have captured every sample
            if row['sample'] not in self.data:
                logger.warning(
                    'Sample "%s" not found in conseqs but was in contigs?!' %
                    sample)
            # Else if sample is already in self.data
            else:
                if passed and self.data[sample]['contig_passed']:
                    self.setFailed(
                        row, error='Sample has multiple passed sequences')
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
            'sample', 'run', 'passed', 'error', 'reference', 'seqtype',
            'seqlen', 'sequence'
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
            for sample in self.data:
                self.data[sample] = {
                    k: v
                    for k, v in self.data[sample].items() if k in fieldnames
                }
                writer.writerow(self.data[sample])

    def reduce(self):
        # Perform some final filtering steps
        for sample in self.data:
            # Set passed if either contig or conseq passes
            self.data[sample]['passed'] = self.data[sample][
                'conseq_passed'] or self.data[sample]['contig_passed']

            # Zabrina's request to remove all failed seqs if sample passed
            if self.data[sample]['passed']:
                self.data[sample]['failed'] = []

            # Zabrina's request to simply display a single error if all failures are due to not aligning to HIV
            count_non_hiv = 0
            is_hiv_indicies = []
            for i, fail in enumerate(self.data[sample]['failed']):
                if fail[f'fail_error_{i}'] == 'Sample does not align to HIV':
                    count_non_hiv += 1
                else:
                    is_hiv_indicies.append(i)

            # If the number of failures is equal to the count of non-hiv failures then all failures are due to non-hiv
            if count_non_hiv == len(
                    self.data[sample]
                ['failed']) and not self.data[sample]['passed'] and len(
                    self.data[sample]['failed']) > 0:
                self.data[sample]['error'] = 'Sample does not align to HIV'
                self.data[sample]['failed'] = []

            # Otherwise at least one sequence was not non-hiv so we should display only the error for that sequence
            else:
                new_failed = []
                for i in is_hiv_indicies:
                    new_failed.append(self.data[sample]['failed'][i])
                self.data[sample]['failed'] = new_failed
            for fail in self.data[sample]['failed']:
                for k, v in fail.items():
                    self.data[sample][k] = v

            # Update the max number of failures
            nfailed = len(self.data[sample]['failed'])
            if nfailed > self.max_failed:
                self.max_failed = nfailed

    def create(self, conseqs_df, contigs_df):
        # Normalize nan to None
        contigs_df = contigs_df.where(pd.notnull(contigs_df), None)
        conseqs_df = conseqs_df.where(pd.notnull(conseqs_df), None)

        self.processConseqs(conseqs_df)
        self.processContigs(contigs_df)

        self.reduce()

        self.write()