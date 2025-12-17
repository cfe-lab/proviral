import csv
import logging

import pandas as pd
from cfeproviral.version import get_version, get_cfeintact_version
from cfeproviral.primer_finder_errors import PrimerFinderErrors
from cfeproviral.helpers.proviral_helper import ProviralHelper

logger = logging.getLogger(__name__)


class OutcomeSummary:
    def __init__(self, conseqs_df, contigs_df, outpath, force_all_proviral=False):
        self.proviral_helper = ProviralHelper(force_all_proviral)
        self.errors = PrimerFinderErrors()
        self.data = {}
        self.path = outpath / "outcome_summary.csv"
        # To know how many columns to produce
        self.max_failed = 0
        self.create(conseqs_df, contigs_df)

    def add_sample(self, row):
        self.data[row["sample"]] = {
            "sample": row["sample"],
            # Since I add suffix "_conseq", remove it
            "run": row["run_name"].rsplit("_", 1)[0],
            "conseq_passed": False,
            "contig_passed": False,
            "reference": "",
            "seqlen": "",
            "seqtype": "",
            "sequence": "",
            "failed": [],
            "error": "",
        }

    def set_passed(self, row):
        seqtype = row["seqtype"]
        if seqtype.endswith("s"):
            seqtype = seqtype[:-1]
        passed_field = seqtype + "_passed"
        self.data[row["sample"]]["reference"] = row["reference"]
        self.data[row["sample"]]["seqlen"] = row["seqlen"]
        self.data[row["sample"]][passed_field] = True
        self.data[row["sample"]]["sequence"] = row["sequence"]
        self.data[row["sample"]]["seqtype"] = seqtype
        self.data[row["sample"]]["is_rev_comp"] = row["is_rev_comp"]

    def set_failed(self, row, error):
        logger.critical('Sample "%s" already has a passed sequence!' % row["sample"])
        self.data[row["sample"]]["conseq_passed"] = False
        self.data[row["sample"]]["contig_passed"] = False
        self.data[row["sample"]]["sequence"] = ""
        self.data[row["sample"]]["seqtype"] = ""
        self.data[row["sample"]]["error"] = error

    def add_failure(self, row, seqtype):
        nfailed = len(self.data[row["sample"]]["failed"])
        self.data[row["sample"]]["failed"].append(
            {
                f"fail_error_{nfailed}": row["error"],
                f"fail_fwd_err_{nfailed}": row["fwd_error"],
                f"fail_rev_err_{nfailed}": row["rev_error"],
                f"fail_seqtype_{nfailed}": seqtype,
                f"fail_is_rev_comp_{nfailed}": row["is_rev_comp"],
                f"fail_seqlen_{nfailed}": row["seqlen"],
                f"fail_sequence_{nfailed}": row["sequence"],
                f"fail_ref_{nfailed}": row["reference"],
            }
        )

    def handle_edge_cases(self, row):
        # Determine type of error for certain cases
        # Raises a specific exception to tell main loop to continue

        # If conseq is not max, do not record it
        if (
            row["error"] == self.errors.not_max
            or row["error"] == self.errors.non_proviral
        ):
            raise IndexError
        # If the sample had no contigs/conseqs (i.e. remap == 0)
        elif row["error"] == self.errors.no_sequence:
            self.data[row["sample"]]["error"] = row["error"]
            raise IndexError
        # If any reverse or unknown in reference, sample is non-HIV
        elif any([x in row["reference"] for x in ("reverse", "unknown")]):
            row["error"] = self.errors.non_hiv

    def process_conseqs(self, conseqs_df):
        for index, row in conseqs_df.iterrows():
            seqtype = "conseq"
            sample = row["sample"]
            passed = not row["error"] and not row["fwd_error"] and not row["rev_error"]
            # If sample not in data yet
            if row["sample"] not in self.data:
                self.add_sample(row)
                if passed:
                    self.set_passed(row)
                else:
                    try:
                        self.handle_edge_cases(row)
                    except IndexError:
                        continue
                    self.add_failure(row, seqtype=seqtype)
            # Else if sample is already in self.data
            else:
                # If we have already seen a conseq that passed, set them both to
                # fail since this means multiple passing conseqs
                if passed and self.data[sample]["conseq_passed"]:
                    self.set_failed(row, error=self.errors.multiple_passed)
                # If we have not seen a conseq that passed, set to pass
                elif passed:
                    self.set_passed(row)
                # If not passed
                else:
                    try:
                        self.handle_edge_cases(row)
                    except IndexError:
                        continue
                    self.add_failure(row, seqtype=seqtype)

    def process_contigs(self, contigs_df):
        # Go through all the contigs
        for index, row in contigs_df.iterrows():
            sample = row["sample"]
            conseq_data = self.data.get(sample)
            if conseq_data and conseq_data["conseq_passed"]:
                continue

            seqtype = "contig"
            passed = not row["error"] and not row["fwd_error"] and not row["rev_error"]
            # If sample not in data yet, print a warning because we already went
            # through all of the conseqs and so we should have captured every
            # sample
            if row["sample"] not in self.data:
                logger.warning(
                    'Sample "%s" not found in conseqs but was in contigs?!' % sample
                )
                self.add_sample(row)
                if passed:
                    self.set_passed(row)
                else:
                    try:
                        self.handle_edge_cases(row)
                    except IndexError:
                        continue
                    self.add_failure(row, seqtype=seqtype)
            # Else if sample is already in self.data
            else:
                if passed and self.data[sample]["contig_passed"]:
                    self.set_failed(row, error=self.errors.multiple_passed)
                # If the contig passes
                elif passed:
                    self.set_passed(row)
                else:
                    try:
                        self.handle_edge_cases(row)
                    except IndexError:
                        continue
                    self.add_failure(row, seqtype=seqtype)

    def write(self):
        fieldnames = [
            "sample",
            "run",
            "passed",
            "error",
            "reference",
            "seqtype",
            "is_rev_comp",
            "seqlen",
            "sequence",
            "fwd_err",
            "rev_err",
            "cfeproviral_version",
            "cfeintact_version",
        ]

        # Write the rows
        with open(self.path, "w", newline="") as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            for sample in sorted(
                self.data, key=lambda x: x and int(x.rsplit("_")[-1][1:])
            ):
                row_data = {
                    k: v for k, v in self.data[sample].items() if k in fieldnames
                }
                row_data["cfeproviral_version"] = get_version()
                row_data["cfeintact_version"] = get_cfeintact_version()
                # Convert seqlen from float to int if it's a valid number
                if "seqlen" in row_data and row_data["seqlen"] != "":
                    try:
                        # Convert to float first (in case it's a string), then to int
                        seqlen_value = float(row_data["seqlen"])
                        if not pd.isna(seqlen_value):
                            row_data["seqlen"] = int(seqlen_value)
                    except (ValueError, TypeError):
                        # Keep original value if conversion fails
                        pass
                writer.writerow(row_data)

    def reduce(self):
        # Perform some final filtering steps
        for sample in self.data:
            # Set passed if either contig or conseq passes
            sample_passed = self.data[sample]["passed"] = (
                self.data[sample]["conseq_passed"] or self.data[sample]["contig_passed"]
            )

            # Zabrina's request to remove all failed seqs if sample passed
            if sample_passed:
                self.data[sample]["failed"] = []

            # Zabrina's request to simply display a single error if all failures are due to not aligning to HIV
            count_non_hiv = 0
            is_hiv_indicies = []
            for i, fail in enumerate(self.data[sample]["failed"]):
                if fail[f"fail_error_{i}"] == self.errors.non_hiv:
                    count_non_hiv += 1
                else:
                    is_hiv_indicies.append(i)

            # If the number of failures is equal to the count of non-hiv failures then all failures are due to non-hiv
            if (
                count_non_hiv == len(self.data[sample]["failed"])
                and not self.data[sample]["passed"]
                and len(self.data[sample]["failed"]) > 0
            ):
                self.data[sample]["error"] = self.errors.non_hiv
                self.data[sample]["failed"] = []
                continue

            # Otherwise at least one failed sequence was hiv so we should
            # display only the error for that sequence
            elif not sample_passed and self.data[sample]["failed"]:
                new_failed = []
                for i in is_hiv_indicies:
                    new_failed.append(self.data[sample]["failed"][i])
                self.data[sample]["failed"] = new_failed
                for i, fail in enumerate(self.data[sample]["failed"]):
                    newfail = {}
                    for k, v in fail.items():
                        # Reset the number
                        k = k.rsplit("_", 1)[0] + f"_{i}"
                        newfail[k] = v
                        # Why do I need this?
                        # self.data[sample][k] = v
                    self.data[sample]["failed"][i] = newfail

                # Update the max number of failures
                nfailed = len(self.data[sample]["failed"])
                if nfailed > self.max_failed:
                    self.max_failed = nfailed
                self.data[sample]["error"] = self.errors.hiv_but_failed

            # When generating output summary for a sample:
            # Any contigs or conseqs that didn't BLAST to HIV have been removed
            # by earlier code.
            # 1. If conseqs passed, produce output summary from them. (earlier)
            # 2. Else if contigs passed, produce output summary from them. (earlier)
            # 3. Else if there are no conseqs, and contig error is "non-HIV" or
            #    "no sequence", report that.
            # 4. Else if there are no conseqs, generate error "low coverage".
            # 5. Else if sample has only one conseq and it had primer failure,
            #    report primer error
            # 6. A sample that yielded only one conseq, and it had low coverage,
            #    report low coverage
            # 7. A sample that yielded multiple contigs, all of which SOLELY
            #    suffered primer failure -> MULTIPLE CONTIGS
            # 8. A sample that yielded multiple contigs, all of which SOLELY
            #    suffered low coverage ->  MULTIPLE CONTIGS
            # 9. A sample that yielded multiple contigs, all of which failed
            #    either due to primer, or to low coverage ->  MULTIPLE CONTIGS
            conseq_failures = []
            for i, row in enumerate(self.data[sample]["failed"]):
                sequence_type = row[f"fail_seqtype_{i}"]
                if sequence_type == "conseq":
                    j = len(conseq_failures)
                    new_row = {}
                    for old_name, value in row.items():
                        name_parts = old_name.split("_")
                        name_parts[-1] = str(j)
                        new_row["_".join(name_parts)] = value
                    conseq_failures.append(new_row)
            conseq_failure_count = len(conseq_failures)

            # Case 1 and 2
            if sample_passed:
                self.data[sample]["error"] = None
            elif conseq_failure_count == 0:
                contig_errors = self.data[sample]["failed"]
                if not contig_errors:
                    assert self.data[sample]["error"]
                else:
                    first_error = (
                        contig_errors[0]["fail_fwd_err_0"]
                        or contig_errors[0]["fail_rev_err_0"]
                    )
                    if first_error in (self.errors.non_hiv, self.errors.no_sequence):
                        # Case 3
                        self.data[sample]["error"] = first_error
                    else:
                        # Case 4
                        self.data[sample]["error"] = self.errors.low_cov
            elif conseq_failure_count == 1:
                error_row = conseq_failures[0]
                for field_name in (
                    "ref",
                    "seqtype",
                    "is_rev_comp",
                    "seqlen",
                    "sequence",
                    "fwd_err",
                    "rev_err",
                ):
                    target = "reference" if field_name == "ref" else field_name
                    field_name = f"fail_{field_name}_0"
                    self.data[sample][target] = error_row[field_name]
                # Case 5
                primer_errors = (
                    self.errors.no_primer,
                    self.errors.failed_validation,
                    self.errors.low_end_cov,
                )
                if (
                    error_row["fail_fwd_err_0"] in primer_errors
                    or error_row["fail_rev_err_0"] in primer_errors
                ):
                    self.data[sample]["error"] = self.errors.primer_error
                # Case 6
                elif error_row["fail_error_0"] == self.errors.low_internal_cov:
                    self.data[sample]["error"] = self.errors.low_internal_cov
            else:
                # Case 7, 8, 9
                self.data[sample]["error"] = self.errors.multiple_contigs

    def create(self, conseqs_df, contigs_df):
        # Normalize nan to None
        contigs_df = contigs_df.where(pd.notnull(contigs_df), "")
        conseqs_df = conseqs_df.where(pd.notnull(conseqs_df), "")

        self.process_conseqs(conseqs_df)
        self.process_contigs(contigs_df)

        self.reduce()

        self.write()
