import shutil
import os
import subprocess
from pathlib import Path

import utils
from logger import logger


class Alignment:
    def __init__(self,
                 target: str,
                 query: str,
                 outpath: Path,
                 aligner_params: dict = {'-a': ''},
                 aligner_path: str = 'minimap2') -> None:
        self.target = target
        self.query = query
        if not outpath:
            self.outpath = Path('/tmp' / utils.Random.upper())
        else:
            self.outpath = outpath
        self.aligner_path = aligner_path
        self.aligner_params = aligner_params
        if self.aligner_available():
            self.align(self.aligner_params)

    def aligner_available(self):
        cmd = [self.aligner_path]
        process = subprocess.run(cmd)
        if process.returncode == 0:
            return True
        # If process is not successful
        logger.error(
            f'Aligner "{self.aligner_path}" not in PATH! Unable to align!')
        return False

    def align(self, aligner_params: dict, clean: bool = True):
        """Align a query sequence to a target sequence

        Args:
            aligner_params (dict): Parameters of the aligner where keys are the flag and values are the value, i.e. `-s 10` becomes: {'-s': 10}. For parameters requiring "=" symbols include them in the key, i.e.: {'--setting=': 10}

            clean (bool): Set this to False to prevent the removal of intermediate files: (query.fasta, target.fasta). Defaults to True

        Returns:
            Path: A path to the alignment file generated
        """
        outpath = self.outpath / 'minimap2_aln'
        if os.path.isdir(outpath):
            shutil.rmtree(outpath)
        os.makedirs(outpath)
        # Write the query fasta
        query_fasta_path = utils.write_fasta({'query': self.query},
                                             outpath / 'query.fasta')
        # Write the target fasta
        target_fasta_path = utils.write_fasta({'target': self.target},
                                              outpath / 'target.fasta')
        aligner_params = [str(k) + str(v) for k, v in aligner_params.items()]
        cmd = [
            self.aligner_path, *aligner_params, target_fasta_path,
            query_fasta_path
        ]
        alignment_path = outpath / 'alignment.sam'
        with open(alignment_path, 'w') as alignment:
            self.process = subprocess.run(cmd, stdout=alignment, errors=True)
        if self.process.returncode != 0:
            logger.error('Alignment failed! Error: %s' % self.process.stderr)
        if clean:

        return alignment_path