import logging
import shutil
import os
import subprocess
from pathlib import Path

import requests
import zipfile
import io

logger = logging.getLogger(__name__)
SOURCE_URL = 'https://github.com/guineverelee/HIVSeqinR/raw/master/HIVSeqinR_ver2.7.1.zip'


class Hivseqinr:
    def __init__(self, source_path: Path, outpath: Path, fasta: Path):
        self.source_path = source_path.resolve()
        self.blast_db = self.source_path / 'hxb2_blast_db'
        self.outpath = outpath
        self.fasta = fasta

    def copy_fasta(self):
        raw_fastas_path = self.outpath / 'RAW_FASTA'
        try:
            os.makedirs(raw_fastas_path)
        except FileExistsError:
            pass
        logger.debug('Attempting to copy "%s" to "%s"' %
                     (self.fasta, raw_fastas_path))
        shutil.copy(self.fasta, raw_fastas_path)
        return True

    def make_blast_dir(self):
        self.blast_db.mkdir(exist_ok=True)
        hxb2_path = self.blast_db / 'HXB2.fasta'
        shutil.copyfile(self.source_path / 'R_HXB2.fasta', hxb2_path)
        cmd = [
            'makeblastdb', '-in', hxb2_path, '-parse_seqids', '-dbtype', 'nucl'
        ]
        log_path = self.blast_db / 'blast.log'
        with log_path.open('w') as log_file:
            subprocess.run(cmd, stdout=log_file, check=True)
        log_path.unlink()

    def download(self):
        license_path = self.source_path / 'LICENSE.txt'
        if license_path.exists():
            return
        response = requests.get(SOURCE_URL)
        file_like_object = io.BytesIO(response.content)
        zipfile_obj = zipfile.ZipFile(file_like_object)
        if not os.path.isdir(self.source_path):
            os.makedirs(self.source_path)
        zipfile_obj.extractall(self.source_path)
        self.make_blast_dir()
        self.fix_wds()

    def fix_wds(self):
        rscript_path = self.source_path / 'R_HIVSeqinR_Combined_ver09_ScrambleFix.R'
        new_rscript_path = self.source_path / 'modified.R'
        with open(new_rscript_path, 'w') as outfile:
            with open(rscript_path, 'r') as infile:
                for line in infile:
                    if False and 'MyWD <- getwd()' in line:
                        line = line.replace('MyWD <- getwd()',
                                            f'MyWD = "{self.outpath}"\n')
                    elif line.startswith('MyBlastnDir <-'):
                        line = f'MyBlastnDir = "{str(self.blast_db) + os.path.sep*2}"\n'
                    outfile.write(line)

    def run(self):
        self.download()
        self.copy_fasta()
        outpath = Path(self.outpath)
        cmd = ['Rscript', self.source_path / 'modified.R']
        log_path = outpath / 'hivseqinr.log'
        error_path = outpath / 'hivseqinr_error.log'
        with log_path.open('w') as log_file, error_path.open('w') as error_file:
            subprocess.run(cmd,
                           stdout=log_file,
                           stderr=error_file,
                           cwd=self.outpath,
                           check=True)
        self.finalize()

    @staticmethod
    def clean_output(output):
        return ''.join(
            [i if ord(i) < 128 else '\n' for i in output.decode('utf')])

    def finalize(self):
        path = self.outpath / 'HIVSEQINR_COMPLETE'
        path.touch()
