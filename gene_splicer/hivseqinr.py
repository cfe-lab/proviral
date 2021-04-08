import shutil
import os
import subprocess
import requests
import zipfile
import sys
import io
from gene_splicer.logger import logger


class Hivseqinr:
    def __init__(self, outpath, fasta):
        self.url = 'https://github.com/guineverelee/HIVSeqinR/raw/master/HIVSeqinR_ver2.7.1.zip'
        self.outpath = outpath
        self.fasta = fasta
        self.download()
        self.make_blast_dir()
        self.fix_wds()
        self.copy_fasta()
        self.job = self.run()
        self.finalize()

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
        dbdir = self.outpath / 'hxb2_blast_db'
        try:
            os.mkdir(dbdir)
        except FileExistsError:
            pass
        hxb2_path = dbdir / 'HXB2.fasta'
        shutil.copyfile(self.outpath / 'R_HXB2.fasta', hxb2_path)
        cmd = [
            'makeblastdb', '-in', hxb2_path, '-parse_seqids', '-dbtype', 'nucl'
        ]
        self.dbdir = dbdir
        job = subprocess.run(cmd)
        return job

    def download(self):
        response = requests.get(self.url)
        file_like_object = io.BytesIO(response.content)
        zipfile_obj = zipfile.ZipFile(file_like_object)
        if not os.path.isdir(self.outpath):
            os.makedirs(self.outpath)
        zipfile_obj.extractall(self.outpath)
        return True

    def fix_wds(self):
        rscript_path = os.path.join(
            self.outpath, 'R_HIVSeqinR_Combined_ver09_ScrambleFix.R')
        new_rscript_path = os.path.join(self.outpath, 'modified.R')
        self.new_rscript_path = new_rscript_path
        with open(self.new_rscript_path, 'w') as outfile:
            with open(rscript_path, 'r') as infile:
                for line in infile:
                    if 'MyWD <- getwd()' in line:
                        line = line.replace('MyWD <- getwd()',
                                            f'MyWD = "{self.outpath}"\n')
                    elif line.startswith('MyBlastnDir <-'):
                        line = f'MyBlastnDir = "{str(self.dbdir) + os.path.sep*2}"\n'
                    outfile.write(line)

    def run(self):
        cwd = os.getcwd()
        os.chdir(self.outpath)
        cmd = ['Rscript', './modified.R']
        job = subprocess.run(cmd,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        print(self.clean_output(job.stdout))
        print(self.clean_output(job.stderr), file=sys.stderr)
        os.chdir(cwd)
        return job

    @staticmethod
    def clean_output(output):
        return ''.join(
            [i if ord(i) < 128 else '\n' for i in output.decode('utf')])

    def finalize(self):
        path = self.outpath / 'HIVSEQINR_COMPLETE'
        path.touch()
