import yaml
import os
import csv
from pathlib import Path
import gene_splicer.helpers.yaml_helper
from gene_splicer.logger import logger


# force_all_proviral is for testing, forces samples to be considered proviral
class ProviralHelper:
    def __init__(self, force_all_proviral=False) -> None:
        self.cwd = self.getcwd()
        self.force_all_proviral = force_all_proviral
        self.config = self.load_config()
        self.proviral_samples = self.load_proviral_samples()

    @staticmethod
    def getcwd():
        return Path(os.path.realpath(__file__)).parent

    def load_config(self):
        with open(self.cwd.parent / 'config.yaml') as o:
            return yaml.load(o, Loader=yaml.Loader)

    def load_proviral_samples(self):
        if self.force_all_proviral:
            return
        apath = self.config['RESOURCES']['PROVIRAL_SAMPLES']
        proviral_samples = {}
        with open(apath, newline='') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                if row['is_proviral'] == 'TRUE':
                    proviral_samples[row['sample']] = row
        return proviral_samples

    def is_proviral(self, sample):
        if self.force_all_proviral:
            return True
        elif sample in self.proviral_samples:
            return True
        elif self.config['RESOURCES']['PROVIRAL_PROJECT_CODE'] in sample:
            logger.debug(
                'Sample was deemed proviral by project code but "%s" not in "%s"'
                % (sample, self.config['RESOURCES']['PROVIRAL_SAMPLES']))
            return True

    def pid_from_sample(self, sample):
        try:
            return self.proviral_samples[sample]['pid']
        except KeyError:
            logger.warning('Sample %s not in proviral samples file' % (sample))
