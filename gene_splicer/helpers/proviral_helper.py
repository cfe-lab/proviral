import yaml
import os
from pathlib import Path
import helpers.yaml_helper


# force_all_proviral is for testing, forces samples to be considered proviral
class ProviralHelper:
    def __init__(self, force_all_proviral=False) -> None:
        self.force_all_proviral = force_all_proviral
        self.cwd = self.getcwd()
        self.load_config()

    @staticmethod
    def getcwd():
        return Path(os.path.realpath(__file__)).parent

    def load_config(self):
        with open(self.cwd.parent / 'config.yaml') as o:
            self.config = yaml.load(o, Loader=yaml.Loader)

    def get_proviral_samples(self):
        return self.config['RESOURCES']['PROVIRAL_SAMPLES']

    def is_proviral(self, sample):
        return sample in self.config['RESOURCES'][
            'PROVIRAL_SAMPLES'] or self.force_all_proviral

    def get_sample_pid_mapping(self):
        return self.config['RESOURCES']['SAMPLE_PID_MAPPING']
