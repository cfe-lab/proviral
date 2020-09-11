import os
import io
import shutil
from contextlib import contextmanager
from pathlib import Path


class File:
    def __init__(self, path: Path = None, mode='r', clean: bool = False):
        try:
            self.path = Path(path)
        except TypeError:
            self.path = None
        self.mode = mode
        self.clean = clean
        self.file = self.get_file()

    def get_file(self):
        file_object = None
        if self.path is None:
            file_object = io.StringIO()
        else:
            try:
                self.path = self.path.resolve()
            except FileNotFoundError:
                # If the parent folder does not exist
                if not os.path.isdir(self.path.parent):
                    self.parent_created = True
                    os.makedirs(self.path.parent)
            file_object = open(self.path, self.mode)
        return file_object

    def cleanup(self):
        if self.path:
            if self.parent_created:
                shutil.rmtree(self.path.parent)
            else:
                os.remove(self.path)

    @contextmanager
    def open(self):
        self.file.seek(0)
        try:
            yield self.file
        finally:
            if self.path is not None:
                self.file.close()
                if self.clean:
                    self.cleanup()
            else:
                self.file.seek(0)

    def write(self, *args, **kwargs):
        return self.file.write(*args, **kwargs)

    def read(self, *args, **kwargs):
        return self.file.read(*args, **kwargs)

    def readlines(self, *args, **kwargs):
        return self.file.readlines(*args, **kwargs)

    def seek(self, *args, **kwargs):
        return self.file.seek(*args, **kwargs)