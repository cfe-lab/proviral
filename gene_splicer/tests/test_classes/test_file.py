import os
import sys
from pathlib import Path

cwd = Path(os.path.realpath(__file__)).parent
sys.path.append(str(cwd.parent.parent))

from file import File


def test_file_inmem():
    # Should create an in-memory file
    afile = File()
    # Test write
    afile.write('wow\n')
    afile.seek(0)
    # Test readlines
    lines = afile.readlines()
    assert lines == ['wow\n']
    # Reset cursor to top of file
    # afile.seek(0)
    # Test context manager
    with afile.open() as f:
        for line in f:
            assert line == 'wow\n'


def test_file_real():
    # This verifies that:
    ## context manager works
    ## writing to file works
    afile = File('./tempfile', 'w')
    with afile.open() as f:
        f.write('wow')
    afile = File('./tempfile')
    lines = afile.readlines()
    assert lines == ['wow']