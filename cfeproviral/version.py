from importlib.metadata import version, PackageNotFoundError

def get_version():
    try:
        return version('cfeproviral')
    except PackageNotFoundError:
        return 'unknown'
