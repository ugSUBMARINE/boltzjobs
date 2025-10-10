from .jobs import Job

try:
    from importlib.metadata import version, PackageNotFoundError
except ImportError:
    # Python < 3.8 fallback
    from importlib_metadata import version, PackageNotFoundError  # type: ignore

try:
    __version__ = version("boltzjobs")
except PackageNotFoundError:
    # Package not installed, use development version
    __version__ = "0.2.0-dev"

__all__ = ["Job", "__version__"]
