"""Command line interface for boltzjobs."""

import sys
from . import __version__


def main() -> None:
    """Main CLI entry point."""
    if len(sys.argv) > 1 and sys.argv[1] in ("--version", "-V"):
        print(f"boltzjobs {__version__}")
    else:
        print(f"boltzjobs {__version__}")
        print("A Python package for creating Boltz-2 YAML input files.")
        print()
        print("Usage:")
        print("  python -m boltzjobs --version    Show version")
        print("  python -c 'import boltzjobs; print(boltzjobs.__version__)'")


if __name__ == "__main__":
    main()
