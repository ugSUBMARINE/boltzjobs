"""Tests for version functionality."""

import pytest
from boltzjobs import __version__


class TestVersion:
    """Test version information access."""

    @pytest.mark.unit
    def test_version_exists(self):
        """Test that __version__ attribute exists."""
        assert __version__ is not None
        assert isinstance(__version__, str)

    @pytest.mark.unit
    def test_version_format(self):
        """Test that version follows semantic versioning format."""
        # Should be in format X.Y.Z or X.Y.Z-dev
        import re

        pattern = r"^\d+\.\d+\.\d+(-dev)?$"
        assert re.match(pattern, __version__), (
            f"Version '{__version__}' doesn't match expected format"
        )

    @pytest.mark.unit
    def test_version_accessible_from_module(self):
        """Test that version can be accessed from the module."""
        import boltzjobs

        assert hasattr(boltzjobs, "__version__")
        assert boltzjobs.__version__ == __version__

    @pytest.mark.unit
    def test_version_in_all(self):
        """Test that __version__ is properly exported."""
        import boltzjobs

        assert "__version__" in boltzjobs.__all__
