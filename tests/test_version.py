"""Module to check version is a valid semantic version

https://semver.org/#is-there-a-suggested-regular-expression-regex-to-check-a-semver-string
"""
import re

from ga4gh.vrsatile.pydantic.version import __version__


def test_version():
    """Test that version is a semantic version"""
    regex = r"^(?P<major>0|[1-9]\d*)\.(?P<minor>0|[1-9]\d*)\.(?P<patch>0|[1-9]\d*)(?:-(?P<prerelease>(?:0|[1-9]\d*|\d*[a-zA-Z-][0-9a-zA-Z-]*)(?:\.(?:0|[1-9]\d*|\d*[a-zA-Z-][0-9a-zA-Z-]*))*))?(?:\+(?P<buildmetadata>[0-9a-zA-Z-]+(?:\.[0-9a-zA-Z-]+)*))?$"  # noqa: E501
    semvar_re = re.compile(regex)
    assert semvar_re.match(__version__), "Not a valid semantic version"
