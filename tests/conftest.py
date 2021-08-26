"""Module for pytest config tools."""
import pytest
from ga4gh.models.vrs_model import SequenceInterval, \
    CytobandInterval, SequenceLocation, DerivedSequenceExpression, Number, \
    IndefiniteRange, DefiniteRange, Allele, LiteralSequenceExpression, Gene


@pytest.fixture(scope="session")
def number():
    """Create test fixture for Number."""
    return Number(value=3)


@pytest.fixture(scope="session")
def indefinite_range():
    """Create test fixture for Indefinite Range."""
    return IndefiniteRange(value=3, comparator=">=")


@pytest.fixture(scope="session")
def definite_range():
    """Create test fixture for Definite Range."""
    return DefiniteRange(min=22, max=33)


@pytest.fixture(scope="session")
def sequence_interval():
    """Create test fixture for Sequence Interval."""
    return SequenceInterval(
        start=Number(value=44908821),
        end=Number(value=44908822)
    )


@pytest.fixture(scope="session")
def cytoband_interval():
    """Create test fixture for Cytoband Interval."""
    return CytobandInterval(
        start="q13.32", end="q13.32"
    )


@pytest.fixture(scope="session")
def sequence_location(sequence_interval):
    """Create test fixture for Sequence Location."""
    return SequenceLocation(
        sequence_id="ga4gh:SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl",
        interval=sequence_interval
    )


@pytest.fixture(scope="session")
def derived_sequence_expression(sequence_location):
    """Create test fixture for Derived Sequence Expression.."""
    return DerivedSequenceExpression(
        location=sequence_location,
        reverse_complement=False,
    )


@pytest.fixture(scope="session")
def allele(sequence_location):
    """Create test fixture for Allele."""
    return Allele(
        location=sequence_location,
        state=LiteralSequenceExpression(sequence="C")
    )


@pytest.fixture(scope="session")
def gene():
    """Create test fixture for Gene."""
    return Gene(gene_id="ncbigene:348")
