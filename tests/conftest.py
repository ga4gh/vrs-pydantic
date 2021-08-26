"""Module for pytest config tools."""
import pytest
from ga4gh.models.vrs_model import SequenceInterval, Number, \
    CytobandInterval, SequenceLocation, DerivedSequenceExpression, Number, \
    IndefiniteRange, DefiniteRange, Allele, LiteralSequenceExpression, Gene


@pytest.fixture(scope="session")
def number():
    return Number(value=3)


@pytest.fixture(scope="session")
def indefinite_range():
    return IndefiniteRange(value=3, comparator=">=")


@pytest.fixture(scope="session")
def definite_range():
    return DefiniteRange(min=22, max=33)


@pytest.fixture(scope="session")
def sequence_interval():
    return SequenceInterval(
        start=Number(value=44908821),
        end=Number(value=44908822)
    )


@pytest.fixture(scope="session")
def cytoband_interval():
    return CytobandInterval(
        start="q13.32", end="q13.32"
    )


@pytest.fixture(scope="session")
def sequence_location(sequence_interval):
    return SequenceLocation(
        sequence_id="ga4gh:SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl",
        interval=sequence_interval
    )


@pytest.fixture(scope="session")
def derived_sequence_expression(sequence_location):
    return DerivedSequenceExpression(
        location=sequence_location,
        reverse_complement=False,
    )


@pytest.fixture(scope="session")
def allele(sequence_location):
    return Allele(
        location=sequence_location,
        state=LiteralSequenceExpression(sequence="C")
    )


@pytest.fixture(scope="session")
def gene():
    return Gene(gene_id="ncbigene:348")
