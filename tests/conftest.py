"""Module for pytest config tools."""
import pytest
from ga4gh.vrsatile.pydantic.vrs_models import SequenceInterval, \
    CytobandInterval, SequenceLocation, DerivedSequenceExpression, Number, \
    IndefiniteRange, DefiniteRange, Allele, LiteralSequenceExpression, Gene, \
    ChromosomeLocation
from ga4gh.vrsatile.pydantic.vrsatile_models import Extension, Expression, \
    SequenceDescriptor, LocationDescriptor, GeneDescriptor, VCFRecord


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
def chromosome_location(cytoband_interval):
    """Create test fixture for Chromosome Location."""
    return ChromosomeLocation(
        chr="19",
        interval=cytoband_interval,
        species_id="taxonomy:9606"
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


@pytest.fixture(scope="session")
def extension():
    """Create test fixture for Extension."""
    return Extension(name="name", value=["value1", "value2"])


@pytest.fixture(scope="session")
def expression():
    """Create test fixture for Expression."""
    return Expression(syntax="hgvs.p", value="NP_005219.2:p.Leu858Arg",
                      version="1.0")


@pytest.fixture(scope="session")
def sequence_descriptor():
    """Create test fixture for Sequence Descriptor."""
    return SequenceDescriptor(id="vod:id", sequence_id="sequence:id")


@pytest.fixture(scope="session")
def location_descriptor(chromosome_location):
    """Create test fixture for Location Descriptor."""
    return LocationDescriptor(id="vod:id", location_id="gene:a",
                              location=chromosome_location)


@pytest.fixture(scope="session")
def gene_descriptor(gene):
    """Create test fixture for Gene Descriptor."""
    return GeneDescriptor(id="vod:id", gene_id="gene:abl1")


@pytest.fixture(scope="session")
def vcf_record():
    """Create test fixture for VCF Record."""
    return VCFRecord(genome_assembly="grch38", chrom="9", pos=123,
                     ref="A", alt="C")


@pytest.fixture(scope="session")
def braf_v600e_variation():
    """Create test fixture for BRAF V600E variation"""
    return {
        "id": "ga4gh:VA.8JkgnqIgYqufNl-OV_hpRG_aWF9UFQCE",
        "type": "Allele",
        "location": {
            "id": "ga4gh:VSL.AqrQ-EkAvTrXOFn70_8i3dXF5shBBZ5i",
            "type": "SequenceLocation",
            "sequence_id": "ga4gh:SQ.WaAJ_cXXn9YpMNfhcq9lnzIvaB9ALawo",
            "interval": {
                "type": "SequenceInterval",
                "start": {"type": "Number", "value": 639},
                "end": {"type": "Number", "value": 640}
            }
        },
        "state": {"type": "LiteralSequenceExpression", "sequence": "E"}
    }


@pytest.fixture(scope="session")
def braf_v600e_vd(braf_v600e_variation):
    """Create test fixture for BRAF V600E variation descriptor"""
    return {
        "id": "normalize.variation:braf%20v600e",
        "type": "VariationDescriptor",
        "variation_id": "ga4gh:VA.8JkgnqIgYqufNl-OV_hpRG_aWF9UFQCE",
        "variation": braf_v600e_variation,
        "molecule_context": "protein",
        "structural_type": "SO:0001606",
        "vrs_ref_allele_seq": "V"
    }
