"""Module for pytest config tools."""
import pytest
from ga4gh.vrsatile.pydantic.core_models import CombinationTherapeuticCollection, \
    Disease, Extension, Gene, Phenotype, SubstituteTherapeuticCollection, Therapeutic
from ga4gh.vrsatile.pydantic.vrs_models import Haplotype, RepeatedSequenceExpression,\
    SequenceLocation, DerivedSequenceExpression, Number, IndefiniteRange,\
    DefiniteRange, Allele, LiteralSequenceExpression, ChromosomeLocation
from ga4gh.vrsatile.pydantic.vrsatile_models import Expression, \
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
def chromosome_location():
    """Create test fixture for Chromosome Location."""
    return ChromosomeLocation(
        chr="19",
        start="q13.32",
        end="q13.32",
        species_id="taxonomy:9606"
    )


@pytest.fixture(scope="session")
def sequence_location():
    """Create test fixture for Sequence Location."""
    return SequenceLocation(
        sequence_id="ga4gh:SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl",
        start=Number(value=44908821),
        end=Number(value=44908822)
    )


@pytest.fixture(scope="session")
def literal_sequence_expression():
    """Create test fixture for Literal Sequence Expression"""
    return LiteralSequenceExpression(sequence="ACGT")


@pytest.fixture(scope="session")
def derived_sequence_expression(sequence_location):
    """Create test fixture for Derived Sequence Expression.."""
    return DerivedSequenceExpression(
        location=sequence_location,
        reverse_complement=False,
    )


@pytest.fixture(scope="session")
def repeated_sequence_expression(derived_sequence_expression, number):
    """Create test fixture for Repeated Sequence Expression"""
    return RepeatedSequenceExpression(
        seq_expr=derived_sequence_expression, count=number
    )


@pytest.fixture(scope="session")
def allele(sequence_location):
    """Create test fixture for Allele."""
    return Allele(
        location=sequence_location,
        state=LiteralSequenceExpression(sequence="C")
    )


@pytest.fixture(scope="session")
def haplotype():
    """Create test fixture for Haplotype"""
    return Haplotype(members=["ga4gh:VA.-kUJh47Pu24Y3Wdsk1rXEDKsXWNY-68x",
                              "ga4gh:VA.Z_rYRxpUvwqCLsCBO3YLl70o2uf9_Op1"])


@pytest.fixture(scope="session")
def gene():
    """Create test fixture for Gene."""
    return Gene(id="ncbigene:348")


@pytest.fixture(scope="module")
def phenotype():
    """Create test fixture for Phenotype"""
    return Phenotype(id="HP:0000002", type="Phenotype")


@pytest.fixture(scope="module")
def disease():
    """Create test fixture for Disease"""
    return Disease(id="ncit:C4989", type="Disease")


@pytest.fixture(scope="session")
def therapeutic1():
    """Create test fixture for therapeutic"""
    return Therapeutic(id="rxcui:282388")


@pytest.fixture(scope="session")
def therapeutic2():
    """Create test fixture for therapeutic"""
    return Therapeutic(id="rxcui:1147220")


@pytest.fixture(scope="session")
def combination_therapeutic_collection(therapeutic1, therapeutic2):
    """Create test fixture for combination therapeutic collection"""
    return CombinationTherapeuticCollection(members=[therapeutic1, therapeutic2])


@pytest.fixture(scope="session")
def substitute_therapeutic_collection(therapeutic1, therapeutic2):
    """Create test fixture for substitute therapeutic collection"""
    return SubstituteTherapeuticCollection(members=[therapeutic1, therapeutic2])


@pytest.fixture(scope="session")
def extension():
    """Create test fixture for Extension."""
    return Extension(name="name", value=["value1", "value2"])


@pytest.fixture(scope="session")
def expression():
    """Create test fixture for Expression."""
    return Expression(syntax="hgvs.p", value="NP_005219.2:p.Leu858Arg",
                      syntax_version="1.0")


@pytest.fixture(scope="session")
def sequence_descriptor():
    """Create test fixture for Sequence Descriptor."""
    return SequenceDescriptor(id="vod:id", sequence="sequence:id")


@pytest.fixture(scope="session")
def location_descriptor(chromosome_location):
    """Create test fixture for Location Descriptor."""
    return LocationDescriptor(id="vod:id", location=chromosome_location)


@pytest.fixture(scope="session")
def gene_descriptor(gene):
    """Create test fixture for Gene Descriptor."""
    return GeneDescriptor(id="vod:id", gene="gene:abl1")


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
            "start": {"type": "Number", "value": 639},
            "end": {"type": "Number", "value": 640}
        },
        "state": {"type": "LiteralSequenceExpression", "sequence": "E"}
    }


@pytest.fixture(scope="session")
def braf_v600e_vd(braf_v600e_variation):
    """Create test fixture for BRAF V600E variation descriptor"""
    return {
        "id": "normalize.variation:braf%20v600e",
        "type": "VariationDescriptor",
        "variation": braf_v600e_variation,
        "molecule_context": "protein",
        "structural_type": "SO:0001606",
        "vrs_ref_allele_seq": "V"
    }
