"""Module for testing the VRSATILE model."""
import pydantic
import pytest

from ga4gh.vrsatile.pydantic.vrs_models import Number, SequenceLocation, \
    SequenceInterval, LiteralSequenceExpression, Allele, Sequence
from ga4gh.vrsatile.pydantic.vrsatile_models import  \
    MoleculeContext, Extension, Expression, ValueObjectDescriptor, \
    SequenceDescriptor, LocationDescriptor, GeneDescriptor, \
    VariationDescriptor, VCFRecord, CanonicalVariation, ComplexVariation, \
    ComplexVariationOperator, CategoricalVariationDescriptor, VariationMember


@pytest.fixture(scope="module")
def variation_member():
    """Provide example of an individual VariationMember value."""
    return {
        "type": "VariationMember",
        "expressions": [
            {
                "type": "Expression",
                "syntax": "hgvs.g",
                "syntax_version": "1.5",
                "value": "NC_000013.10:g.20763488del",
            },
            {
                "type": "Expression",
                "syntax": "hgvs.c",
                "value": "LRG_1350t1:c.235del"
            }
        ],
        "variation_id": "ga4gh:VA.KGopzor-bEw8Ot5sAQQ5o5SVx4o7TuLN"
    }


@pytest.fixture(scope="module")
def simple_repeating_del(variation_member):
    """Provide example of a complete Categorical Variation Descriptor."""
    return {
        "id": "clinvar:17014",
        "version": "2022-02-10",
        "type": "CategoricalVariationDescriptor",
        "categorical_variation": {
            "_id": "unk:?",
            "complement": False,
            "type": "CanonicalVariation",
            "variation": {
                "_id": "ga4gh:VA.KGopzor-bEw8Ot5sAQQ5o5SVx4o7TuLN",
                "type": "Allele",
                "location": {
                    "_id": "ga4gh:VSL.p4e9kMEY9PrKZ1BbNRuFr6n30DkwXWlX",
                    "type": "SequenceLocation",
                    "sequence_id": "ga4gh:SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT",
                    "interval": {
                        "type": "SequenceInterval",
                        "start": {
                            "type": "Number",
                            "value": 20189346
                        },
                        "end": {
                            "type": "Number",
                            "value": 20189349
                        }
                    }
                },
                "state": {
                    "type": "LiteralSequenceExpression",
                    "sequence": "GG"
                }
            }
        },
        "members": [variation_member]
    }


def test_molecule_context():
    """Test that Molecule Context model works correctly."""
    assert [key for key in MoleculeContext.__members__.keys()] == \
           ["GENOMIC", "TRANSCRIPT", "PROTEIN"]
    assert MoleculeContext.GENOMIC == 'genomic'
    assert MoleculeContext.TRANSCRIPT == "transcript"
    assert MoleculeContext.PROTEIN == "protein"


def test_extension(extension):
    """Test that Extension model works correctly."""
    assert extension.name == "name"
    assert len(extension.value) == 2
    assert extension.value[0] == "value1"
    assert extension.value[1] == "value2"
    assert extension.type == "Extension"

    e = Extension(name="example", value="value", type="Extension")
    assert e.name == "example"
    assert e.value == "value"
    assert e.type == "Extension"

    invalid_params = [
        {"name": 1, "value": extension.value},
        {"name": "example", "value": extension.value, "type": "Expression"}
    ]

    for invalid_param in invalid_params:
        with pytest.raises(pydantic.error_wrappers.ValidationError):
            Extension(**invalid_param)


def test_expression(expression):
    """Test that Expression model works correctly."""
    assert expression.syntax == "hgvs.p"
    assert expression.value == "NP_005219.2:p.Leu858Arg"
    assert expression.syntax_version == "1.0"
    assert expression.type == "Expression"

    e = Expression(syntax="hgvs.g", value="NC_000007.13:g.55259515T>G",
                   type="Expression")
    assert e.syntax == "hgvs.g"
    assert e.value == "NC_000007.13:g.55259515T>G"
    assert e.type == "Expression"

    invalid_params = [
        {"syntax": "valid:curie", "value": expression.value, "type": "S"},
        {"syntax": "curie", "value": expression.value},
        {"syntax": expression.syntax},
        {"value": expression.value},
        {"val": expression.value},
        {"syntax": expression.syntax, "value": 1},
        {"syntax": expression.syntax, "value": expression.value,
         "syntax_version": 1.0}
    ]

    for invalid_param in invalid_params:
        with pytest.raises(pydantic.error_wrappers.ValidationError):
            Expression(**invalid_param)


def test_value_object_descriptor(extension, expression):
    """Test that Value Object Descriptor model works correctly."""
    vod = ValueObjectDescriptor(id="value:id", type="VariationDescriptor")
    assert vod.id == "value:id"
    assert vod.type == "VariationDescriptor"

    vod = ValueObjectDescriptor(id="value:id", type="VariationDescriptor",
                                label="label", description="description",
                                xrefs=["hgnc:4"], alternate_labels=["a", "b"],
                                extensions=[extension])
    assert vod.id == "value:id"
    assert vod.type == "VariationDescriptor"
    assert vod.label == "label"
    assert vod.description == "description"
    assert vod.xrefs == ["hgnc:4"]
    assert vod.alternate_labels == ["a", "b"]
    assert vod.extensions == [extension]

    assert ValueObjectDescriptor(id="normalize.disease:id",
                                 type="DiseaseDescriptor",
                                 disease_id="ncit:C53",
                                 label="Disease label")

    invalid_params = [
        {"id": "vod:", "type": "GeneDescriptor"},
        {"id": "vod:1", "type": "SequenceDescriptor", "label": [1]},
        {"id": vod.id, "type": vod.type, "xrefs": ["xref", "xrefs"]},
        {"id": vod.id, "type": vod.type, "alternate_labels": ["xref", 1]},
        {"id": vod.id, "type": vod.type, "extensions": [extension, expression]}
    ]

    for invalid_param in invalid_params:
        with pytest.raises(pydantic.error_wrappers.ValidationError):
            ValueObjectDescriptor(**invalid_param)


def test_sequence_descriptor(sequence_descriptor, gene):
    """Test that Sequence Descriptor model works correctly."""
    assert sequence_descriptor.id == "vod:id"
    assert sequence_descriptor.sequence_id == "sequence:id"
    assert sequence_descriptor.type == "SequenceDescriptor"

    s = SequenceDescriptor(id=sequence_descriptor.id, sequence="AC",
                           type="SequenceDescriptor")
    assert s.id == "vod:id"
    assert s.sequence == "AC"
    assert s.type == "SequenceDescriptor"

    invalid_params = [
        {"id": "s:1", "sequence_id": sequence_descriptor.sequence_id,
         "type": "VariationDescriptor"},
        {"id": "sequence", "sequence_id": sequence_descriptor.sequence_id},
        {"id": sequence_descriptor.id},
        {"id": sequence_descriptor.id, "sequence": gene}
    ]

    for invalid_param in invalid_params:
        with pytest.raises(pydantic.error_wrappers.ValidationError):
            SequenceDescriptor(**invalid_param)


def test_location_descriptor(location_descriptor, sequence_location, gene):
    """Test that Location Descriptor model works correctly."""
    assert location_descriptor.id == "vod:id"
    assert location_descriptor.type == "LocationDescriptor"
    assert location_descriptor.location_id == "gene:a"
    assert location_descriptor.location.chr == "19"
    assert location_descriptor.location.interval.start == "q13.32"
    assert location_descriptor.location.interval.end == "q13.32"
    assert location_descriptor.location.species_id == "taxonomy:9606"
    assert location_descriptor.location.type == "ChromosomeLocation"

    ld = LocationDescriptor(id="vod:id2", location_id="gene:b",
                            location=sequence_location,
                            type="LocationDescriptor")
    assert ld.id == "vod:id2"
    assert ld.location_id == "gene:b"
    assert ld.location.sequence_id == \
           "ga4gh:SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl"
    assert ld.location.interval.type == "SequenceInterval"
    assert ld.location.interval.start.type == "Number" ==\
           ld.location.interval.end.type
    assert ld.location.interval.start.value == 44908821
    assert ld.location.interval.end.value == 44908822
    assert ld.location.type == "SequenceLocation"

    invalid_params = [
        {"id": "s:1", "location_id": location_descriptor.location_id,
         "type": "SequenceDescriptor"},
        {"id": "sequence", "location_id": location_descriptor.location_id},
        {"id": location_descriptor.id},
        {"id": location_descriptor.id, "sequence": gene}
    ]

    for invalid_param in invalid_params:
        with pytest.raises(pydantic.error_wrappers.ValidationError):
            LocationDescriptor(**invalid_param)


def test_gene_descriptor(gene_descriptor, gene):
    """Test that Gene Descriptor model works correctly."""
    assert gene_descriptor.id == "vod:id"
    assert gene_descriptor.gene_id == "gene:abl1"
    assert gene_descriptor.type == "GeneDescriptor"

    g = GeneDescriptor(id="vod:gene", gene=gene, gene_id="gene:348",
                       type="GeneDescriptor")
    assert g.id == "vod:gene"
    assert g.gene.gene_id == "ncbigene:348"
    assert g.gene.type == "Gene"
    assert g.gene_id == "gene:348"

    invalid_params = [
        {"id": "g:1", "gene_id": gene_descriptor.gene_id, "type": "Gene"},
        {"id": "sequence", "gene_id": gene_descriptor.gene_id},
        {"id": gene_descriptor.id},
        {"id": gene_descriptor.id, "gene": "BRAF"}
    ]

    for invalid_param in invalid_params:
        with pytest.raises(pydantic.error_wrappers.ValidationError):
            GeneDescriptor(**invalid_param)


def test_vcf_record(vcf_record):
    """Test that VCF Record model works correctly."""
    assert vcf_record.genome_assembly == "grch38"
    assert vcf_record.chrom == "9"
    assert vcf_record.pos == 123
    assert vcf_record.ref == "A"
    assert vcf_record.alt == "C"

    invalid_params = [
        {"ix": 1},
        {"genome_assembly": "grch38", "chrom": 9, "pos": 1,
         "ref": "A", "alt": "C"},
        {"genome_assembly": "grch38", "chrom": "9", "pos": 1,
         "ref": "A", "alt": "C", "qual": ["s"]},
    ]

    for invalid_param in invalid_params:
        with pytest.raises(pydantic.error_wrappers.ValidationError):
            VCFRecord(**invalid_param)


def test_variation_descriptor(allele, gene_descriptor, vcf_record, expression,
                              extension, braf_v600e_vd):
    """Test that Variation Descriptor model works correctly."""
    vd = VariationDescriptor(id="var:id", variation_id="variation:id")
    assert vd.id == "var:id"
    assert vd.variation_id == "variation:id"
    assert vd.type == "VariationDescriptor"

    vd = VariationDescriptor(**braf_v600e_vd)
    assert vd.variation.type == "Allele"
    assert vd.variation.location.type == "SequenceLocation"
    assert vd.variation.location.interval.type == "SequenceInterval"

    vd = VariationDescriptor(id="var:id", variation=allele,
                             type="VariationDescriptor",
                             gene_context=gene_descriptor,
                             vcf_record=vcf_record,
                             molecule_context="genomic",
                             expressions=[expression],
                             structural_type="SO:0001537",
                             vrs_ref_allele_seq="C",
                             allelic_state="GENO:00000875")
    assert vd.variation.location.type == "SequenceLocation"
    assert vd.variation.location.sequence_id == \
           "ga4gh:SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl"
    assert vd.variation.location.interval.type == "SequenceInterval"
    assert vd.variation.location.interval.start.type == "Number" ==\
           vd.variation.location.interval.end.type
    assert vd.variation.location.interval.start.value == 44908821
    assert vd.variation.location.interval.end.value == 44908822
    assert vd.variation.location.type == "SequenceLocation"
    assert vd.variation.state.sequence == "C"
    assert vd.gene_context.id == "vod:id"
    assert vd.gene_context.gene_id == "gene:abl1"
    assert vd.vcf_record.genome_assembly == "grch38"
    assert vd.vcf_record.chrom == "9"
    assert vd.vcf_record.pos == 123
    assert vd.vcf_record.ref == "A"
    assert vd.vcf_record.alt == "C"
    assert vd.molecule_context == "genomic"
    assert len(vd.expressions) == 1
    assert vd.expressions[0].syntax == "hgvs.p"
    assert vd.expressions[0].value == "NP_005219.2:p.Leu858Arg"
    assert vd.expressions[0].syntax_version == "1.0"
    assert vd.structural_type == "SO:0001537"
    assert vd.vrs_ref_allele_seq == "C"
    assert vd.allelic_state == "GENO:00000875"

    invalid_params = [
        {"id": "vod:id", "variation_id": "var:id", "type": "GeneDescriptor"},
        {"id": "vod:id", "variation_id": "var:id", "molecule_context": "g"},
        {"id": "vod:id", "variation_id": "var:id", "structural_type": "g"},
        {"id": "vod:id", "variation_id": "var:id", "expressions": [extension]},
        {"id": "vod:id", "variation_id": "var:id", "expressions": expression},
        {"id": "vod:id", "variation_id": "var:id", "vcf_record": expression},
        {"id": "vod:id", "variation_id": "var:id", "gene_context": extension},
        {"id": "vod:id", "variation_id": "var:id", "vrs_ref_allele_seq": "A!"},
        {"id": "vod:id", "variation_id": "var:id", "allelic_state": "ACT"},
    ]

    for invalid_param in invalid_params:
        with pytest.raises(pydantic.error_wrappers.ValidationError):
            VariationDescriptor(**invalid_param)


def test_canonical_variation(allele):
    """Test creation and usage of canonical variations."""
    cv = CanonicalVariation(
        _id="clinvar:13961",
        complement=True,
        variation=allele
    )
    assert cv.complement is True

    assert cv.variation.type == "Allele"
    assert cv.variation.location.type == "SequenceLocation"
    assert cv.variation.location.sequence_id == \
        "ga4gh:SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl"

    assert cv.variation.state.type.value == "LiteralSequenceExpression"
    assert cv.variation.state.sequence == "C"

    assert set(cv.schema()["properties"].keys()) == {"_id", "type",
                                                     "complement", "variation"}

    # check forbid extra
    with pytest.raises(pydantic.error_wrappers.ValidationError):
        cv = CanonicalVariation(
            _id="clinvar:13961",
            complement=False,
            variation=allele,
            variation_id="clinvar:13961"
        )


def test_categorical_variation(allele):
    """Test creation and usage of categorical variations."""
    cv = ComplexVariation(
        _id="complex:variation",
        operands=[
            CanonicalVariation(
                _id="clinvar:13961",
                complement=False,
                variation=allele
            ),
            CanonicalVariation(
                _id="clinvar:375941",
                complement=False,
                variation=Allele(
                    location=SequenceLocation(
                        sequence_id="ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",  # noqa: E501
                        interval=SequenceInterval(
                            start=Number(value=140753336),
                            end=Number(value=140753337)
                        )
                    ),
                    state=LiteralSequenceExpression(
                        sequence=Sequence(__root__="T")
                    )
                )
            )
        ],
        operator=ComplexVariationOperator.AND,
        complement=False
    )

    assert cv.id == "complex:variation"
    assert cv.complement is False
    assert cv.operator == "AND"

    assert len(cv.operands) == 2

    op_0: CanonicalVariation = cv.operands[0]
    assert op_0.id == "clinvar:13961"
    assert op_0.complement is False
    assert op_0.variation is not None
    assert op_0.variation.type == "Allele"
    assert op_0.variation.location.type == "SequenceLocation"
    assert op_0.variation.location.sequence_id == \
        "ga4gh:SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl"
    assert op_0.variation.location.interval.start.value == 44908821
    assert op_0.variation.location.interval.end.value == 44908822
    assert op_0.variation.state.sequence == "C"

    op_1: CanonicalVariation = cv.operands[1]
    assert op_1.id == "clinvar:375941"
    assert op_1.complement is False
    assert op_1.variation is not None
    assert op_1.variation.type == "Allele"
    assert op_1.variation.location.type == "SequenceLocation"
    assert op_1.variation.location.sequence_id == \
        "ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul"
    assert op_1.variation.location.interval.start.value == 140753336
    assert op_1.variation.location.interval.end.value == 140753337
    assert op_1.variation.state.sequence == "T"

    with pytest.raises(pydantic.ValidationError):
        ComplexVariation(
            _id="complex:variation",
            operands=[
                CanonicalVariation(
                    _id="clinvar:13961",
                    complement=False,
                    variation=allele
                )
            ],
            operator=ComplexVariationOperator.OR,
            complement=False
        )


def test_variation_member(variation_member):
    """Test variation member descriptor"""
    vm = VariationMember(**variation_member)
    assert vm.type == "VariationMember"
    assert vm.variation_id == "ga4gh:VA.KGopzor-bEw8Ot5sAQQ5o5SVx4o7TuLN"
    expr_1 = next(i for i in vm.expressions if i.syntax == "hgvs.g")
    assert expr_1.type == "Expression"
    assert expr_1.value == "NC_000013.10:g.20763488del"
    expr_2 = next(i for i in vm.expressions if i.syntax == "hgvs.c")
    assert expr_2.type == "Expression"
    assert expr_2.value == "LRG_1350t1:c.235del"

    with pytest.raises(pydantic.ValidationError):
        VariationMember(**{
            "type": "VariationMember",
            "expressions": [],
            "variation_id": "ga4gh:VA.KGopzor-bEw8Ot5sAQQ5o5SVx4o7TuLN"
        })


def test_categorical_variation_descriptor(simple_repeating_del):
    """Test categorical variation descriptor"""
    cvd = CategoricalVariationDescriptor(**simple_repeating_del)
    assert cvd.id == "clinvar:17014"
    assert cvd.version == "2022-02-10"
    assert cvd.type == "CategoricalVariationDescriptor"
    assert cvd.categorical_variation.id == "unk:?"
    assert cvd.categorical_variation.complement is False
    assert cvd.categorical_variation.type == "CanonicalVariation"
    assert cvd.categorical_variation.variation.id == \
        "ga4gh:VA.KGopzor-bEw8Ot5sAQQ5o5SVx4o7TuLN"
    location = cvd.categorical_variation.variation.location
    assert location.id == "ga4gh:VSL.p4e9kMEY9PrKZ1BbNRuFr6n30DkwXWlX"
    assert location.sequence_id == "ga4gh:SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT"
    assert location.interval.type == "SequenceInterval"
    assert cvd.categorical_variation.variation.state.sequence == "GG"

    assert len(cvd.members) == 1
    member = cvd.members[0]
    assert member.variation_id == "ga4gh:VA.KGopzor-bEw8Ot5sAQQ5o5SVx4o7TuLN"
    assert len(member.expressions) == 2
    expressions = sorted(member.expressions, key=lambda e: e.value)
    assert expressions[0].type == "Expression"
    assert expressions[0].syntax == "hgvs.c"
    assert expressions[0].value == "LRG_1350t1:c.235del"
    assert expressions[1].type == "Expression"
    assert expressions[1].syntax == "hgvs.g"
    assert expressions[1].value == "NC_000013.10:g.20763488del"
