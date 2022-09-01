"""Module for testing the VRSATILE model."""
import pydantic
from ga4gh.vrsatile.pydantic.core_models import Condition, Disease, Extension, Phenotype
import pytest

from ga4gh.vrsatile.pydantic.vrs_models import Number, SequenceLocation, \
    LiteralSequenceExpression, Allele, Sequence
from ga4gh.vrsatile.pydantic.vrsatile_models import  \
    CanonicalVariationDescriptor, ConditionDescriptor, DiseaseDescriptor, \
    MoleculeContext, Expression, PhenotypeDescriptor, TherapeuticCollectionDescriptor, \
    TherapeuticDescriptor, ValueObjectDescriptor, SequenceDescriptor, \
    LocationDescriptor, GeneDescriptor, VariationDescriptor, VCFRecord, \
    CanonicalVariation, ComplexVariation, ComplexVariationOperator, \
    CategoricalVariationDescriptor, VariationMember


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
        "type": "CategoricalVariationDescriptor",
        "categorical_variation": {
            "id": "unk:?",
            "type": "CanonicalVariation",
            "canonical_context": {
                "id": "ga4gh:VA.KGopzor-bEw8Ot5sAQQ5o5SVx4o7TuLN",
                "type": "Allele",
                "location": {
                    "id": "ga4gh:VSL.p4e9kMEY9PrKZ1BbNRuFr6n30DkwXWlX",
                    "type": "SequenceLocation",
                    "sequence_id": "ga4gh:SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT",
                    "start": {"type": "Number", "value": 20189346},
                    "end": {"type": "Number", "value": 20189349}
                },
                "state": {
                    "type": "LiteralSequenceExpression",
                    "sequence": "GG"
                }
            }
        },
        "members": [variation_member]
    }


@pytest.fixture(scope="module")
def phenotype_descriptor():
    """Create test fixture for phenotype descriptor"""
    return PhenotypeDescriptor(
        id="phenotype:1",
        phenotype=Phenotype(id="HP:0000002"),
        label="Abnormality of body height"
    )


@pytest.fixture(scope="module")
def disease_descriptor():
    """Create test fixture for disease descriptor"""
    return DiseaseDescriptor(
        id="disease:1",
        disease=Disease(id="ncit:C2926"),
        label="Lung Non-Small Cell Carcinoma",
        alternate_labels=["NSCLC", "Non-Small Cell Carcinoma of the Lung"],
        xrefs=["mondo:0005233", "oncotree:NSCLC", "DOID:3908"],
        extensions=[Extension(name="associated_with", value=["HP:0030358"])]
    )


def test_molecule_context():
    """Test that Molecule Context model works correctly."""
    assert [key for key in MoleculeContext.__members__.keys()] == \
           ["GENOMIC", "TRANSCRIPT", "PROTEIN"]
    assert MoleculeContext.GENOMIC == 'genomic'
    assert MoleculeContext.TRANSCRIPT == "transcript"
    assert MoleculeContext.PROTEIN == "protein"


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


def test_value_object_descriptor(extension, expression, gene):
    """Test that Value Object Descriptor model works correctly."""
    vod = ValueObjectDescriptor(id="value:id", type="VariationDescriptor",
                                value="value:id")
    assert vod.id == "value:id"
    assert vod.type == "VariationDescriptor"
    assert vod.value == "value:id"

    vod = ValueObjectDescriptor(id="value:id", type="VariationDescriptor",
                                value=gene, label="label", description="description",
                                xrefs=["hgnc:4"], alternate_labels=["a", "b"],
                                extensions=[extension])
    assert vod.id == "value:id"
    assert vod.type == "VariationDescriptor"
    assert vod.label == "label"
    assert vod.description == "description"
    assert vod.xrefs == ["hgnc:4"]
    assert vod.alternate_labels == ["a", "b"]
    assert vod.extensions == [extension]
    assert vod.value == gene

    invalid_params = [
        {"id": "vod:", "type": "GeneDescriptor"},
        {"id": "vod:1", "type": "SequenceDescriptor", "label": [1]},
        {"id": vod.id, "type": vod.type, "xrefs": ["xref", "xrefs"]},
        {"id": vod.id, "type": vod.type, "alternate_labels": ["xref", 1]},
        {"id": vod.id, "type": vod.type, "extensions": [extension, expression]},
        {"id": vod.id, "type": vod.type, "value": "value:id",
         "xrefs": ["test:1", "test:1"]}
    ]

    for invalid_param in invalid_params:
        with pytest.raises(pydantic.error_wrappers.ValidationError):
            ValueObjectDescriptor(**invalid_param)


def test_sequence_descriptor(sequence_descriptor, gene):
    """Test that Sequence Descriptor model works correctly."""
    assert sequence_descriptor.id == "vod:id"
    assert sequence_descriptor.sequence == "sequence:id"
    assert sequence_descriptor.type == "SequenceDescriptor"

    s = SequenceDescriptor(id=sequence_descriptor.id, sequence="AC",
                           type="SequenceDescriptor")
    assert s.id == "vod:id"
    assert s.sequence == "AC"
    assert s.type == "SequenceDescriptor"

    invalid_params = [
        {"id": "s:1", "sequence_id": "test:1",
         "type": "VariationDescriptor"},
        {"id": sequence_descriptor.id},
        {"id": sequence_descriptor.id, "sequence": gene}
    ]

    for invalid_param in invalid_params:
        with pytest.raises(pydantic.error_wrappers.ValidationError):
            SequenceDescriptor(**invalid_param)


def test_location_descriptor(location_descriptor, gene):
    """Test that Location Descriptor model works correctly."""
    assert location_descriptor.id == "vod:id"
    assert location_descriptor.type == "LocationDescriptor"
    assert location_descriptor.location.chr == "19"
    assert location_descriptor.location.start == "q13.32"
    assert location_descriptor.location.end == "q13.32"
    assert location_descriptor.location.species_id == "taxonomy:9606"
    assert location_descriptor.location.type == "ChromosomeLocation"

    ld = LocationDescriptor(id="vod:id2", location="gene:b", type="LocationDescriptor")
    assert ld.id == "vod:id2"
    assert ld.location == "gene:b"

    invalid_params = [
        {"id": "s:1", "location_id": "location:1",
         "type": "SequenceDescriptor"},
        {"id": location_descriptor.id},
        {"id": location_descriptor.id, "sequence": gene}
    ]

    for invalid_param in invalid_params:
        with pytest.raises(pydantic.error_wrappers.ValidationError):
            LocationDescriptor(**invalid_param)


def test_gene_descriptor(gene_descriptor, gene):
    """Test that Gene Descriptor model works correctly."""
    assert gene_descriptor.id == "vod:id"
    assert gene_descriptor.gene == "gene:abl1"
    assert gene_descriptor.type == "GeneDescriptor"

    g = GeneDescriptor(id="vod:gene", gene=gene, type="GeneDescriptor")
    assert g.id == "vod:gene"
    assert g.gene.id == "ncbigene:348"
    assert g.gene.type == "Gene"

    invalid_params = [
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
    vd = VariationDescriptor(id="var:id", variation="variation:id")
    assert vd.id == "var:id"
    assert vd.variation == "variation:id"
    assert vd.type == "VariationDescriptor"

    vd = VariationDescriptor(**braf_v600e_vd)
    assert vd.variation.type == "Allele"
    assert vd.variation.location.type == "SequenceLocation"

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
    assert vd.variation.location.start.type == "Number" == \
           vd.variation.location.end.type
    assert vd.variation.location.start.value == 44908821
    assert vd.variation.location.end.value == 44908822
    assert vd.variation.location.type == "SequenceLocation"
    assert vd.variation.state.sequence == "C"
    assert vd.gene_context.id == "vod:id"
    assert vd.gene_context.gene == "gene:abl1"
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
    cv = CanonicalVariation(id="clinvar:13961", canonical_context=allele)

    assert cv.canonical_context.type == "Allele"
    assert cv.canonical_context.location.type == "SequenceLocation"
    assert cv.canonical_context.location.sequence_id == \
        "ga4gh:SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl"

    assert cv.canonical_context.state.type.value == "LiteralSequenceExpression"
    assert cv.canonical_context.state.sequence == "C"

    assert set(cv.schema()["properties"].keys()) == {"id", "type", "canonical_context"}

    # check forbid extra
    with pytest.raises(pydantic.error_wrappers.ValidationError):
        CanonicalVariation(
            id="clinvar:13961",
            complement=False,
            canonical_context=allele
        )


def test_complex_variation(allele, braf_v600e_variation):
    """Test that ComplexVariation model works correctly"""
    canonical_v1 = CanonicalVariation(id="canonical:1", canonical_context=allele)
    canonical_v2 = CanonicalVariation(id="canonical:2",
                                      canonical_context=braf_v600e_variation)
    complex_variation = ComplexVariation(
        id="complex:1",
        operands=[canonical_v1, canonical_v2],
        operator="AND"
    )
    assert complex_variation.id == "complex:1"
    assert complex_variation.operands == [canonical_v1, canonical_v2]
    assert complex_variation.operator == "AND"
    assert complex_variation.type == "ComplexVariation"

    invalid_params = [
        {"id": "value:1", "operands": [canonical_v1]},
        {"id": "value:1", "operands": [canonical_v1, canonical_v2], "operator": "NOT"}
    ]
    for invalid_param in invalid_params:
        with pytest.raises(pydantic.error_wrappers.ValidationError):
            ComplexVariation(**invalid_param)


def test_categorical_variation(allele):
    """Test creation and usage of categorical variations."""
    cv = ComplexVariation(
        id="complex:variation",
        operands=[
            CanonicalVariation(
                id="clinvar:13961",
                canonical_context=allele
            ),
            CanonicalVariation(
                id="clinvar:375941",
                canonical_context=Allele(
                    location=SequenceLocation(
                        sequence_id="ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",
                        start=Number(value=140753336),
                        end=Number(value=140753337)
                    ),
                    state=LiteralSequenceExpression(
                        sequence=Sequence(__root__="T")
                    )
                )
            )
        ],
        operator=ComplexVariationOperator.AND
    )

    assert cv.id == "complex:variation"
    assert cv.operator == "AND"

    assert len(cv.operands) == 2

    op_0: CanonicalVariation = cv.operands[0]
    assert op_0.id == "clinvar:13961"
    assert op_0.canonical_context is not None
    assert op_0.canonical_context.type == "Allele"
    assert op_0.canonical_context.location.type == "SequenceLocation"
    assert op_0.canonical_context.location.sequence_id == \
        "ga4gh:SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl"
    assert op_0.canonical_context.location.start.value == 44908821
    assert op_0.canonical_context.location.end.value == 44908822
    assert op_0.canonical_context.state.sequence == "C"

    op_1: CanonicalVariation = cv.operands[1]
    assert op_1.id == "clinvar:375941"
    assert op_1.canonical_context is not None
    assert op_1.canonical_context.type == "Allele"
    assert op_1.canonical_context.location.type == "SequenceLocation"
    assert op_1.canonical_context.location.sequence_id == \
        "ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul"
    assert op_1.canonical_context.location.start.value == 140753336
    assert op_1.canonical_context.location.end.value == 140753337
    assert op_1.canonical_context.state.sequence == "T"

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
    assert cvd.type == "CategoricalVariationDescriptor"
    assert cvd.categorical_variation.id == "unk:?"
    assert cvd.categorical_variation.type == "CanonicalVariation"
    assert cvd.categorical_variation.canonical_context.id == \
        "ga4gh:VA.KGopzor-bEw8Ot5sAQQ5o5SVx4o7TuLN"
    location = cvd.categorical_variation.canonical_context.location
    assert location.id == "ga4gh:VSL.p4e9kMEY9PrKZ1BbNRuFr6n30DkwXWlX"
    assert location.sequence_id == "ga4gh:SQ._0wi-qoDrvram155UmcSC-zA5ZK4fpLT"
    assert cvd.categorical_variation.canonical_context.state.sequence == "GG"

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


def test_disease_descriptor(disease_descriptor):
    """Test that DiseaseDescriptor model works correctly"""
    assert disease_descriptor.id == "disease:1"
    assert disease_descriptor.disease.id == "ncit:C2926"
    assert disease_descriptor.disease.type == "Disease"
    assert disease_descriptor.label == "Lung Non-Small Cell Carcinoma"
    assert disease_descriptor.alternate_labels == \
        ["NSCLC", "Non-Small Cell Carcinoma of the Lung"]
    assert disease_descriptor.xrefs == ["mondo:0005233", "oncotree:NSCLC", "DOID:3908"]
    assert len(disease_descriptor.extensions) == 1
    ext = disease_descriptor.extensions[0]
    assert ext.name == "associated_with"
    assert ext.value == ["HP:0030358"]

    disease_descr = DiseaseDescriptor(disease="ncit:C2926")
    assert disease_descr.disease == "ncit:C2926"

    with pytest.raises(pydantic.error_wrappers.ValidationError):
        DiseaseDescriptor(id="disease:1")


def test_phenotype_descriptor(phenotype_descriptor):
    """Test that PhenotypeDescriptor model works correctly"""
    assert phenotype_descriptor.id == "phenotype:1"
    assert phenotype_descriptor.phenotype.id == "HP:0000002"
    assert phenotype_descriptor.phenotype.type == "Phenotype"
    assert phenotype_descriptor.label == "Abnormality of body height"
    assert phenotype_descriptor.type == "PhenotypeDescriptor"

    phenotype_descr = PhenotypeDescriptor(phenotype="HP:0000002")
    assert phenotype_descr.phenotype == "HP:0000002"

    with pytest.raises(pydantic.error_wrappers.ValidationError):
        PhenotypeDescriptor(id="phenotype:1")


def test_condition_descriptor(phenotype_descriptor, disease_descriptor, phenotype,
                              disease):
    """Test that ConditionDescriptor model works correctly"""
    condition_descr = ConditionDescriptor(
        id="condition:1",
        condition=Condition(members=[disease, phenotype], type="Condition"),
        member_descriptors=[disease_descriptor, phenotype_descriptor],
        label="condition descriptor"
    )
    assert condition_descr.type == "ConditionDescriptor"
    assert condition_descr.condition.members == [disease, phenotype]
    assert condition_descr.condition.type == "Condition"
    assert condition_descr.member_descriptors == [disease_descriptor,
                                                  phenotype_descriptor]
    assert condition_descr.label == "condition descriptor"

    condition_descr = ConditionDescriptor(
        condition="condition:1",
        member_descriptors=[phenotype_descriptor]
    )
    assert condition_descr.condition == "condition:1"
    assert condition_descr.member_descriptors == [phenotype_descriptor]

    invalid_params = [
        {"member_descriptors": [phenotype_descriptor]},
        {"id": "condition:1", "condition_id": "condition:1"},
        {"member_descriptors": phenotype_descriptor, "condition_id": "condition:1"},
        {"condition": {"member": ["disease:1", "phenotype:1"]},
         "member_descriptors": [disease_descriptor, phenotype_descriptor]}
    ]
    for invalid_param in invalid_params:
        with pytest.raises(pydantic.error_wrappers.ValidationError):
            ConditionDescriptor(**invalid_param)


def test_canonical_variation_descriptor(allele):
    """Test that CanonicalVariationDescriptor model works correctly"""
    canonical_variation = CanonicalVariation(id="clinvar:13961",
                                             canonical_context=allele)
    subject_descriptor = VariationDescriptor(variation="clinvar:13961",
                                             label="test")
    canonical_vd = CanonicalVariationDescriptor(
        subject_variation_descriptor=subject_descriptor,
        canonical_variation=canonical_variation
    )
    assert canonical_vd.type == "CanonicalVariationDescriptor"
    assert canonical_vd.subject_variation_descriptor == subject_descriptor
    assert canonical_vd.canonical_variation == canonical_variation

    canonical_vd = CanonicalVariationDescriptor(
        subject_variation_descriptor=subject_descriptor,
        canonical_variation="canonical_variation:1"
    )
    assert canonical_vd.type == "CanonicalVariationDescriptor"
    assert canonical_vd.subject_variation_descriptor == subject_descriptor
    assert canonical_vd.canonical_variation == "canonical_variation:1"

    invalid_params = [
        {"canonical_variation_id": "curie:1"},
        {"subject_variation_descriptor": subject_descriptor}
    ]
    for invalid_param in invalid_params:
        with pytest.raises(pydantic.error_wrappers.ValidationError):
            CanonicalVariationDescriptor(**invalid_param)


def test_therapeutic_descriptor(therapeutic1, disease):
    """Test that Therapeutic Descriptor model works correctly"""
    td = TherapeuticDescriptor(therapeutic=therapeutic1, label="imatinib")
    assert td.type == "TherapeuticDescriptor"
    assert td.therapeutic.id == "rxcui:282388"
    assert td.therapeutic.type == "Therapeutic"
    assert td.label == "imatinib"

    td = TherapeuticDescriptor(**{"therapeutic": "rxcui:282388"})
    assert td.type == "TherapeuticDescriptor"
    assert td.therapeutic == "rxcui:282388"

    invalid_params = [
        {"therapy_id": "rxcui:282388"},
        {"therapeutic": disease}
    ]
    for invalid_param in invalid_params:
        with pytest.raises(pydantic.error_wrappers.ValidationError):
            TherapeuticDescriptor(**invalid_param)


def test_therapeutic_collection_descriptor(combination_therapeutic_collection,
                                           substitute_therapeutic_collection):
    """Test that Therapeutic Collection Descriptor model works correctly"""
    tcd = TherapeuticCollectionDescriptor(
        therapeutic_collection=combination_therapeutic_collection)
    assert tcd.type == "TherapeuticsCollectionDescriptor"
    assert tcd.therapeutic_collection == combination_therapeutic_collection
    assert tcd.member_descriptors is None

    tcd = TherapeuticCollectionDescriptor(**{
        "therapeutic_collection": "therapeutic_collection_descriptor:1",
        "member_descriptors": [
            TherapeuticDescriptor(therapeutic="rxcui:282388"),
            TherapeuticDescriptor(therapeutic="rxcui:1147220"),
        ]
    })
    assert tcd.type == "TherapeuticsCollectionDescriptor"
    assert tcd.therapeutic_collection == "therapeutic_collection_descriptor:1"
    assert len(tcd.member_descriptors) == 2
    assert tcd.member_descriptors[0].type == "TherapeuticDescriptor"
    assert tcd.member_descriptors[0].therapeutic == "rxcui:282388"
    assert tcd.member_descriptors[1].type == "TherapeuticDescriptor"
    assert tcd.member_descriptors[1].therapeutic == "rxcui:1147220"

    invalid_params = [
        {"therapeutic_collection_id": "invalid"}
    ]
    for invalid_param in invalid_params:
        with pytest.raises(pydantic.error_wrappers.ValidationError):
            TherapeuticCollectionDescriptor(**invalid_param)
