"""Module for testing the VRS model."""
import pytest

import pydantic

from ga4gh.vrsatile.pydantic.vrs_models import ComposedSequenceExpression, Genotype,\
    GenotypeMember, Number, Comparator, IndefiniteRange, DefiniteRange, Text, \
    DerivedSequenceExpression, LiteralSequenceExpression, \
    RepeatedSequenceExpression, SequenceLocation, VariationSet, Haplotype, \
    AbsoluteCopyNumber, Allele, ChromosomeLocation, Feature, SystemicVariation, \
    RelativeCopyNumber


def test_number(number):
    """Test that Number model works correctly."""
    assert number.value == 3
    assert number.type == "Number"

    assert Number(value=2, type="Number")

    invalid_params = [
        {"value": '2'},
        {"value": 2.0},
        {"value": 2, "n": 1},
        {"value": 2, "type": "number"}
    ]

    for invalid_param in invalid_params:
        with pytest.raises(pydantic.error_wrappers.ValidationError):
            Number(**invalid_param)


def test_comparator():
    """Test that Comparator model works correctly."""
    assert Comparator.LT_OR_EQUAL == "<="
    assert Comparator.GT_OR_EQUAL == ">="


def test_indefinite_range(indefinite_range):
    """Test that Indefinite Range model works correctly."""
    assert indefinite_range.value == 3
    assert indefinite_range.comparator == Comparator.GT_OR_EQUAL

    assert IndefiniteRange(value=2, comparator="<=", type="IndefiniteRange")

    invalid_params = [
        {"value": "3", "comparator": Comparator.LT_OR_EQUAL},
        {"values": 2, "comparator": Comparator.LT_OR_EQUAL},
        {"value": 3, "comparator": "=="},
        {"value": 2, "comparator": Comparator.LT_OR_EQUAL, "type": "s"}
    ]

    for invalid_param in invalid_params:
        with pytest.raises(pydantic.error_wrappers.ValidationError):
            IndefiniteRange(**invalid_param)


def test_definite_range(definite_range):
    """Test that Definite Range model works correctly."""
    assert definite_range.min == 22
    assert definite_range.max == 33
    assert definite_range.type == "DefiniteRange"

    assert DefiniteRange(min=0, max=2, type="DefiniteRange")

    invalid_params = [
        {"min": 22, "max": 33, "type": "IndefiniteRange"},
        {"min": 22.0, "max": 33}
    ]

    for invalid_param in invalid_params:
        with pytest.raises(pydantic.error_wrappers.ValidationError):
            DefiniteRange(**invalid_param)


def test_text():
    """Test that Text model works correctly."""
    definition = "APOE_LOSS"
    t = Text(id="ga4gh:id", definition=definition)
    assert t.definition == definition
    assert t.type == "Text"
    assert t.id == "ga4gh:id"

    params = {"definition": definition, "id": "ga4gh:id"}
    assert Text(**params)

    invalid_params = [
        {"definition": definition, "type": "Definition"},
        {"definition": definition, "_id": "id"}
    ]

    for invalid_param in invalid_params:
        with pytest.raises(pydantic.error_wrappers.ValidationError):
            Text(**invalid_param)


def test_literal_sequence_expression(literal_sequence_expression):
    """Test that Literal Sequence Expression model works correctly."""
    assert literal_sequence_expression.sequence == "ACGT"
    assert literal_sequence_expression.type == "LiteralSequenceExpression"

    assert LiteralSequenceExpression(sequence="ACGT",
                                     type="LiteralSequenceExpression")

    invalid_params = [
        {"sequence": "actg"},
        {"sequence": "ACTx"},
        {"sequence": "ACT", "type": "RepeatedSequenceExpression"}
    ]

    for invalid_param in invalid_params:
        with pytest.raises(pydantic.error_wrappers.ValidationError):
            LiteralSequenceExpression(**invalid_param)


def test_chromosome_location(chromosome_location):
    """Test that Chromosome Location model works correctly."""
    assert chromosome_location.chr == "19"
    assert chromosome_location.start == chromosome_location.end == "q13.32"
    assert chromosome_location.type == "ChromosomeLocation"

    assert ChromosomeLocation(
        chr="X",
        start="q13.32",
        end="q13.32",
        species_id="taxonomy:9606",
        type="ChromosomeLocation"
    )

    invalid_params = [
        {"chr": "1",
         "interval": {"start": "q13.32", "end": "q13.32"},
         "species_id": "taxonomy:9606"
         },
        {"chr": "1",
         "start": "q13.32",
         "end": "q13.32",
         "species": "taxonomy:9606"
         },
        {"chr": "23",
         "start": "q13.32",
         "end": "q13.32",
         "species_id": "taxonomy:9606",
         },
        {"chr": 1,
         "start": "q13.32",
         "end": "q13.32",
         "species_id": "taxonomy:9606",
         }
    ]

    for invalid_param in invalid_params:
        with pytest.raises(pydantic.error_wrappers.ValidationError):
            ChromosomeLocation(**invalid_param)


def test_sequence_location(sequence_location):
    """Test that Sequence Location model works correctly."""
    assert sequence_location.sequence_id == \
           "ga4gh:SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl"
    assert sequence_location.start.value == 44908821
    assert sequence_location.end.value == 44908822
    assert sequence_location.type == "SequenceLocation"

    s = SequenceLocation(id="sequence:id", sequence_id="refseq:NC_000007.13",
                         start=Number(value=44908821), end=Number(value=44908822),
                         type="SequenceLocation")
    assert s.id == "sequence:id"
    assert s.sequence_id == "refseq:NC_000007.13"
    assert sequence_location.type == "SequenceLocation"

    params = {
        "id": "sequence:id",
        "sequence_id": "refseq:NC_000007.13",
        "start": {"value": 44908821, "type": "Number"},
        "end": {"value": 44908822, "type": "Number"}
    }
    assert SequenceLocation(**params)

    params = {
        "id": "sequence:id",
        "sequence_id": "refseq:NC_000007.13",
        "start": {"value": 44908821, "type": "Number"},
        "end": {"value": 44908822, "type": "Number"}
    }
    assert SequenceLocation(**params)

    invalid_params = [
        {
            "_id": "sequence",
            "sequence_id": "NC_000007.13",
            "start": {"value": 44908821, "type": "Number"},
            "end": {"value": 44908822, "type": "Number"}
        },
        {
            "id": "sequence:1",
            "sequence_id": "NC_000007.13",
            "start": {"value": 44908821, "type": "Number"},
            "end": {"value": 44908822, "type": "Number"},
            "type": "ChromosomeLocation"
        }
    ]

    for invalid_param in invalid_params:
        with pytest.raises(pydantic.error_wrappers.ValidationError):
            SequenceLocation(**invalid_param)


def test_derived_sequence_expression(sequence_location,
                                     derived_sequence_expression):
    """Test that Derived Sequence Expression model works correctly."""
    assert derived_sequence_expression.reverse_complement is False

    d = DerivedSequenceExpression(
        location=sequence_location,
        reverse_complement=True,
    )
    assert d.reverse_complement is True
    assert d.type == "DerivedSequenceExpression"

    assert DerivedSequenceExpression(
        location=sequence_location,
        reverse_complement=False,
        type="DerivedSequenceExpression"
    )

    invalid_params = [
        {"location": sequence_location, "reverse_complement": 0},
        {"location": sequence_location, "reverse_complement": False,
         "type": "DSE"}
    ]

    for invalid_param in invalid_params:
        with pytest.raises(pydantic.error_wrappers.ValidationError):
            DerivedSequenceExpression(**invalid_param)


def test_repeated_sequence_expression(repeated_sequence_expression,
                                      derived_sequence_expression,
                                      definite_range, indefinite_range):
    """Test that Repeated Sequence Expression model works correctly."""

    def _check_seq_expr(r):
        """Test that seq_expr has intended values."""
        assert r.seq_expr.location.sequence_id == \
               "ga4gh:SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl"
        assert r.seq_expr.location.start.value == 44908821
        assert r.seq_expr.location.end.value == 44908822
        assert r.seq_expr.reverse_complement is False

    _check_seq_expr(repeated_sequence_expression)
    assert repeated_sequence_expression.count.value == 3
    assert repeated_sequence_expression.count.type == "Number"
    assert repeated_sequence_expression.type == "RepeatedSequenceExpression"

    r = RepeatedSequenceExpression(
        seq_expr=derived_sequence_expression,
        count=definite_range,
        type="RepeatedSequenceExpression"
    )
    _check_seq_expr(r)
    assert r.count.min == 22
    assert r.count.max == 33
    assert r.count.type == "DefiniteRange"

    r = RepeatedSequenceExpression(
        seq_expr=derived_sequence_expression,
        count=indefinite_range
    )
    _check_seq_expr(r)
    assert r.count.value == 3
    assert r.count.comparator == Comparator.GT_OR_EQUAL
    assert r.count.type == "IndefiniteRange"

    invalid_params = [
        {"type": "RepeatedSequenceExpression"},
        {"seq_expr": derived_sequence_expression, "count": 2},
        {"seq_expr": derived_sequence_expression, "count": definite_range,
         "type": "RSE"}
    ]

    for invalid_param in invalid_params:
        with pytest.raises(pydantic.error_wrappers.ValidationError):
            RepeatedSequenceExpression(**invalid_param)


def test_composed_sequence_expression(derived_sequence_expression,
                                      repeated_sequence_expression,
                                      literal_sequence_expression):
    """Test that ComposedSequenceExpression model works correctly"""
    def _composed_seq_expr_checks(composed_seq_expr, comp):
        assert composed_seq_expr.components == comp
        assert composed_seq_expr.type == "ComposedSequenceExpression"

    components = [literal_sequence_expression, derived_sequence_expression,
                  repeated_sequence_expression]
    cse = ComposedSequenceExpression(components=components)
    _composed_seq_expr_checks(cse, components)

    cse = ComposedSequenceExpression(**{"components": components})
    _composed_seq_expr_checks(cse, components)

    invalid_params = [
        {"type": "ComposedSequenceExpression",
         "components": [literal_sequence_expression]},
        {"type": "ComposedSequenceExpression",
         "components": [literal_sequence_expression, literal_sequence_expression]}
    ]
    for invalid_param in invalid_params:
        with pytest.raises(pydantic.error_wrappers.ValidationError):
            ComposedSequenceExpression(**invalid_param)


def test_allele(allele, sequence_location, derived_sequence_expression):
    """Test that Allele model works correctly."""
    assert allele.type == "Allele"
    assert allele.location.type == "SequenceLocation"
    assert allele.location.sequence_id == \
           "ga4gh:SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl"

    a = Allele(location="ga4gh:location", state=derived_sequence_expression,
               type="Allele")
    assert a.location == "ga4gh:location"

    invalid_params = [
        {"location": sequence_location, "state": sequence_location},
        {"location": "loc:1", "state": derived_sequence_expression,
         "type": "Haplotype"}
    ]

    for invalid_param in invalid_params:
        with pytest.raises(pydantic.error_wrappers.ValidationError):
            Allele(**invalid_param)


def test_haplotype(allele, haplotype):
    """Test that Haplotype model works correctly."""
    assert len(haplotype.members) == 2

    invalid_params = [
        {"members": [allele], "type": "Allele"},
        {"members": allele},
        {"members": [allele, "ga4ghVA"]},
        {"members": [allele], "id_": "ga4gh:VA"}
    ]

    for invalid_param in invalid_params:
        with pytest.raises(pydantic.error_wrappers.ValidationError):
            Haplotype(**invalid_param)


def test_absolute_copy_number(number, definite_range, gene, allele, sequence_location,
                              chromosome_location):
    """Test that Absolute Copy Number model works correctly."""
    c = AbsoluteCopyNumber(location="location:curie", copies=number)
    assert c.location == "location:curie"
    assert c.copies.value == 3
    assert c.type == "AbsoluteCopyNumber"

    c = AbsoluteCopyNumber(location=sequence_location, copies=number)
    assert c.location.type == "SequenceLocation"

    c = AbsoluteCopyNumber(location=chromosome_location, copies=number)
    assert c.location.type == "ChromosomeLocation"

    invalid_params = [
        {"locations": number, "copies": number},
        {"ID": "ga4gh:id", "location": gene, "copies": number},
        {"location": [allele], "copies": number},
        {"location": gene, "copies": definite_range}
    ]

    for invalid_param in invalid_params:
        with pytest.raises(pydantic.error_wrappers.ValidationError):
            AbsoluteCopyNumber(**invalid_param)


def test_relative_copy_number(number, sequence_location, gene, allele,
                              chromosome_location):
    """Test that Relative Copy Number model works correctly."""
    c = RelativeCopyNumber(location="location:curie",
                           relative_copy_class="EFO:0030069")
    assert c.location == "location:curie"
    assert c.relative_copy_class == "EFO:0030069"
    assert c.type == "RelativeCopyNumber"

    c = RelativeCopyNumber(location=sequence_location,
                           relative_copy_class="EFO:0030073")
    assert c.location.type == "SequenceLocation"
    assert c.relative_copy_class == "EFO:0030073"

    c = RelativeCopyNumber(location=chromosome_location,
                           relative_copy_class="EFO:0030068")
    assert c.location.type == "ChromosomeLocation"
    assert c.relative_copy_class == "EFO:0030068"

    for relative_copy_class in {"EFO:0030070", "EFO:0030072", "EFO:0030067",
                                "EFO:0030069"}:
        assert RelativeCopyNumber(**{
            "location": sequence_location,
            "relative_copy_class": relative_copy_class
        })

    invalid_params = [
        {"location": gene, "copies": number, "relative_copy_class": "EFO:0030068"},
        {"location": number, "copies": number},
        {"ID": "ga4gh:id", "location": gene, "relative_copy_class": "EFO:0030068"},
        {"location": allele},
        {"location": "fake:curie", "relative_copy_class": "EFO:0030068", "extra": 0},
        {"location": allele, "relative_copy_class": "partial loss"},
        {"location": sequence_location, "relative_copy_class": "complete loss"},
        {"location": sequence_location, "relative_copy_class": "fake:curie"}
    ]

    for invalid_param in invalid_params:
        with pytest.raises(pydantic.error_wrappers.ValidationError):
            RelativeCopyNumber(**invalid_param)


def test_variation_set(allele, sequence_location):
    """Test that Variation Set model works correctly."""
    v = VariationSet(members=[])
    assert len(v.members) == 0
    assert v.type == "VariationSet"

    v = VariationSet(members=[allele], type="VariationSet")
    assert len(v.members) == 1
    assert v.members[0].type == "Allele"
    assert v.members[0].location.type == "SequenceLocation"
    assert v.members[0].location.sequence_id == \
           "ga4gh:SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl"

    invalid_params = [
        {"members": [1]},
        {"members": [allele, sequence_location]},
        {"members": [allele], "type": "VS"}
    ]

    for invalid_param in invalid_params:
        with pytest.raises(pydantic.error_wrappers.ValidationError):
            VariationSet(**invalid_param)


def test_feature(gene):
    """Test Feature class."""
    schema = Feature.schema()
    assert schema["title"] == "Feature"
    assert schema["description"] == "A named entity that can be mapped to a Location. Genes, protein domains,\nexons, and chromosomes are some examples of common biological entities\nthat may be Features."  # noqa: E501
    assert schema
    assert schema["anyOf"][0]["$ref"] == "#/components/schemas/Gene"

    assert Feature(__root__=gene)


def test_systemic_variation(sequence_location, number):
    """Test SystemicVariation class."""
    c = AbsoluteCopyNumber(location=sequence_location, copies=number)
    assert SystemicVariation(__root__=c)


def test_genotype_member(indefinite_range, allele, haplotype, definite_range):
    """Test GenotypeMember model works correctly"""
    genotype_member = GenotypeMember(count=indefinite_range, variation=allele)
    assert genotype_member.count == indefinite_range
    assert genotype_member.variation == allele
    assert genotype_member.type == "GenotypeMember"

    genotype_member = GenotypeMember(**{"count": definite_range,
                                        "variation": haplotype})
    assert genotype_member.count == definite_range
    assert genotype_member.variation == haplotype
    assert genotype_member.type == "GenotypeMember"

    invalid_params = [
        {"count": 1, "variation": allele},
        {"count": definite_range, "variation": "vrs:digest"}
    ]
    for invalid_param in invalid_params:
        with pytest.raises(pydantic.error_wrappers.ValidationError):
            GenotypeMember(**invalid_param)


def test_genotype(number, allele, haplotype, definite_range, indefinite_range):
    """Test Genotype model works correctly"""
    genotype_member1 = GenotypeMember(count=number, variation=allele)
    genotype = Genotype(count=definite_range, members=[genotype_member1])
    assert genotype.count == definite_range
    assert len(genotype.members) == 1
    assert genotype.members[0] == genotype_member1
    assert genotype.type == "Genotype"

    genotype_member2 = GenotypeMember(count=definite_range, variation=haplotype)
    genotype = Genotype(**{"count": indefinite_range,
                           "members": [genotype_member1, genotype_member2]})
    assert genotype.count == indefinite_range
    assert len(genotype.members) == 2
    assert genotype.members[0] == genotype_member1
    assert genotype.members[1] == genotype_member2
    assert genotype.type == "Genotype"

    invalid_params = [
        {"count": number, "members": [genotype_member1, genotype_member1]},
        {"count": 1, "members": [genotype_member1]}
    ]
    for invalid_param in invalid_params:
        with pytest.raises(pydantic.error_wrappers.ValidationError):
            Genotype(**invalid_param)
