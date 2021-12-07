"""Module for testing the VRS model."""
import pydantic
import pytest
from ga4gh.vrsatile.pydantic.vrs_models import Number, Comparator, \
    IndefiniteRange, DefiniteRange, Text, \
    SequenceInterval, CytobandInterval, DerivedSequenceExpression, \
    LiteralSequenceExpression, RepeatedSequenceExpression, Gene, \
    SequenceLocation, VariationSet, Haplotype, \
    CopyNumber, Allele, ChromosomeLocation


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
    t = Text(definition=definition)
    assert t.definition == definition
    assert t.type == "Text"

    assert Text(definition=definition, type="Text")

    t = Text(_id="ga4gh:id", definition=definition)
    assert t.definition == definition
    assert t.type == "Text"
    assert t.id == "ga4gh:id"
    t_dict = t.dict(by_alias=True)
    assert t_dict['_id'] == 'ga4gh:id'
    assert 'id' not in t_dict.keys()

    params = {"definition": definition, "id": "ga4gh:id"}
    assert Text(**params)

    invalid_params = [
        {"definition": definition, "type": "Definition"},
        {"definition": definition, "_id": "id"}
    ]

    for invalid_param in invalid_params:
        with pytest.raises(pydantic.error_wrappers.ValidationError):
            Text(**invalid_param)


def test_sequence_interval(sequence_interval):
    """Test that Sequence Interval model works correctly."""
    assert sequence_interval.start.value == 44908821
    assert sequence_interval.end.value == 44908822
    assert sequence_interval.type == "SequenceInterval"

    assert SequenceInterval(
        start=Number(value=44908821),
        end=Number(value=44908822),
        type="SequenceInterval"
    )

    invalid_params = [
        {"start": 44908821, "end": 44908822},
        {"start": {"value": 2}, "end": {"value": 3}, "type": "M"}
    ]

    for invalid_param in invalid_params:
        with pytest.raises(pydantic.error_wrappers.ValidationError):
            SequenceInterval(**invalid_param)


def test_cytoband_interval(cytoband_interval):
    """Test that Cytoband Interval model works correctly."""
    human_cytoband = "q13.32"
    assert cytoband_interval.start == human_cytoband == cytoband_interval.end
    assert cytoband_interval.type == "CytobandInterval"

    assert CytobandInterval(
        start="q13.32", end="q13.32", type="CytobandInterval"
    )

    invalid_params = [
        {"start": "x9", "end": human_cytoband},
        {"start": human_cytoband, "end": human_cytoband, "type": "CI"}
    ]

    for invalid_param in invalid_params:
        with pytest.raises(pydantic.error_wrappers.ValidationError):
            CytobandInterval(**invalid_param)


def test_literal_sequence_expression():
    """Test that Literal Sequence Expression model works correctly."""
    lse = LiteralSequenceExpression(sequence="ACGT")
    assert lse.sequence == "ACGT"
    assert lse.type == "LiteralSequenceExpression"

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


def test_gene():
    """Test that Gene model works correctly."""
    g = Gene(gene_id="hgnc:5")
    assert g.gene_id == "hgnc:5"
    assert g.type == "Gene"

    assert Gene(gene_id="hgnc:5", type="Gene")

    invalid_params = [
        {"gene_id": "hgnc5"},
        {"gene_id": "test:1", "type": "GeneDescriptor"}
    ]

    for invalid_param in invalid_params:
        with pytest.raises(pydantic.error_wrappers.ValidationError):
            Gene(**invalid_param)


def test_chromosome_location(chromosome_location, cytoband_interval):
    """Test that Chromosome Location model works correctly."""
    assert chromosome_location.chr == "19"
    assert chromosome_location.interval.start == chromosome_location.interval.end == "q13.32"  # noqa: E501
    assert chromosome_location.type == "ChromosomeLocation"

    assert ChromosomeLocation(
        chr="19",
        interval=cytoband_interval,
        species_id="taxonomy:9606",
        type="ChromosomeLocation"
    )

    invalid_params = [
        {"chr": "1",
         "interval": {"start": "q13.32!", "end": "q13.32"},
         "species_id": "taxonomy:9606"
         },
        {"chr": "1",
         "interval": {"start": "q13.32", "end": "q13.32"},
         "species": "taxonomy:9606"
         }
    ]

    for invalid_param in invalid_params:
        with pytest.raises(pydantic.error_wrappers.ValidationError):
            ChromosomeLocation(**invalid_param)


def test_sequence_location(sequence_location, sequence_interval):
    """Test that Sequence Location model works correctly."""
    assert sequence_location.sequence_id == \
           "ga4gh:SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl"
    assert sequence_location.interval.start.value == 44908821
    assert sequence_location.interval.end.value == 44908822
    assert sequence_location.type == "SequenceLocation"

    s = SequenceLocation(_id="sequence:id", sequence_id="refseq:NC_000007.13",
                         interval=sequence_interval, type="SequenceLocation")
    assert s.id == "sequence:id"
    assert s.sequence_id == "refseq:NC_000007.13"
    assert sequence_location.type == "SequenceLocation"

    params = {
        "id": "sequence:id",
        "sequence_id": "refseq:NC_000007.13",
        "interval": sequence_interval
    }
    assert SequenceLocation(**params)

    invalid_params = [
        {
            "_id": "sequence",
            "sequence_id": "NC_000007.13",
            "interval": sequence_interval
        },
        {
            "_id": "sequence:1",
            "sequence_id": "NC_000007.13",
            "interval": sequence_interval,
            "type": "ChromosomeLocation"
        },
        {
            "id": "sequence:id",
            "sequence_id": "test:1",
            "interval": {
                "type": "SimpleInterval",
                "start": 1,
                "end": 2
            }
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


def test_repeated_sequence_expression(derived_sequence_expression, number,
                                      definite_range, indefinite_range):
    """Test that Repeated Sequence Expression model works correctly."""

    def _check_seq_expr(r):
        """Test that seq_expr has intended values."""
        assert r.seq_expr.location.sequence_id == \
               "ga4gh:SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl"
        assert r.seq_expr.location.interval.start.value == 44908821
        assert r.seq_expr.location.interval.end.value == 44908822
        assert r.seq_expr.reverse_complement is False

    r = RepeatedSequenceExpression(
        seq_expr=derived_sequence_expression,
        count=number
    )
    _check_seq_expr(r)
    assert r.count.value == 3
    assert r.count.type == "Number"
    assert r.type == "RepeatedSequenceExpression"

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


def test_haplotype(allele):
    """Test that Haplotype model works correctly."""
    h = Haplotype(members=["ga4gh:VA.-kUJh47Pu24Y3Wdsk1rXEDKsXWNY-68x",
                           "ga4gh:VA.Z_rYRxpUvwqCLsCBO3YLl70o2uf9_Op1"])
    assert len(h.members) == 2

    h = Haplotype(members=[allele], type="Haplotype")
    assert len(h.members) == 1

    invalid_params = [
        {"members": [allele], "type": "Alleles"},
        {"members": allele},
        {"members": [allele, "ga4ghVA"]},
        {"members": [allele], "id_": "ga4gh:VA"}
    ]

    for invalid_param in invalid_params:
        with pytest.raises(pydantic.error_wrappers.ValidationError):
            Haplotype(**invalid_param)


def test_copy_number(number, definite_range, indefinite_range, gene,
                     allele):
    """Test that Copy Number model works correctly."""
    c = CopyNumber(subject=gene, copies=number)
    assert c.subject.gene_id == "ncbigene:348"
    assert c.copies.value == 3
    assert c.type == "CopyNumber"

    c = CopyNumber(subject=gene, copies=definite_range, type="CopyNumber")
    assert c.subject.gene_id == "ncbigene:348"
    assert c.copies.min == 22
    assert c.copies.max == 33

    c = CopyNumber(subject=gene, copies=indefinite_range)
    assert c.subject.gene_id == "ncbigene:348"
    assert c.copies.value == 3
    assert c.copies.comparator == ">="

    invalid_params = [
        {"subjects": number, "copies": number},
        {"ID": "ga4gh:id", "subject": gene, "copies": number},
        {"subject": [allele], "copies": number}
    ]

    for invalid_param in invalid_params:
        with pytest.raises(pydantic.error_wrappers.ValidationError):
            CopyNumber(**invalid_param)


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
