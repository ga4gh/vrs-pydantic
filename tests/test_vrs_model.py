"""Module for testing the VRS model."""
import pydantic
import pytest
from ga4gh.models.vrs_model import Number, Comparator, \
    IndefiniteRange, DefiniteRange, SequenceState, SimpleInterval, Text, \
    SequenceInterval, CytobandInterval, DerivedSequenceExpression, \
    LiteralSequenceExpression, RepeatedSequenceExpression, Gene, \
    ChromosomeLocation, SequenceLocation, VariationSet, Haplotype, \
    CopyNumber, Allele


def test_number(number):
    assert number.value == 3
    assert number.type == "Number"

    invalid_params = [
        {"value": '2'},
        {"value": 2.0},
        {"value": 2, "n": 1}
    ]

    for invalid_param in invalid_params:
        with pytest.raises(pydantic.error_wrappers.ValidationError):
            Number(**invalid_param)


def test_comparator():
    assert Comparator.LT_OR_EQUAL == "<="
    assert Comparator.GT_OR_EQUAL == ">="


def test_indefinite_range(indefinite_range):
    assert indefinite_range.value == 3
    assert indefinite_range.comparator == Comparator.GT_OR_EQUAL

    invalid_params = [
        {"value": "3", "comparator": Comparator.LT_OR_EQUAL},
        {"values": 2, "comparator": Comparator.LT_OR_EQUAL},
        {"value": 3, "comparator": "=="}
    ]

    for invalid_param in invalid_params:
        with pytest.raises(pydantic.error_wrappers.ValidationError):
            IndefiniteRange(**invalid_param)


def test_definite_range(definite_range):
    assert definite_range.min == 22
    assert definite_range.max == 33

    invalid_params = [
        {"min": 22.0, "max": 33}
    ]

    for invalid_param in invalid_params:
        with pytest.raises(pydantic.error_wrappers.ValidationError):
            DefiniteRange(**invalid_param)


def test_sequence_state():
    s = SequenceState(sequence="T")
    assert s.sequence.__root__ == 'T'
    assert s.type == "SequenceState"

    invalid_params = [
        {"sequence": "t"},
        {"sequence": "hello,world"}
    ]

    for invalid_param in invalid_params:
        with pytest.raises(pydantic.error_wrappers.ValidationError):
            SequenceState(**invalid_param)


def test_simple_interval():
    s = SimpleInterval(start=2, end=2)
    assert s.start == 2
    assert s.end == 2
    assert s.type == "SimpleInterval"

    invalid_params = [
        {"start": 2.0, "end": 2},
        {"start": 2, "end": '2'}
    ]

    for invalid_param in invalid_params:
        with pytest.raises(pydantic.error_wrappers.ValidationError):
            SimpleInterval(**invalid_param)


def test_text():
    definition = "APOE_LOSS"
    t = Text(definition=definition)
    assert t.definition == definition
    assert t.type == "Text"

    t = Text(_id="ga4gh:id", definition=definition)
    assert t.definition == definition
    assert t.type == "Text"
    assert t.id.__root__ == "ga4gh:id"
    t_dict = t.dict(by_alias=True)
    assert t_dict['_id'] == 'ga4gh:id'
    assert 'id' not in t_dict.keys()

    invalid_params = [
        {"definition": definition, "_id": "id"},
        {"definition": definition, "id": "ga4gh:id"}
    ]

    for invalid_param in invalid_params:
        with pytest.raises(pydantic.error_wrappers.ValidationError):
            Text(**invalid_param)


def test_sequence_interval(sequence_interval):
    assert sequence_interval.start.value == 44908821
    assert sequence_interval.end.value == 44908822
    assert sequence_interval.type == "SequenceInterval"

    invalid_params = [
        {"start": 44908821, "end": 44908822}
    ]

    for invalid_param in invalid_params:
        with pytest.raises(pydantic.error_wrappers.ValidationError):
            SequenceInterval(**invalid_param)


def test_cytoband_interval(cytoband_interval):
    human_cytoband = "q13.32"
    assert cytoband_interval.start.__root__ == human_cytoband == cytoband_interval.end.__root__

    invalid_params = [
        {"start": "x9", "end": human_cytoband}
    ]

    for invalid_param in invalid_params:
        with pytest.raises(pydantic.error_wrappers.ValidationError):
            CytobandInterval(**invalid_param)


def test_literal_sequence_expression():
    l = LiteralSequenceExpression(sequence="ACGT")
    assert l.sequence.__root__ == "ACGT"

    invalid_params = [
        {"sequence": "actg"},
        {"sequence": "ACTx"}
    ]

    for invalid_param in invalid_params:
        with pytest.raises(pydantic.error_wrappers.ValidationError):
            LiteralSequenceExpression(**invalid_param)


def test_gene():
    g = Gene(gene_id="hgnc:5")
    assert g.gene_id.__root__ == "hgnc:5"
    assert g.type == "Gene"

    invalid_params = [
        {"gene_id": "hgnc5"}
    ]

    for invalid_param in invalid_params:
        with pytest.raises(pydantic.error_wrappers.ValidationError):
            Gene(**invalid_param)


def test_chromosome_location(cytoband_interval):
    human_cytoband = "q13.32"
    human_taxonomy = "taxonomy:9606"
    c = ChromosomeLocation(
        chr="19",
        interval=cytoband_interval,
        species_id=human_taxonomy)
    assert c.chr == "19"
    assert c.interval.start.__root__ == c.interval.end.__root__ == "q13.32"

    invalid_params = [
        {"chr": "K",
         "interval": {"start": human_cytoband, "end": human_cytoband},
         "species_id": human_taxonomy
         }
    ]

    for invalid_param in invalid_params:
        with pytest.raises(pydantic.error_wrappers.ValidationError):
            CytobandInterval(**invalid_param)


def test_sequence_location(sequence_location, sequence_interval):
    assert sequence_location.sequence_id.__root__ == \
           "ga4gh:SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl"
    assert sequence_location.interval.start.value == 44908821
    assert sequence_location.interval.end.value == 44908822
    assert sequence_location.type == "SequenceLocation"

    s = SequenceLocation(_id="sequence:id", sequence_id="refseq:NC_000007.13",
                         interval=sequence_interval)
    assert s.id.__root__ == "sequence:id"
    assert s.sequence_id.__root__ == "refseq:NC_000007.13"

    invalid_params = [
        {
            "_id": "sequence",
            "sequence_id": "NC_000007.13",
            "interval": sequence_interval
         },
        {
            "id": "sequence:id",
            "sequence_id": "refseq:NC_000007.13",
            "interval": sequence_interval
        }
    ]

    for invalid_param in invalid_params:
        with pytest.raises(pydantic.error_wrappers.ValidationError):
            SequenceLocation(**invalid_param)


def test_derived_sequence_expression(sequence_location,
                                     derived_sequence_expression):
    assert derived_sequence_expression.reverse_complement is False

    d = DerivedSequenceExpression(
        location=sequence_location,
        reverse_complement=True,
    )
    assert d.reverse_complement is True

    invalid_params = [
        {"location": sequence_location, "reverse_complement": 0}
    ]

    for invalid_param in invalid_params:
        with pytest.raises(pydantic.error_wrappers.ValidationError):
            DerivedSequenceExpression(**invalid_param)


def test_repeated_sequence_expression(derived_sequence_expression, number,
                                      definite_range, indefinite_range):

    def _check_seq_expr(r):
        assert r.seq_expr.location.sequence_id.__root__ == \
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

    r = RepeatedSequenceExpression(
        seq_expr=derived_sequence_expression,
        count=definite_range
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
        {"seq_expr": derived_sequence_expression, "count": 2}
    ]

    for invalid_param in invalid_params:
        with pytest.raises(pydantic.error_wrappers.ValidationError):
            RepeatedSequenceExpression(**invalid_param)


def test_allele(allele, sequence_location, derived_sequence_expression):
    assert allele.type == "Allele"
    assert allele.location.__root__.type == "SequenceLocation"
    assert allele.location.__root__.sequence_id.__root__ == \
           "ga4gh:SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl"

    a = Allele(location="ga4gh:location", state=derived_sequence_expression)
    assert a.location.__root__ == "ga4gh:location"

    invalid_params = [
        {"location": sequence_location, "state": sequence_location}
    ]

    for invalid_param in invalid_params:
        with pytest.raises(pydantic.error_wrappers.ValidationError):
            Allele(**invalid_param)


def test_haplotype(allele):
    h = Haplotype(members=["ga4gh:VA.-kUJh47Pu24Y3Wdsk1rXEDKsXWNY-68x",
                           "ga4gh:VA.Z_rYRxpUvwqCLsCBO3YLl70o2uf9_Op1"])
    assert len(h.members) == 2

    h = Haplotype(members=[allele])
    assert len(h.members) == 1

    invalid_params = [
        {"members": allele},
        {"members": [allele, "ga4ghVA"]},
        {"members": [allele], "id_": "ga4gh:VA"}
    ]

    for invalid_param in invalid_params:
        with pytest.raises(pydantic.error_wrappers.ValidationError):
            Haplotype(**invalid_param)


def test_copy_number_variation(number, definite_range, indefinite_range, gene, allele):
    c = CopyNumber(subject=gene, copies=number)
    assert c.subject.__root__.gene_id.__root__ == "ncbigene:348"
    assert c.copies.value == 3

    c = CopyNumber(subject=gene, copies=definite_range)
    assert c.subject.__root__.gene_id.__root__ == "ncbigene:348"
    assert c.copies.min == 22
    assert c.copies.max == 33

    c = CopyNumber(subject=gene, copies=indefinite_range)
    assert c.subject.__root__.gene_id.__root__ == "ncbigene:348"
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
    v = VariationSet(members=[])
    v = VariationSet(members=[allele, allele])

    invalid_params = [
        {"members": [1]},
        {"members": [allele, sequence_location]}
    ]

    for invalid_param in invalid_params:
        with pytest.raises(pydantic.error_wrappers.ValidationError):
            VariationSet(**invalid_param)
