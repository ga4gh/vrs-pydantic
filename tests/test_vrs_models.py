"""Module for testing the VRS model."""
import pydantic
import pytest

from ga4gh.vrsatile.pydantic.vrs_models import ComposedSequenceExpression, Number, \
    Comparator, IndefiniteRange, DefiniteRange, Text, SequenceInterval, \
    CytobandInterval, Gene, Haplotype, Allele, DerivedSequenceExpression, \
    SimpleInterval, SequenceState, Feature, VariationSet, \
    LiteralSequenceExpression, RepeatedSequenceExpression, SystemicVariation, \
    SequenceLocation, ChromosomeLocation, CopyNumberChange, CopyNumberCount


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
        with pytest.raises(pydantic.ValidationError):
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
    assert IndefiniteRange(value=2.0, comparator="<=", type="IndefiniteRange")
    assert IndefiniteRange(value=2.342342, comparator="<=", type="IndefiniteRange")

    invalid_params = [
        {"value": "3", "comparator": Comparator.LT_OR_EQUAL},
        {"value": "3.0", "comparator": Comparator.LT_OR_EQUAL},
        {"value": "3.2323", "comparator": Comparator.LT_OR_EQUAL},
        {"values": 2, "comparator": Comparator.LT_OR_EQUAL},
        {"value": 3, "comparator": "=="},
        {"value": 2, "comparator": Comparator.LT_OR_EQUAL, "type": "s"}
    ]

    for invalid_param in invalid_params:
        with pytest.raises(pydantic.ValidationError):
            IndefiniteRange(**invalid_param)


def test_definite_range(definite_range):
    """Test that Definite Range model works correctly."""
    assert definite_range.min == 22
    assert definite_range.max == 33
    assert definite_range.type == "DefiniteRange"

    assert DefiniteRange(min=0, max=2, type="DefiniteRange")
    assert DefiniteRange(min=0.0, max=2.0)
    assert DefiniteRange(min=0.01, max=2.1241241)

    invalid_params = [
        {"min": 22, "max": 33, "type": "IndefiniteRange"},
        {"min": "22.0", "max": 33},
        {"min": "22.0", "max": 33},
        {"min": "22.230", "max": 33},
        {"min": 22.0, "max": "33"},
        {"min": 22.0, "max": "33.0"},
        {"min": 22.0, "max": "33.230"}
    ]

    for invalid_param in invalid_params:
        with pytest.raises(pydantic.ValidationError):
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
    t_dict = t.model_dump(by_alias=True)
    assert t_dict['_id'] == 'ga4gh:id'
    assert 'id' not in t_dict.keys()

    params = {"definition": definition, "_id": "ga4gh:id"}
    assert Text(**params)

    invalid_params = [
        {"definition": definition, "type": "Definition"},
        {"definition": definition, "_id": "id"}
    ]

    for invalid_param in invalid_params:
        with pytest.raises(pydantic.ValidationError):
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

    assert SequenceInterval(
        start=IndefiniteRange(value=2, comparator="<=", type="IndefiniteRange"),
        end=Number(value=44908822),
        type="SequenceInterval"
    )

    assert SequenceInterval(
        start=DefiniteRange(min=0, max=200),
        end=Number(value=44908822),
        type="SequenceInterval"
    )

    invalid_params = [
        {"start": 44908821, "end": 44908822},
        {"start": {"value": 2}, "end": {"value": 3}, "type": "M"},
        {
            "start": {"value": -1, "type": "Number"},
            "end": {"value": 44908822, "type": "Number"}
        },
        {
            "start": {"value": 44908823, "type": "Number"},
            "end": {"value": 44908822, "type": "Number"}
        },
        {
            "start": {"min": 0, "max": 44908823, "type": "DefiniteRange"},
            "end": {"value": 44908822, "type": "Number"}
        }
    ]

    for invalid_param in invalid_params:
        with pytest.raises(pydantic.ValidationError):
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
        with pytest.raises(pydantic.ValidationError):
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
        with pytest.raises(pydantic.ValidationError):
            LiteralSequenceExpression(**invalid_param)


def test_composed_sequence_expression(derived_sequence_expression, number):
    """Test that Composed Sequence Expression model works correctly"""
    assert ComposedSequenceExpression(
        components=[
            LiteralSequenceExpression(sequence="ACGT",
                                      type="LiteralSequenceExpression"),
            RepeatedSequenceExpression(
                seq_expr=derived_sequence_expression,
                count=number
            )
        ]
    )

    with pytest.raises(pydantic.ValidationError):
        ComposedSequenceExpression(components=[
            LiteralSequenceExpression(sequence="ACGT",
                                      type="LiteralSequenceExpression"),
            LiteralSequenceExpression(sequence="AT", type="LiteralSequenceExpression"),
        ])


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
        with pytest.raises(pydantic.ValidationError):
            Gene(**invalid_param)


def test_chromosome_location(chromosome_location, cytoband_interval):
    """Test that Chromosome Location model works correctly."""
    assert chromosome_location.chr == "19"
    assert chromosome_location.interval.start == chromosome_location.interval.end == "q13.32"  # noqa: E501
    assert chromosome_location.type == "ChromosomeLocation"

    chrs = [str(i) for i in range(1, 23)] + ["X", "Y"]
    for chr in chrs:
        cl = ChromosomeLocation(
            chr=chr,
            interval=cytoband_interval,
        )
        assert cl
        assert cl.type == "ChromosomeLocation"
        assert cl.chr == chr
        assert cl.species_id == "taxonomy:9606"

    invalid_params = [
        {
            "chr": "1",
            "interval": {"start": "q13.32!", "end": "q13.32"},
            "species_id": "taxonomy:9606"
        },
        {
            "chr": "1",
            "interval": {"start": "q13.32", "end": "q13.32"},
            "species": "taxonomy:9606"
        },
        {
            "chr": "Z",
            "interval": {"start": "q13.32", "end": "q13.32"}
        }
    ]

    for invalid_param in invalid_params:
        with pytest.raises(pydantic.ValidationError):
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
        "_id": "sequence:id",
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
        }
    ]

    for invalid_param in invalid_params:
        with pytest.raises(pydantic.ValidationError):
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
        with pytest.raises(pydantic.ValidationError):
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
        with pytest.raises(pydantic.ValidationError):
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
        with pytest.raises(pydantic.ValidationError):
            Allele(**invalid_param)


def test_haplotype(allele):
    """Test that Haplotype model works correctly."""
    h = Haplotype(members=["ga4gh:VA.-kUJh47Pu24Y3Wdsk1rXEDKsXWNY-68x",
                           "ga4gh:VA.Z_rYRxpUvwqCLsCBO3YLl70o2uf9_Op1"])
    assert len(h.members) == 2

    invalid_params = [
        {"members": [allele], "type": "Haplotype"},  # members must have len >= 2
        {"members": allele},  # members is a list
        {"members": [allele, "ga4ghVA"]},  # Not a CURIE
        {"members": [allele], "id_": "ga4gh:VA"},  # invalid property name id_
        {"members": [allele, allele]},  # not unique
        {"members": ["ga4gh:VA.-kUJh47Pu24Y3Wdsk1rXEDKsXWNY-68x",
                     "ga4gh:VA.-kUJh47Pu24Y3Wdsk1rXEDKsXWNY-68x"]},  # not unique
    ]

    for invalid_param in invalid_params:
        with pytest.raises(pydantic.ValidationError):
            Haplotype(**invalid_param)


def test_copy_number_count(number, definite_range, indefinite_range, gene, allele,
                           sequence_location):
    """Test that Copy Number Count model works correctly."""
    c = CopyNumberCount(subject=gene, copies=number)
    assert c.subject.gene_id == "ncbigene:348"
    assert c.copies.value == 3
    assert c.type == "CopyNumberCount"

    c = CopyNumberCount(
        subject=gene, copies=definite_range, type="CopyNumberCount")
    assert c.subject.gene_id == "ncbigene:348"
    assert c.copies.min == 22
    assert c.copies.max == 33

    c = CopyNumberCount(subject=gene, copies=indefinite_range)
    assert c.subject.gene_id == "ncbigene:348"
    assert c.copies.value == 3
    assert c.copies.comparator == ">="

    c = CopyNumberCount(subject=sequence_location, copies=number)
    assert c.subject.type == "SequenceLocation"

    invalid_params = [
        {"subjects": number, "copies": number},
        {"ID": "ga4gh:id", "subject": gene, "copies": number},
        {"subject": [allele], "copies": number}
    ]

    for invalid_param in invalid_params:
        with pytest.raises(pydantic.ValidationError):
            CopyNumberCount(**invalid_param)


def test_copy_number_change(number, sequence_location, gene, allele):
    """Test that Copy Number Change model works correctly."""
    c = CopyNumberChange(subject=gene, copy_change="efo:0030069")
    assert c.subject.gene_id == "ncbigene:348"
    assert c.copy_change == "efo:0030069"
    assert c.type == "CopyNumberChange"

    c = CopyNumberChange(subject=sequence_location, copy_change="efo:0030071")
    assert c.subject.type == "SequenceLocation"
    assert c.copy_change == "efo:0030071"

    invalid_params = [
        {"subject": number, "copies": number},
        {"ID": "ga4gh:id", "subject": gene, "copy_change": "efo:0030069"},
        {"subject": allele},
        {"subject": "fake:curie", "copy_change": "efo:0030071", "extra": 0},
        {"subject": "fake:curie", "copy_change": "low-level"},
        {"subject": allele, "copy_change": "partial loss"},
    ]

    for invalid_param in invalid_params:
        with pytest.raises(pydantic.ValidationError):
            CopyNumberChange(**invalid_param)


def test_variation_set(allele, sequence_location):
    """Test that Variation Set model works correctly."""
    v = VariationSet(members=[])
    assert len(v.members) == 0
    assert v.type == "VariationSet"

    assert VariationSet(members=[allele, "fake:curie", Text(definition="def")])

    v = VariationSet(members=[allele], type="VariationSet")
    assert len(v.members) == 1
    assert v.members[0].type == "Allele"
    assert v.members[0].location.type == "SequenceLocation"
    assert v.members[0].location.sequence_id == \
           "ga4gh:SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl"

    invalid_params = [
        {"members": [1]},
        {"members": [allele, sequence_location]},
        {"members": [allele], "type": "VS"},
        {"members": [allele, allele]},  # not unique
        {"members": ["fake:curie", "fake:curie"]}  # not unique
    ]

    for invalid_param in invalid_params:
        with pytest.raises(pydantic.ValidationError):
            VariationSet(**invalid_param)


def test_feature(gene):
    """Test Feature class."""
    assert Feature(root=gene)


def test_systemic_variation(gene, number):
    """Test SystemicVariation class."""
    c = CopyNumberCount(subject=gene, copies=number)
    assert SystemicVariation(root=c)


def test_deprecated_objects(caplog, deprecated_allele):
    """Test that deprecated objects work and log appropriately."""
    seqstate_deprecated_msg = "SequenceState is deprecated. Use default='LiteralSequenceExpression' instead."  # noqa: E501
    simpleint_deprecated_msg = "SimpleInterval is deprecated. Use default='SequenceInterval' instead."  # noqa: E501
    seqstate = SequenceState(**deprecated_allele["state"])
    assert seqstate.type == "SequenceState"
    assert seqstate.sequence == "T"
    assert seqstate_deprecated_msg in caplog.text

    invalid_params = [
        {"sequence": "t"},
        {"sequence": "T", "type": "Sequence"},
        {"sequence": "hello,world"}
    ]
    for invalid_param in invalid_params:
        with pytest.raises(pydantic.ValidationError):
            SequenceState(**invalid_param)

    simpleint = SimpleInterval(**deprecated_allele["location"]["interval"])
    assert simpleint.type == "SimpleInterval"
    assert simpleint.start == 140753335
    assert simpleint.end == 140753336
    assert simpleint_deprecated_msg in caplog.text

    invalid_params = [
        {"start": 2.0, "end": 2},
        {"start": 2, "end": 2, "type": "CytobandInterval"},
        {"start": 2, "end": '2'}
    ]
    for invalid_param in invalid_params:
        with pytest.raises(pydantic.ValidationError):
            SimpleInterval(**invalid_param)

    allele = Allele(**deprecated_allele)
    assert allele.state.type == "SequenceState"
    assert allele.state.sequence == "T"
    assert seqstate_deprecated_msg in caplog.text
    assert allele.location.interval.type == "SimpleInterval"
    assert allele.location.interval.start == 140753335
    assert allele.location.interval.end == 140753336
    assert simpleint_deprecated_msg in caplog.text

    # should default to non-deprecated option when possible
    allele = Allele(state={"sequence": "T"},
                    location=deprecated_allele["location"])
    assert allele.state.type == "LiteralSequenceExpression"
