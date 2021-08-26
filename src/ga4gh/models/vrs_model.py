from __future__ import annotations
from enum import Enum
from typing import List, Optional, Union
from pydantic import BaseModel, Extra, Field, constr, StrictInt, StrictStr, \
    StrictBool


class Number(BaseModel):
    """A simple number value as a VRS class."""
    class Config:
        extra = Extra.forbid

    type = "Number"
    value: StrictInt


class Comparator(str, Enum):
    """A range comparator."""
    LT_OR_EQUAL = "<="
    GT_OR_EQUAL = ">="


class IndefiniteRange(BaseModel):
    """An indefinite range represented as a number and associated comparator.
    The bound operator is interpreted as follows: `>=` are all values greater
    than and including the value, `<=` are all numbers less than and including
    the value."""
    class Config:
        extra = Extra.forbid

    type = "IndefiniteRange"
    value: StrictInt
    comparator: Comparator


class DefiniteRange(BaseModel):
    """A bounded, inclusive range of numbers."""
    class Config:
        extra = Extra.forbid

    type = "DefiniteRange"
    min: StrictInt
    max: StrictInt


class Sequence(BaseModel):
    """A character string of residues that represents a biological sequence
    using the conventional sequence order (5’-to-3’ for nucleic acid sequences,
     and amino-to-carboxyl for amino acid sequences). IUPAC ambiguity codes
     are permitted in Sequences."""

    class Config:
        extra = Extra.forbid

    __root__: constr(regex=r"^[A-Z*\-]*$") = Field(..., example="ACTG")


class CURIE(BaseModel):
    """A string that refers to an object uniquely.  The lifetime and scope of
    an id is defined by the sender. VRS does not impose any contraints on
    strings used as ids in messages. However, to maximize sharability of data,
    VRS RECOMMENDS that implementations use [W3C Compact URI (CURIE)]
    (https://www.w3.org/TR/curie/) syntax. String CURIEs are represented as
    `prefix`:`reference` (W3C terminology), but often referred to as `
    namespace`:`accession` or `namespace`:`local id` colloquially. VRS also
    RECOMMENDS that `prefix` be defined in identifiers.org. The `reference`
    component is an unconstrained string. A CURIE is a URI.  URIs may *locate*
    objects (i.e., specify where to retrieve them) or *name* objects
    conceptually.  VRS uses CURIEs primarily as a naming mechanism.
    Implementations MAY provide CURIE resolution mechanisms for prefixes to
    make these objects locatable. Using internal ids in public messages is
    strongly discouraged."""
    class Config:
        extra = Extra.forbid

    __root__: constr(regex=r"^\w[^:]*:.+$") = \
        Field(..., example="ensembl:ENSG00000139618")


class HumanCytoband(BaseModel):
    """A interval on a stained metaphase chromosome specified by cytobands.
    CytobandIntervals include the regions described by the start and end
    cytobands."""
    class Config:
        extra = Extra.forbid

    __root__: constr(regex=r"^cen|[pq](ter|([1-9][0-9]*(\.[1-9][0-9]*)?))$") \
        = Field(..., example="q22.3")


class SequenceState(BaseModel):
    class Config:
        extra = Extra.forbid

    type = "SequenceState"
    sequence: Sequence


class SimpleInterval(BaseModel):
    class Config:
        extra = Extra.forbid

    type = "SimpleInterval"
    start: StrictInt
    end: StrictInt


class Text(BaseModel):
    """An textual representation of variation intended to capture variation
    descriptions that cannot be parsed, but still treated as variation."""
    class Config:
        extra = Extra.forbid

    id: Optional[CURIE] = Field(alias='_id')
    type = "Text"
    definition: StrictStr


class SequenceInterval(BaseModel):
    class Config:
        extra = Extra.forbid

    type = "SequenceInterval"
    start: Union[Number, IndefiniteRange, DefiniteRange]
    end: Union[Number, IndefiniteRange, DefiniteRange]


class CytobandInterval(BaseModel):
    class Config:
        extra = Extra.forbid

    type = "CytobandInterval"
    start: HumanCytoband
    end: HumanCytoband


class LiteralSequenceExpression(BaseModel):
    class Config:
        extra = Extra.forbid

    type = "LiteralSequenceExpression"
    sequence: Sequence


class Gene(BaseModel):
    class Config:
        extra = Extra.forbid

    type = "Gene"
    gene_id: CURIE


class ChromosomeLocation(BaseModel):
    class Config:
        extra = Extra.forbid

    type = "ChromosomeLocation"
    id: Optional[CURIE] = Field(alias='_id')
    species_id: CURIE
    chr: StrictStr
    interval: CytobandInterval


class SequenceLocation(BaseModel):
    class Config:
        extra = Extra.forbid

    id: Optional[CURIE] = Field(alias='_id')
    type = "SequenceLocation"
    sequence_id: CURIE
    interval: Union[SequenceInterval, SimpleInterval]


class DerivedSequenceExpression(BaseModel):
    class Config:
        extra = Extra.forbid

    type= "DerivedSequenceExpression"
    location: SequenceLocation
    reverse_complement: StrictBool


class RepeatedSequenceExpression(BaseModel):
    class Config:
        extra = Extra.forbid

    type = "RepeatedSequenceExpression"
    seq_expr: Union[LiteralSequenceExpression, DerivedSequenceExpression]
    count: Union[Number, IndefiniteRange, DefiniteRange]


class Feature(BaseModel):
    """A named entity that can be mapped to a Location. Genes, protein domains,
     exons, and chromosomes are some examples of common biological entities
     that may be Features."""
    __root__: Gene


class Location(BaseModel):
    """A contiguous segment of a biological sequence."""
    __root__: Union[ChromosomeLocation, SequenceLocation]


class SequenceExpression(BaseModel):
    """One of a set of sequence representation syntaxes."""
    __root__: Union[
        LiteralSequenceExpression,
        DerivedSequenceExpression,
        RepeatedSequenceExpression
    ]


class Allele(BaseModel):
    class Config:
        extra = Extra.forbid

    id: Optional[CURIE] = Field(alias='_id')
    type = "Allele"
    location: Union[CURIE, Location]
    state: Union[SequenceState, SequenceExpression]


class Haplotype(BaseModel):
    class Config:
        extra = Extra.forbid

    id: Optional[CURIE] = Field(alias='_id')
    type = "Haplotype"
    members: List[Union[Allele, CURIE]] = Field(..., min_items=1)


class MolecularVariation(BaseModel):
    """A variation on a contiguous molecule."""
    __root__: Union[Allele, Haplotype]


class CopyNumber(BaseModel):
    class Config:
        extra = Extra.forbid

    id: Optional[CURIE] = Field(alias='_id')
    type = "CopyNumber"
    subject: Union[MolecularVariation, Feature, SequenceExpression, CURIE]
    copies: Union[Number, IndefiniteRange, DefiniteRange]


class SystemicVariation(BaseModel):
    """A Variation of multiple molecules in the context of a system,
    e.g. a genome, sample, or homologous chromosomes."""
    __root__: CopyNumber


class Variation(BaseModel):
    """A representation of the state of one or more biomolecules."""
    __root__: Union[MolecularVariation, SystemicVariation, UtilityVariation]


class UtilityVariation(BaseModel):
    """A collection of :ref:`Variation` subclasses that cannot be constrained
     to a specific class of biological variation, but are necessary for some
     applications of VRS."""
    __root__: Union[Text, VariationSet]


class VariationSet(BaseModel):
    class Config:
        extra = Extra.forbid

    id: Optional[CURIE] = Field(alias='_id')
    type = "VariationSet"
    members: List[Union[CURIE, Variation]]


Variation.update_forward_refs()
UtilityVariation.update_forward_refs()
