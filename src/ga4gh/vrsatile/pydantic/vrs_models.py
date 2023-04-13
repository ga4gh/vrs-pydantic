"""Define Pydantic Class models for VRS models."""
from __future__ import annotations

import re
from enum import Enum
from typing import List, Optional, Union, Literal

from pydantic import BaseModel, Field, constr, StrictInt, StrictStr, StrictBool, \
    validator, StrictFloat, root_validator
from pydantic.error_wrappers import ValidationError

from ga4gh.vrsatile.pydantic import return_value, BaseModelForbidExtra, \
    BaseModelDeprecated


class VRSTypes(str, Enum):
    """Define types used in VRS."""

    ALLELE = "Allele"
    HAPLOTYPE = "Haplotype"
    TEXT = "Text"
    VARIATION_SET = "VariationSet"
    COPY_NUMBER_COUNT = "CopyNumberCount"
    COPY_NUMBER_CHANGE = "CopyNumberChange"
    GENOTYPE = "Genotype"
    CHROMOSOME_LOCATION = "ChromosomeLocation"
    SEQUENCE_LOCATION = "SequenceLocation"
    SEQUENCE_INTERVAL = "SequenceInterval"
    CYTOBAND_INTERVAL = "CytobandInterval"
    LITERAL_SEQUENCE_EXPRESSION = "LiteralSequenceExpression"
    DERIVED_SEQUENCE_EXPRESSION = "DerivedSequenceExpression"
    REPEATED_SEQUENCE_EXPRESSION = "RepeatedSequenceExpression"
    COMPOSED_SEQUENCE_EXPRESSION = "ComposedSequenceExpression"
    GENOTYPE_MEMBER = "GenotypeMember"
    GENE = "Gene"
    NUMBER = "Number"
    DEFINITE_RANGE = "DefiniteRange"
    INDEFINITE_RANGE = "IndefiniteRange"
    SEQUENCE_STATE = "SequenceState"  # DEPRECATED
    SIMPLE_INTERVAL = "SimpleInterval"  # DEPRECATED


class Comparator(str, Enum):
    """A range comparator."""

    LT_OR_EQUAL = "<="
    GT_OR_EQUAL = ">="


class CopyChange(str, Enum):
    """The copy change (https://www.ebi.ac.uk/efo/)"""

    COMPLETE_GENOMIC_LOSS = "efo:0030069"
    HIGH_LEVEL_LOSS = "efo:0020073"
    LOW_LEVEL_LOSS = "efo:0030068"
    LOSS = "efo:0030067"
    REGIONAL_BASE_PLOIDY = "efo:0030064"
    GAIN = "efo:0030070"
    LOW_LEVEL_GAIN = "efo:0030071"
    HIGH_LEVEL_GAIN = "efo:0030072"


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# BASIC TYPES (STRUCTURES)
# These types do NOT have a VRS `type` attribute
# These types are used solely within other definitions.


class CURIE(BaseModelForbidExtra):
    """A `W3C Compact URI <https://www.w3.org/TR/curie/>` formatted string.  A CURIE
    string has the structure ``prefix``:``reference``, as defined by the W3C syntax.
    """

    __root__: constr(regex=r"^\w[^:]*:.+$") = Field(  # noqa: F722
        ...,
        description=("A `W3C Compact URI <https://www.w3.org/TR/curie/>`_ formatted "
                     "string. A CURIE string has the structure "
                     "``prefix``:``reference``, as defined by the W3C syntax."),
        example="ensembl:ENSG00000139618"
    )


class HumanCytoband(BaseModelForbidExtra):
    """A character string representing cytobands derived from the *International System
    for Human Cytogenomic Nomenclature* (ISCN)
    [guidelines](http://doi.org/10.1159/isbn.978-3-318-06861-0).
    """

    __root__: constr(regex=r"^cen|[pq](ter|([1-9][0-9]*(\.[1-9][0-9]*)?))$") = Field(  # noqa: F722, E501
        ...,
        description=("A character string representing cytobands derived from the "
                     "*International System for Human Cytogenomic Nomenclature* (ISCN) "
                     "`guidelines <http://doi.org/10.1159/isbn.978-3-318-06861-0>`_."),
        example="q22.3"
    )


class Residue(BaseModelForbidExtra):
    """A character representing a specific residue (i.e., molecular species) or
    groupings of these ("ambiguity codes"), using `one-letter IUPAC abbreviations
    <https://en.wikipedia.org/wiki/International_Union_of_Pure_and_Applied_Chemistry#Amino_acid_and_nucleotide_base_codes>`
    for nucleic acids and amino acids.
    """

    __root__: constr(regex=r"[A-Z*\-]") = Field(  # noqa: F722
        ...,
        description=("A character representing a specific residue (i.e., molecular "
                     "species) or groupings of these ('ambiguity codes'), using "
                     "`one-letter IUPAC abbreviations "
                     "<https://en.wikipedia.org/wiki/International_Union_of_Pure_and_Applied_Chemistry#Amino_acid_and_nucleotide_base_codes>`_ "  # noqa: E501
                     "for nucleic acids and amino acids.")
    )


class Sequence(BaseModelForbidExtra):
    """A character string of Residues that represents a biological sequence using the
    conventional sequence order (5'-to-3' for nucleic acid sequences, and
    amino-to-carboxyl for amino acid  sequences). IUPAC ambiguity codes are permitted in
    Sequences.
    """

    __root__: constr(regex=r"^[A-Z*\-]*$") = Field(  # noqa: F722
        ...,
        description=("A character string of Residues that represents a biological "
                     "sequence using the conventional sequence order (5'-to-3' for "
                     "nucleic acid sequences, and amino-to-carboxyl for amino acid "
                     "sequences). IUPAC ambiguity codes are permitted in Sequences."),
        example="ACTG"
    )


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Numerics, Comparators, and Ranges


class Number(BaseModelForbidExtra):
    """A simple integer value as a VRS class."""

    type: Literal[VRSTypes.NUMBER] = VRSTypes.NUMBER
    value: StrictInt = Field(..., description="The value represented by Number")


class DefiniteRange(BaseModelForbidExtra):
    """A bounded, inclusive range of numbers."""

    type: Literal[VRSTypes.DEFINITE_RANGE] = VRSTypes.DEFINITE_RANGE
    min: Union[StrictFloat, StrictInt] = Field(
        ...,
        description="The minimum value; inclusive"
    )
    max: Union[StrictFloat, StrictInt] = Field(
        ...,
        description="The maximum value; inclusive"
    )


class IndefiniteRange(BaseModelForbidExtra):
    """A half-bounded range of numbers represented as a number bound and associated
    comparator. The bound operator is interpreted as follows: '>=' are all numbers
    greater than and including `value`, '<=' are all numbers less than and including
    `value`.
    """

    type: Literal[VRSTypes.INDEFINITE_RANGE] = VRSTypes.INDEFINITE_RANGE
    value: Union[StrictFloat, StrictInt] = Field(
        ...,
        description="The bounded value; inclusive"
    )
    comparator: Comparator = Field(
        ...,
        description=("MUST be one of '<=' or '>=', indicating which direction the "
                     "range is indefinite"),
    )


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Intervals


class SequenceInterval(BaseModelForbidExtra):
    """A SequenceInterval represents a span on a Sequence. Positions are always
    represented by contiguous spans using interbase coordinates or coordinate ranges.
    """

    type: Literal[VRSTypes.SEQUENCE_INTERVAL] = VRSTypes.SEQUENCE_INTERVAL
    start: Union[DefiniteRange, IndefiniteRange, Number] = Field(
        ...,
        description=("The start coordinate or range of the interval. The minimum value "
                     "of this coordinate or range is 0. MUST represent a coordinate or "
                     "range less than the value of `end`.")
    )
    end: Union[DefiniteRange, IndefiniteRange, Number] = Field(
        ...,
        description=("The end coordinate or range of the interval. The minimum value "
                     "of this coordinate or range is 0. MUST represent a coordinate or "
                     "range greater than the value of `start`.")
    )

    @root_validator(pre=False, skip_on_failure=True)
    def check_start_end_value(cls, v):
        """Check that start is less than or equal to end and that both have minimum
        value of 0
        """
        start = v.get("start")
        if start.type in {VRSTypes.NUMBER, VRSTypes.INDEFINITE_RANGE}:
            start_values = [start.value]
        else:
            start_values = [start.min, start.max]

        for start_value in start_values:
            assert start_value >= 0, "`start` minimum value is 0"

        end = v.get("end")
        if end.type in {VRSTypes.NUMBER, VRSTypes.INDEFINITE_RANGE}:
            end_values = [end.value]
        else:
            end_values = [end.min, end.max]

        for end_value in end_values:
            assert end_value >= 0, "`end` minimum value is 0"

        for sv in start_values:
            for ev in end_values:
                assert sv <= ev, "`start` must be less than or equal to `end"

        return v


class CytobandInterval(BaseModelForbidExtra):
    """A contiguous span on a chromosome defined by cytoband features. The span includes
    the constituent regions described by the start and end cytobands, as well as any
    intervening regions.
    """

    type: Literal[VRSTypes.CYTOBAND_INTERVAL] = VRSTypes.CYTOBAND_INTERVAL
    start: HumanCytoband = Field(
        ...,
        description=("The start cytoband region. MUST specify a region nearer the "
                     "terminal end (telomere) of the chromosome p-arm than `end`."),
        example="q22.2"
    )
    end: HumanCytoband = Field(
        ...,
        description=("The end cytoband region. MUST specify a region nearer the "
                     "terminal end (telomere) of the chromosome q-arm than `start`."),
        example="q22.3"
    )

    _get_start_val = validator("start", allow_reuse=True)(return_value)
    _get_end_val = validator("end", allow_reuse=True)(return_value)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# DEPRECATED Intervals

class SimpleInterval(BaseModelForbidExtra, BaseModelDeprecated):
    """DEPRECATED: A SimpleInterval represents a span of sequence. Positions are always
    represented by contiguous spans using interbase coordinates.
    This class is deprecated. Use SequenceInterval instead.
    """

    type: Literal[VRSTypes.SIMPLE_INTERVAL] = VRSTypes.SIMPLE_INTERVAL
    start: StrictInt = Field(..., description="The start coordinate", example=11)
    end: StrictInt = Field(..., description="The end coordinate", example=22)

    _replace_with = "SequenceInterval"


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Locations


class ChromosomeLocation(BaseModelForbidExtra):
    """A Location on a chromosome defined by a species and chromosome name."""

    id: Optional[CURIE] = Field(
        alias="_id",
        description="Location Id. MUST be unique within document."
    )
    type: Literal[VRSTypes.CHROMOSOME_LOCATION] = VRSTypes.CHROMOSOME_LOCATION
    species_id: CURIE = Field(
        default="taxonomy:9606",
        description=("CURIE representing a species from the `NCBI species taxonomy "
                     "<https://registry.identifiers.org/registry/taxonomy>`_. "
                     "Default: 'taxonomy:9606' (human)")
    )
    chr: StrictStr = Field(
        ...,
        description=("The symbolic chromosome name. For humans, For humans, chromosome "
                     "names MUST be one of 1..22, X, Y (case-sensitive)")
    )
    interval: CytobandInterval = Field(
        ..., description='The chromosome region defined by a CytobandInterval'
    )

    _get_id_val = validator("id", allow_reuse=True)(return_value)
    _get_species_id_val = validator("species_id", allow_reuse=True)(return_value)

    @validator("chr")
    def check_chr_value(cls, v):
        """Check chr value"""
        msg = "`chr` must be 1..22, X, or Y (case-sensitive)"
        assert re.match(r"^(X|Y|([1-9]|1[0-9]|2[0-2]))$", v), msg
        return v


class SequenceLocation(BaseModelForbidExtra):
    """A Location defined by an interval on a referenced Sequence."""

    id: Optional[CURIE] = Field(
        alias="_id",
        description="Variation Id. MUST be unique within document."
    )
    type: Literal[VRSTypes.SEQUENCE_LOCATION] = VRSTypes.SEQUENCE_LOCATION
    sequence_id: CURIE = Field(
        ...,
        description="A VRS Computed Identifier for the reference Sequence."
    )
    interval: Union[SequenceInterval, SimpleInterval] = Field(
        ...,
        description="Reference sequence region defined by a SequenceInterval."
    )

    _get_id_val = validator("id", allow_reuse=True)(return_value)
    _get_sequence_id_val = validator('sequence_id', allow_reuse=True)(return_value)


class Location(BaseModel):
    """A contiguous segment of a biological sequence."""

    __root__: Union[ChromosomeLocation, SequenceLocation] = Field(
        ..., description="A contiguous segment of a biological sequence."
    )


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# SequenceExpression


class LiteralSequenceExpression(BaseModelForbidExtra):
    """An explicit expression of a Sequence."""

    type: Literal[VRSTypes.LITERAL_SEQUENCE_EXPRESSION] = \
        VRSTypes.LITERAL_SEQUENCE_EXPRESSION
    sequence: Sequence = Field(..., description="the literal Sequence expressed")

    _get_sequence_val = validator("sequence", allow_reuse=True)(return_value)


class DerivedSequenceExpression(BaseModelForbidExtra):
    """An approximate expression of a sequence that is derived from a referenced
    sequence location. Use of this class indicates that the derived sequence is
    *approximately equivalent* to the reference indicated, and is typically used for
    describing large regions in contexts where the use of an approximate sequence is
    inconsequential.
    """

    type: Literal[VRSTypes.DERIVED_SEQUENCE_EXPRESSION] = \
        VRSTypes.DERIVED_SEQUENCE_EXPRESSION
    location: SequenceLocation = Field(
        ...,
        description="The location from which the approximate sequence is derived"
    )
    reverse_complement: StrictBool = Field(
        ...,
        description=("A flag indicating if the expressed sequence is the reverse "
                     "complement of the sequence referred to by `location`")
    )


class RepeatedSequenceExpression(BaseModelForbidExtra):
    """An expression of a sequence comprised of a tandem repeating subsequence."""

    type: Literal[VRSTypes.REPEATED_SEQUENCE_EXPRESSION] = \
        VRSTypes.REPEATED_SEQUENCE_EXPRESSION
    seq_expr: Union[DerivedSequenceExpression, LiteralSequenceExpression] = Field(
        ...,
        description="An expression of the repeating subsequence"
    )
    count: Union[DefiniteRange, IndefiniteRange, Number] = Field(
        ..., description="The count of repeated units, as an integer or inclusive range"
    )

    @validator("count")
    def check_count_value(cls, v):
        """Check count value"""
        if v.type in {VRSTypes.NUMBER, VRSTypes.INDEFINITE_RANGE}:
            assert v.value >= 0, "`count.value` minimum is 0"
        else:
            assert v.min >= 0 and v.max >= 0, "`count.min` and `count.max` minimum is 0"
        return v


class ComposedSequenceExpression(BaseModelForbidExtra):
    """An expression of a sequence composed from multiple other Sequence Expressions
    objects. MUST have at least one component that is not a LiteralSequenceExpression.
    CANNOT be composed from nested composed sequence expressions.
    """

    type: Literal[VRSTypes.COMPOSED_SEQUENCE_EXPRESSION] = \
        VRSTypes.COMPOSED_SEQUENCE_EXPRESSION
    components: List[
        Union[
            DerivedSequenceExpression,
            LiteralSequenceExpression,
            RepeatedSequenceExpression
        ]
    ] = Field(
        ...,
        description=("An ordered list of SequenceExpression components comprising "
                     "the expression."),
        min_items=2,
        unique_items=True
    )

    @validator("components", pre=False)
    def ensure_contains_rse_or_dse(cls, v):
        """Ensure that either RepeatedSequenceExpression or DerivedSequenceExpression"""
        e_types = [e.type for e in v
                   if e.type in {VRSTypes.REPEATED_SEQUENCE_EXPRESSION,
                                 VRSTypes.DERIVED_SEQUENCE_EXPRESSION}
                   ]
        assert e_types, ("`components` must contain either "
                         "`RepeatedSequenceExpression` or `DerivedSequenceExpression`")


class SequenceExpression(BaseModel):
    """An expression describing a Sequence."""

    __root__: Union[
        ComposedSequenceExpression,
        DerivedSequenceExpression,
        LiteralSequenceExpression,
        RepeatedSequenceExpression
    ] = Field(..., description="An expression describing a Sequence.")


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# DEPRECATED SequenceExpression


class SequenceState(BaseModelForbidExtra, BaseModelDeprecated):
    """DEPRECATED. A Sequence as a State. This is the State class to use for
    representing "ref-alt" style variation, including SNVs, MNVs, del, ins, and delins.
    This class is deprecated. Use LiteralSequenceExpression instead.
    """

    type: Literal[VRSTypes.SEQUENCE_STATE] = VRSTypes.SEQUENCE_STATE
    sequence: Sequence = Field(..., description="A string of Residues", example="C")

    _replace_with = "LiteralSequenceExpression"

    _get_sequence_val = validator("sequence", allow_reuse=True)(return_value)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Feature


class Gene(BaseModelForbidExtra):
    """A reference to a Gene as defined by an authority. For human genes, the use of
    `hgnc <https://registry.identifiers.org/registry/hgnc>`_ as the gene authority is
    RECOMMENDED.
    """

    type: Literal[VRSTypes.GENE] = VRSTypes.GENE
    gene_id: CURIE = Field(..., description="A CURIE reference to a Gene concept")

    _get_gene_id_val = validator("gene_id", allow_reuse=True)(return_value)


class Feature(BaseModelForbidExtra):
    """A named entity that can be mapped to a Location. Genes, protein domains, exons,
    and chromosomes are some examples of common biological entities that may be
    Features.
    """

    __root__: Gene = Field(
        ...,
        description=("A named entity that can be mapped to a Location. Genes, protein "
                     "domains, exons, and chromosomes are some examples of common "
                     "biological entities that may be Features.")
    )

    class Config:
        """Configure Pydantic attributes."""

        @staticmethod
        def schema_extra(schema, model):
            """Ensure JSON schema output matches original VRS model."""
            if "$ref" in schema:
                del schema["$ref"]
            schema["anyOf"] = [{"$ref": "#/components/schemas/Gene"}]


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Molecular Variation


class Allele(BaseModelForbidExtra):
    """The state of a molecule at a Location."""

    id: Optional[CURIE] = Field(
        alias="_id",
        description="Variation Id. MUST be unique within document."
    )
    type: Literal[VRSTypes.ALLELE] = VRSTypes.ALLELE
    location: Union[
        CURIE,
        ChromosomeLocation,
        SequenceLocation
    ] = Field(..., description="Where Allele is located")
    state: Union[
        ComposedSequenceExpression,
        DerivedSequenceExpression,
        LiteralSequenceExpression,
        RepeatedSequenceExpression,
        SequenceState,
    ] = Field(..., description="An expression of the sequence state")

    _get_id_val = validator("id", allow_reuse=True)(return_value)
    _get_loc_val = validator("location", allow_reuse=True)(return_value)


class Haplotype(BaseModelForbidExtra):
    """A set of non-overlapping Allele members that co-occur on the same molecule."""

    id: Optional[CURIE] = Field(
        alias="_id",
        description="Variation Id. MUST be unique within document."
    )
    type: Literal[VRSTypes.HAPLOTYPE] = VRSTypes.HAPLOTYPE
    members: List[Union[Allele, CURIE]] = Field(
        ...,
        description=("List of Alleles, or references to Alleles, that comprise this "
                     "Haplotype."),
        min_items=2,
        unique_items=True
    )

    _get_id_val = validator("id", allow_reuse=True)(return_value)
    _get_members_val = validator("members", allow_reuse=True)(return_value)


class MolecularVariation(BaseModel):
    """A variation on a contiguous molecule."""

    __root__: Union[Allele, Haplotype] = Field(
        ..., description="A variation on a contiguous molecule."
    )


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# SystemicVariation


class CopyNumberChange(BaseModelForbidExtra):
    """An assessment of the copy number of a Location or a Feature within a system
    (e.g. genome, cell,  etc.) relative to a baseline ploidy.
    """

    id: Optional[CURIE] = Field(
        alias="_id",
        description="Variation Id. MUST be unique within document."
    )
    type: Literal[VRSTypes.COPY_NUMBER_CHANGE] = VRSTypes.COPY_NUMBER_CHANGE
    subject: Union[
        CURIE,
        ChromosomeLocation,
        Gene,
        SequenceLocation
    ] = Field(
        ...,
        description="A location for which the number of systemic copies is described."
    )
    copy_change: CopyChange = Field(
        ...,
        description=("MUST be one of 'efo:0030069' (complete genomic loss), "
                     "'efo:0020073' (high-level loss), 'efo:0030068' (low-level loss),"
                     " 'efo:0030067' (loss), 'efo:0030064' (regional base ploidy), "
                     "'efo:0030070' (gain), 'efo:0030071' (low-level gain), "
                     "'efo:0030072' (high-level gain).")
    )

    _get_id_val = validator("id", allow_reuse=True)(return_value)
    _get_subject_val = validator("subject", allow_reuse=True)(return_value)


class CopyNumberCount(BaseModelForbidExtra):
    """The absolute count of discrete copies of a Location or Feature, within a system
    (e.g. genome, cell, etc.).
    """

    id: Optional[CURIE] = Field(
        alias="_id",
        description="Variation Id. MUST be unique within document."
    )
    type: Literal[VRSTypes.COPY_NUMBER_COUNT] = VRSTypes.COPY_NUMBER_COUNT
    subject: Union[CURIE, ChromosomeLocation, Gene, SequenceLocation] = Field(
        ...,
        description="A location for which the number of systemic copies is described"
    )
    copies: Union[DefiniteRange, IndefiniteRange, Number] = Field(
        ...,
        description="The integral number of copies of the subject in a system"
    )

    _get_id_val = validator("id", allow_reuse=True)(return_value)
    _get_subject_val = validator("subject", allow_reuse=True)(return_value)


class GenotypeMember(BaseModelForbidExtra):
    """A class for expressing the count of a specific MolecularVariation present
    *in-trans* at a genomic locus represented by a Genotype.
    """

    type: Literal[VRSTypes.GENOTYPE_MEMBER] = VRSTypes.GENOTYPE_MEMBER
    count: Union[DefiniteRange, IndefiniteRange, Number] = Field(
        ...,
        description="The number of copies of the `variation` at a Genotype locus."
    )
    variation: Union[Allele, Haplotype] = Field(
        ...,
        description="A MolecularVariation at a Genotype locus."
    )


class Genotype(BaseModelForbidExtra):
    """A quantified set of MolecularVariation associated with a genomic locus."""

    id: Optional[CURIE] = Field(
        alias="_id",
        description="Variation Id. MUST be unique within document."
    )
    type: Literal[VRSTypes.GENOTYPE] = VRSTypes.GENOTYPE
    members: List[GenotypeMember] = Field(
        ...,
        description=("Each GenotypeMember in `members` describes a MolecularVariation "
                     "and the count of that variation at the locus."),
        min_items=1,
        unique_items=True
    )
    count: Union[DefiniteRange, IndefiniteRange, Number] = Field(
        ...,
        description=("The total number of copies of all MolecularVariation at this "
                     "locus, MUST be greater than or equal to the sum of "
                     "GenotypeMember copy counts and MUST be greater than or equal to "
                     "1. If greater than the total of GenotypeMember counts, this "
                     "field describes  additional MolecularVariation that exist but "
                     "are not  explicitly described.")
    )

    _get_id_val = validator("id", allow_reuse=True)(return_value)


class SystemicVariation(BaseModel):
    """A Variation of multiple molecules in the context of a system, e.g. a genome,
    sample, or homologous chromosomes.
    """

    __root__: Union[CopyNumberChange, CopyNumberCount, Genotype] = Field(
        ...,
        description=("A Variation of multiple molecules in the context of a system, "
                     "e.g. a genome, sample, or homologous chromosomes.")
    )


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# UtilityVariation

class Text(BaseModelForbidExtra):
    """A free-text definition of variation."""

    id: Optional[CURIE] = Field(
        alias="_id",
        description="Variation Id. MUST be unique within document."
    )
    type: Literal[VRSTypes.TEXT] = VRSTypes.TEXT
    definition: StrictStr = Field(
        ...,
        description=("A textual representation of variation not representable by "
                     "other subclasses of Variation.")
    )

    _get_id_val = validator("id", allow_reuse=True)(return_value)


class VariationSet(BaseModelForbidExtra):
    """An unconstrained set of Variation members."""

    id: Optional[CURIE] = Field(
        alias="_id",
        description="Variation Id. MUST be unique within document."
    )
    type: Literal[VRSTypes.VARIATION_SET] = VRSTypes.VARIATION_SET
    members: List[
        Union[
            Allele,
            CURIE,
            CopyNumberChange,
            CopyNumberCount,
            Genotype,
            Haplotype,
            Text,
            VariationSet
        ]
    ] = Field(
        ...,
        description=("List of Variation objects or identifiers. Attribute is "
                     "required, but MAY be empty.")
    )

    _get_id_val = validator("id", allow_reuse=True)(return_value)
    _get_members_val = validator("members", allow_reuse=True)(return_value)

    @validator("members", pre=False)
    def ensure_unique_items(cls, members):
        """Ensure items are unique"""
        if not members:
            return members

        for i, value in enumerate(members, start=1):
            if value in members[i:]:
                raise ValidationError("`members` must have unique items")

        return members


class UtilityVariation(BaseModel):
    """A collection of Variation subclasses that cannot be constrained to a specific
    class of biological variation, but are necessary for some applications of VRS.
    """

    __root__: Union[Text, VariationSet] = Field(
        ...,
        description=("A collection of Variation subclasses that cannot be constrained "
                     "to a specific class of biological variation, but are necessary "
                     "for some applications of VRS.")
    )


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Variation


class Variation(BaseModel):
    """A representation of the state of one or more biomolecules."""

    __root__: Union[
        Allele,
        CopyNumberChange,
        CopyNumberCount,
        Genotype,
        Haplotype,
        Text,
        VariationSet,
    ] = Field(
        ..., description="A representation of the state of one or more biomolecules."
    )


VariationSet.update_forward_refs()
