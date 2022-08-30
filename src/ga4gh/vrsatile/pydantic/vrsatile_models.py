"""Define Pydantic Class models for VRSATILE models."""
from enum import Enum
from typing import List, Optional, Union, Literal

from pydantic import StrictInt, StrictStr, validator, Field, BaseModel, Extra

from ga4gh.vrsatile.pydantic import return_value, BaseModelForbidExtra
from ga4gh.vrsatile.pydantic.core_models import CURIE, Therapeutic, ValueEntity, \
    CombinationTherapeuticCollection, Condition, Disease, ExtensibleEntity, Phenotype, \
    SubstituteTherapeuticCollection
from ga4gh.vrsatile.pydantic.vrs_models import Variation, SequenceLocation, \
    ChromosomeLocation, Sequence, Gene


class VODClassName(str, Enum):
    """Define VOD class names."""

    VARIATION_DESCRIPTOR = "VariationDescriptor"
    LOCATION_DESCRIPTOR = "LocationDescriptor"
    SEQUENCE_DESCRIPTOR = "SequenceDescriptor"
    GENE_DESCRIPTOR = "GeneDescriptor"
    CATEGORICAL_VARIATION_DESCRIPTOR = "CategoricalVariationDescriptor"
    CONDITION_DESCRIPTOR = "ConditionDescriptor"
    DISEASE_DESCRIPTOR = "DiseaseDescriptor"
    PHENOTYPE_DESCRIPTOR = "PhenotypeDescriptor"
    CANONICAL_VARIATION_DESCRIPTOR = "CanonicalVariationDescriptor"
    THERAPEUTIC_DESCRIPTOR = "TherapeuticDescriptor"
    THERAPEUTIC_COLLECTION_DESCRIPTOR = "TherapeuticsCollectionDescriptor"


class VRSATILETypes(str, Enum):
    """Define types used in VRSATILE."""

    EXPRESSION = "Expression"


class MoleculeContext(str, Enum):
    """The structural variant type associated with this variant."""

    GENOMIC = "genomic"
    TRANSCRIPT = "transcript"
    PROTEIN = "protein"


class ExpressionSyntax(str, Enum):
    """Possible values for the Expression `syntax` property."""

    HGVS_C = "hgvs.c"
    HGVS_P = "hgvs.p"
    HGVS_G = "hgvs.g"
    HGVS_M = "hgvs.m"
    HGVS_N = "hgvs.n"
    HGVS_R = "hgvs.r"
    ISCN = "iscn"
    GNOMAD = "gnomad"
    SPDI = "spdi"


class CategoricalVariationType(str, Enum):
    """Possible types for Categorical Variations."""

    CANONICAL_VARIATION = "CanonicalVariation"
    COMPLEX_VARIATION = "ComplexVariation"


class CategoricalVariationBase(ValueEntity):
    """Base class for Categorical Variation"""

    pass


class CanonicalVariation(CategoricalVariationBase):
    """A categorical variation domain characterized by a representative
    Variation context to which members lift-over, project, translate, or
    otherwise directly align.
    """

    type: Literal[CategoricalVariationType.CANONICAL_VARIATION] = \
        CategoricalVariationType.CANONICAL_VARIATION
    canonical_context: Variation

    _get_caononical_context_val = \
        validator('canonical_context', allow_reuse=True)(return_value)

    class Config:
        """Class configs."""

        extra = Extra.forbid


class ComplexVariationOperator(str, Enum):
    """Possible values for the Complex Variation's `operator` field."""

    AND = "AND"
    OR = "OR"


class ComplexVariation(CategoricalVariationBase):
    """A categorical variation domain jointly characterized by two or more other
    categorical variation domains.
    """

    type: Literal[CategoricalVariationType.COMPLEX_VARIATION] = \
        CategoricalVariationType.COMPLEX_VARIATION
    operands: List[CanonicalVariation] = Field(..., min_items=2)
    operator: ComplexVariationOperator

    class Config:
        """Class configs."""

        extra = Extra.forbid


class CategoricalVariation(BaseModel):
    """A representation of a categorically-defined `functional domain
    <https://en.wikipedia.org/wiki/Domain_of_a_function>`_ for variation, in
    which individual variation instances may be members.
    """

    __root__: Union[CanonicalVariation, ComplexVariation]


class Expression(BaseModelForbidExtra):
    """Representation of a variation by a specified nomenclature or syntax for
    a Variation object. Common examples of expressions for the description of
    molecular variation include the HGVS and ISCN nomenclatures.
    """

    type: Literal[VRSATILETypes.EXPRESSION] = VRSATILETypes.EXPRESSION
    syntax: ExpressionSyntax
    value: StrictStr
    syntax_version: Optional[StrictStr]


class VCFRecord(BaseModelForbidExtra):
    """This data class is used when it is desirable to pass data as expected
    from a VCF record.
    """

    genome_assembly: StrictStr
    chrom: StrictStr
    pos: StrictInt
    id: Optional[StrictStr]
    ref: StrictStr
    alt: StrictStr
    qual: Optional[StrictStr]
    filter: Optional[StrictStr]
    info: Optional[StrictStr]


class ValueObjectDescriptorBaseModel(ExtensibleEntity):
    """Define Value Object Descriptor Base Model class. Contains fields that ALL
    descriptors use.
    """

    label: Optional[StrictStr]
    description: Optional[StrictStr]
    xrefs: Optional[List[CURIE]]
    alternate_labels: Optional[List[StrictStr]]

    _get_xrefs_val = validator('xrefs', allow_reuse=True)(return_value)

    @validator("xrefs")
    def check_count_value(cls, v):
        """Check xrefs value"""
        if v:
            assert len(v) == len(set(v)), "xrefs must contain unique items"
        return v

    class Config:
        """Class configs."""

        extra = Extra.forbid


class ValueObjectDescriptor(ValueObjectDescriptorBaseModel):
    """The abstract *Value Object Descriptor* parent class. All attributes of this
    parent class are inherited by descendent classes.
    """

    value: Union[CURIE, ValueEntity]

    _get_value = validator('value', allow_reuse=True)(return_value)

    class Config:
        """Class configs."""

        extra = Extra.forbid


class SequenceDescriptor(ValueObjectDescriptorBaseModel):
    """This descriptor is intended to reference VRS Sequence value objects."""

    type: Literal[VODClassName.SEQUENCE_DESCRIPTOR] = VODClassName.SEQUENCE_DESCRIPTOR
    sequence: Union[CURIE, Sequence]
    residue_type: Optional[CURIE]

    _get_sequence_val = validator('sequence', allow_reuse=True)(return_value)
    _get_residue_type_val = validator('residue_type', allow_reuse=True)(return_value)

    class Config:
        """Class configs."""

        extra = Extra.forbid


class LocationDescriptor(ValueObjectDescriptorBaseModel):
    """This descriptor is intended to reference VRS Location value objects."""

    type: Literal[VODClassName.LOCATION_DESCRIPTOR] = VODClassName.LOCATION_DESCRIPTOR
    location: Union[CURIE, SequenceLocation, ChromosomeLocation]

    _get_location_val = validator('location', allow_reuse=True)(return_value)

    class Config:
        """Class configs."""

        extra = Extra.forbid


class GeneDescriptor(ValueObjectDescriptorBaseModel):
    """This descriptor is intended to reference VRS Gene value objects."""

    type: Literal[VODClassName.GENE_DESCRIPTOR] = VODClassName.GENE_DESCRIPTOR
    gene: Union[CURIE, Gene]

    _get_gene_val = validator('gene', allow_reuse=True)(return_value)

    class Config:
        """Class configs."""

        extra = Extra.forbid


class VariationDescriptor(ValueObjectDescriptorBaseModel):
    """This descriptor class is used for describing VRS Variation value objects."""

    type: Literal[VODClassName.VARIATION_DESCRIPTOR] = VODClassName.VARIATION_DESCRIPTOR
    variation: Union[CURIE, Variation]
    molecule_context: Optional[MoleculeContext]
    structural_type: Optional[CURIE]
    expressions: Optional[List[Expression]]
    vcf_record: Optional[VCFRecord]
    gene_context: Optional[Union[CURIE, GeneDescriptor]]
    vrs_ref_allele_seq: Optional[Sequence]
    allelic_state: Optional[CURIE]

    _get_variation_val = validator('variation', allow_reuse=True)(return_value)
    _get_structural_type_val = \
        validator('structural_type', allow_reuse=True)(return_value)
    _get_gene_context_val = validator('gene_context', allow_reuse=True)(return_value)
    _get_vrs_allele_ref_seq_val = \
        validator('vrs_ref_allele_seq', allow_reuse=True)(return_value)
    _get_allelic_state_val = validator('allelic_state', allow_reuse=True)(return_value)

    class Config:
        """Class configs."""

        extra = Extra.forbid


class DiseaseDescriptor(ValueObjectDescriptorBaseModel):
    """This descriptor class is used for describing Disease domain entities."""

    type: Literal[VODClassName.DISEASE_DESCRIPTOR] = VODClassName.DISEASE_DESCRIPTOR
    disease: Union[CURIE, Disease]

    _get_disease_val = validator('disease', allow_reuse=True)(return_value)

    class Config:
        """Class configs."""

        extra = Extra.forbid


class PhenotypeDescriptor(ValueObjectDescriptorBaseModel):
    """This descriptor class is used for describing Phenotype domain entities."""

    type: Literal[VODClassName.PHENOTYPE_DESCRIPTOR] = VODClassName.PHENOTYPE_DESCRIPTOR
    phenotype: Union[CURIE, Phenotype]

    _get_phenotype_val = validator('phenotype', allow_reuse=True)(return_value)

    class Config:
        """Class configs."""

        extra = Extra.forbid


class ConditionDescriptor(ValueObjectDescriptorBaseModel):
    """This descriptor class is used for describing Condition entities."""

    type: Literal[VODClassName.CONDITION_DESCRIPTOR] = VODClassName.CONDITION_DESCRIPTOR
    condition: Union[CURIE, Condition]
    member_descriptors: List[Union[DiseaseDescriptor, PhenotypeDescriptor]]

    _get_condition_val = validator('condition', allow_reuse=True)(return_value)

    class Config:
        """Class configs."""

        extra = Extra.forbid


class TherapeuticDescriptor(ValueObjectDescriptorBaseModel):
    """This descriptor class is used for describing Therapeutic domain entities."""

    type: Literal[VODClassName.THERAPEUTIC_DESCRIPTOR] = \
        VODClassName.THERAPEUTIC_DESCRIPTOR
    therapeutic: Union[CURIE, Therapeutic]

    _get_therapeutic_val = validator("therapeutic", allow_reuse=True)(return_value)

    class Config:
        """Class configs."""

        extra = Extra.forbid


class TherapeuticCollectionDescriptor(ValueObjectDescriptorBaseModel):
    """This descriptor class is used for describing TherapeuticsCollection domain
    entities.
    """

    type: Literal[VODClassName.THERAPEUTIC_COLLECTION_DESCRIPTOR] = \
        VODClassName.THERAPEUTIC_COLLECTION_DESCRIPTOR
    therapeutic_collection: Union[CURIE, CombinationTherapeuticCollection,
                                  SubstituteTherapeuticCollection]
    member_descriptors: Optional[List[TherapeuticDescriptor]]

    _get_therapeutic_collection_val = validator("therapeutic_collection",
                                                allow_reuse=True)(return_value)

    class Config:
        """Class configs."""

        extra = Extra.forbid


class VariationMember(BaseModelForbidExtra):
    """A compact class for representing a variation context that is a member of a
    Categorical Variation. It supports one or more Expressions of a Variation and
    optionally an associated VRS ID.
    """

    type: Literal["VariationMember"] = "VariationMember"
    expressions: List[Expression] = Field(..., min_items=1)
    variation_id: Optional[CURIE]

    _get_variation_id_val = validator("variation_id", allow_reuse=True)(return_value)

    class Config:
        """Class configs."""

        extra = Extra.forbid


class CategoricalVariationDescriptorBaseModel(ValueObjectDescriptorBaseModel):
    """Base model for categorical variation descriptor"""

    type: StrictStr
    members: Optional[List[VariationMember]]


class CategoricalVariationDescriptor(CategoricalVariationDescriptorBaseModel):
    """This descriptor class is used for describing Categorical Variation value
    objects.
    """

    type: Literal[VODClassName.CATEGORICAL_VARIATION_DESCRIPTOR] = \
        VODClassName.CATEGORICAL_VARIATION_DESCRIPTOR
    categorical_variation: Union[CURIE, CanonicalVariation, ComplexVariation]

    _get_categorical_variation_val = \
        validator("categorical_variation", allow_reuse=True)(return_value)

    class Config:
        """Class configs."""

        extra = Extra.forbid


class CanonicalVariationDescriptor(CategoricalVariationDescriptorBaseModel):
    """This descriptor class is used for describing Canonical Variation value
    objects.
    """

    type: Literal[VODClassName.CANONICAL_VARIATION_DESCRIPTOR] = \
        VODClassName.CANONICAL_VARIATION_DESCRIPTOR
    subject_variation_descriptor: VariationDescriptor
    canonical_variation: Union[CURIE, CanonicalVariation]

    _get_canonical_variation_val = \
        validator("canonical_variation", allow_reuse=True)(return_value)

    class Config:
        """Class configs."""

        extra = Extra.forbid
