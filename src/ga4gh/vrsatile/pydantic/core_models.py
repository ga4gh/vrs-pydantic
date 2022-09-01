"""Define models from core."""
from typing import Any, List, Literal, Optional, Union

from pydantic import StrictStr, Field, validator, constr

from ga4gh.vrsatile.pydantic import return_value, BaseModel


class CURIE(BaseModel):
    """A `W3C Compact URI <https://www.w3.org/TR/curie/>` formatted string.  A CURIE
    string has the structure ``prefix``:``reference``, as defined by the W3C syntax.
    """

    __root__: constr(regex=r"^\w[^:]*:.+$") = \
        Field(..., example="ensembl:ENSG00000139618")  # noqa: F722


class Extension(BaseModel):
    """The Extension class provides VODs with a means to extend descriptions with other
    attributes unique to a content provider. These extensions are not expected to be
    natively understood under VRSATILE, but may be used for pre-negotiated exchange of
    message attributes when needed.
    """

    type: Literal["Extension"] = "Extension"
    name: StrictStr
    value: Optional[Any]


class EntityBaseModel(BaseModel):
    """Create base model for Entity model"""

    id: Optional[StrictStr]
    type: StrictStr


class ValueEntity(EntityBaseModel):
    """ValueEntity is the root class for classes that instantiate Value Objects.
    ValueEntity classes are not extensible and MUST NOT have optional properties.
    """

    id: Optional[CURIE]

    _get_id_val = validator('id', allow_reuse=True)(return_value)


class ExtensibleEntity(EntityBaseModel):
    """ExtensibleEntity is the root class for classes that instantiate Extensible
    Objects. Extensible Objects are extensible using the extensions property and MAY
    have optional properties.
    """

    label: Optional[StrictStr]
    extensions: Optional[List[Extension]]


class Entity(BaseModel):
    """Entity is the root class of `core` classes model - those that have identifiers
    and other general metadata like labels, xrefs, urls, descriptions, etc. All core
    classes descend from and inherit its attributes.
    """

    __root__: Union[ExtensibleEntity, ValueEntity]


class DomainEntity(ValueEntity):
    """An abstract :ref:`ValueEntity` class extended to capture specific domain entities
    by reference to an external identifier.
    """

    id: CURIE

    _get_id_val = validator('id', allow_reuse=True)(return_value)


class RecordMetadata(ExtensibleEntity):
    """A re-usable structure that encapsulates provenance metadata that applies to a
    specific concrete record of information as encoded in a particular system, as
    opposed to provenance of the abstract information content/knowledge the record
    represents.
    """

    type: Literal["RecordMetadata"] = "RecordMetadata"
    is_version_of: Optional[CURIE]
    version: Optional[StrictStr]

    _get_is_version_of_val = validator('is_version_of', allow_reuse=True)(return_value)


class Coding(ExtensibleEntity):
    """A `coding` is an extensible entity for labeling or otherwise annotating globally
    namespaced identifiers known as "codes".
    """

    type: Literal["Coding"] = "Coding"
    id: Optional[CURIE]
    record_metadata: Optional[RecordMetadata]

    _get_id_val = validator('id', allow_reuse=True)(return_value)


class Disease(DomainEntity):
    """A reference to a Disease as defined by an authority. For human diseases, the
    use of `MONDO <https://registry.identifiers.org/registry/mondo>` as the disease
    authority is RECOMMENDED.
    """

    type: Literal["Disease"] = "Disease"


class Phenotype(DomainEntity):
    """A reference to a Phenotype as defined by an authority. For human phenotypes, the
    use of `HPO <https://registry.identifiers.org/registry/hpo> as the disease authority
    is RECOMMENDED.
    """

    type: Literal["Phenotype"] = "Phenotype"


class Gene(DomainEntity):
    """A reference to a Gene as defined by an authority. For human genes, the use of
    `hgnc <https://registry.identifiers.org/registry/hgnc>` as the gene authority is
    RECOMMENDED.
    """

    type: Literal["Gene"] = "Gene"


class Condition(ValueEntity):
    """A set of phenotype and/or disease concepts that constitute a condition."""

    members: List[Union[Disease, Phenotype]] = Field(..., min_items=2)


class Therapeutic(DomainEntity):
    """A treatment, therapy, or drug"""

    type: Literal["Therapeutic"] = "Therapeutic"


class TherapeuticCollection(ValueEntity):
    """A collection of therapeutics."""

    members: List[Therapeutic] = Field(..., min_items=2)


class CombinationTherapeuticCollection(TherapeuticCollection):
    """A collection of therapeutics that are taken during a course of treatment."""

    type: Literal["CombinationTherapeutics"] = "CombinationTherapeutics"


class SubstituteTherapeuticCollection(TherapeuticCollection):
    """A collection of therapeutics that are considered as valid alternative
    entities.
    """

    type: Literal["SubstituteTherapeutics"] = "SubstituteTherapeutics"
