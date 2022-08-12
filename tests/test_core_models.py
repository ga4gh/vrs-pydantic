"""Module for testing CORE models"""
import pytest

import pydantic


from ga4gh.vrsatile.pydantic.core_models import CURIE, DomainEntity, ExtensibleEntity, \
    Extension, Entity, RecordMetadata, ValueEntity, Coding, Disease, Phenotype,\
    Gene, Condition


@pytest.fixture(scope="module")
def value_entity():
    """Create test fixture for ValueEntity"""
    return ValueEntity(id="value:id", type="ValueEntity")


@pytest.fixture(scope="module")
def extensible_entity(extension):
    """Create test fixture for ExtensibleEntity"""
    return ExtensibleEntity(id="value:id", label="label", extensions=[extension],
                            type="ExtensibleEntity")


@pytest.fixture(scope="module")
def record_metadata():
    """Create test fixture for RecordMetadata"""
    return RecordMetadata(id="value:id", type="RecordMetadata",
                          is_version_of="version:test", version="1.0.0",
                          label="Test RecordMetadata")


def test_curie():
    """Test that the CURIE model works correctly"""
    curie = CURIE(__root__="valid:curie")
    assert curie.__root__ == "valid:curie"

    with pytest.raises(pydantic.error_wrappers.ValidationError):
        CURIE(__root__="invalid")


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


def test_value_entity(value_entity):
    """Test that ValueEntity model works correctly"""
    def _value_entity_checks(ve):
        assert ve.id == "value:id"
        assert ve.type == "ValueEntity"
    _value_entity_checks(value_entity)

    value_entity = ValueEntity(**{"id": "value:id", "type": "ValueEntity"})
    _value_entity_checks(value_entity)

    invalid_params = [
        {"id": "invalid", "type": "ValueEntity"},
        {"id": "value:id"}
    ]

    for invalid_param in invalid_params:
        with pytest.raises(pydantic.error_wrappers.ValidationError):
            ValueEntity(**invalid_param)


def test_extensible_entity(extensible_entity, extension):
    """Test that ExtensibleEntity model works correctly"""
    def _extensible_entity_checks(ee):
        assert ee.id == "value:id"
        assert ee.label == "label"
        assert len(ee.extensions) == 1
        assert ee.extensions[0].name == "name"
        assert ee.extensions[0].value == ["value1", "value2"]
        assert ee.type == "ExtensibleEntity"

    _extensible_entity_checks(extensible_entity)

    extensible_entity = ExtensibleEntity(**{"id": "value:id", "label": "label",
                                            "extensions": [extension],
                                            "type": "ExtensibleEntity"})
    _extensible_entity_checks(extensible_entity)

    invalid_params = [
        {"id": "value"},
        {"type": "ExtensibleEntity", "extensions": ["label"]}
    ]
    for invalid_param in invalid_params:
        with pytest.raises(pydantic.error_wrappers.ValidationError):
            ExtensibleEntity(**invalid_param)


def test_entity(extensible_entity):
    """Test that Entity model works correctly"""
    entity = Entity(__root__=extensible_entity)
    assert entity.__root__ == extensible_entity

    with pytest.raises(pydantic.error_wrappers.ValidationError):
        Entity(__root__="fail")


def test_domain_entity():
    """Test that DomainEntity model works correctly"""
    def _domain_entity_checks(de):
        assert de.id == "value:id"
        assert de.type == "DomainEntity"

    domain_entity = DomainEntity(id="value:id", type="DomainEntity")
    _domain_entity_checks(domain_entity)

    domain_entity = DomainEntity(**{"id": "value:id", "type": "DomainEntity"})
    _domain_entity_checks(domain_entity)

    invalid_params = [
        {"id": "invalid", "type": "DomainEntity"},
        {"id": "value:id"}
    ]

    for invalid_param in invalid_params:
        with pytest.raises(pydantic.error_wrappers.ValidationError):
            DomainEntity(**invalid_param)


def test_record_metadata(record_metadata):
    """Test that RecordMetadata model works correctly"""
    def _record_metadata_checks(rm):
        assert rm.id == "value:id"
        assert rm.type == "RecordMetadata"
        assert rm.is_version_of == "version:test"
        assert rm.version == "1.0.0"
        assert rm.label == "Test RecordMetadata"

    _record_metadata_checks(record_metadata)

    record_metadata = RecordMetadata(**{"id": "value:id", "type": "RecordMetadata",
                                        "is_version_of": "version:test",
                                        "version": "1.0.0",
                                        "label": "Test RecordMetadata"})
    _record_metadata_checks(record_metadata)

    invalid_params = [
        {"type": "RecordMetadata", "version": 2},
        {"type": "RecordMetadata", "is_version_of": "invalid"}
    ]
    for invalid_param in invalid_params:
        with pytest.raises(pydantic.error_wrappers.ValidationError):
            RecordMetadata(**invalid_param)


def test_coding(record_metadata, extension):
    """Test that Coding model works correctly"""
    def _coding_checks(c):
        assert c.id == "value:id"
        assert c.type == "Coding"
        assert c.record_metadata == record_metadata
        assert c.label == "coding"

    coding = Coding(id="value:id", type="Coding", record_metadata=record_metadata,
                    label="coding")
    _coding_checks(coding)

    coding = Coding(**{"id": "value:id", "type": "Coding", "label": "coding",
                       "record_metadata": record_metadata})
    _coding_checks(coding)

    invalid_params = [
        {"id": "test", "type": "Coding"},
        {"id": "value:id", "type": "Coding", "record_metadata": extension}
    ]
    for invalid_param in invalid_params:
        with pytest.raises(pydantic.error_wrappers.ValidationError):
            Coding(**invalid_param)


def test_disease(disease):
    """Test that Disease model works correctly"""
    def _disease_checks(d):
        assert d.id == "ncit:C4989"
        assert d.type == "Disease"

    _disease_checks(disease)

    disease = Disease(**{"id": "ncit:C4989", "type": "Disease"})
    _disease_checks

    with pytest.raises(pydantic.error_wrappers.ValidationError):
        Disease(**{"id": "invalid"})


def test_phenotype(phenotype):
    """Test that Phenotype model works correctly."""
    def _phenotype_checks(p):
        assert p.id == "HP:0000002"
        assert p.type == "Phenotype"
    _phenotype_checks(phenotype)

    phenotype = Phenotype(**{"id": "HP:0000002", "type": "Phenotype"})
    _phenotype_checks(phenotype)

    with pytest.raises(pydantic.error_wrappers.ValidationError):
        Phenotype(**{"id": "invalid"})


def test_gene():
    """Test that Gene model works correctly."""
    g = Gene(id="hgnc:5")
    assert g.id == "hgnc:5"
    assert g.type == "Gene"

    assert Gene(id="hgnc:5", type="Gene")

    invalid_params = [
        {"id": "hgnc5"},
        {"id": "test:1", "type": "GeneDescriptor"}
    ]

    for invalid_param in invalid_params:
        with pytest.raises(pydantic.error_wrappers.ValidationError):
            Gene(**invalid_param)


def test_condition(phenotype, disease, gene):
    """Test that Condition model works correctly"""

    def _condition_checks(c):
        assert c.members[0] == disease
        assert c.members[1] == phenotype
        assert c.type == "Condition"

    condition = Condition(members=[disease, phenotype], type="Condition")
    _condition_checks(condition)

    condition = Condition(**{"members": [disease, phenotype], "type": "Condition"})
    _condition_checks(condition)

    invalid_params = [
        {"members": [phenotype], "type": "Condition"},
        {"members": [gene], "type": "Condition"}
    ]
    for invalid_param in invalid_params:
        with pytest.raises(pydantic.error_wrappers.ValidationError):
            Condition(**invalid_param)
