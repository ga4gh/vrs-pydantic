# vrsatile-pydantic
Translation of the GA4GH [VRS](https://vrs.ga4gh.org/en/stable/) and [VRSATILE](https://vrsatile.readthedocs.io/en/latest/) schemas to a Pydantic data model

The ga4gh/vrsatile/pydantic repo depends on VRS and VRSATILE models, and therefore each ga4gh.vrsatile.pydantic package on PyPI uses a particular version of VRS and VRSATILE. The correspondences between the packages may be summarized as:

| ga4gh.vrsatile.pydantic branch | ga4gh.vrsatile.pydantic version | VRS version | VRSATILE version |
| ---- | --- | --- | --- |
| vrs-1.2 *(no longer being updated)* | 0.0.X | [1.2.X](https://github.com/ga4gh/vrs/tree/1.2) | [main](https://github.com/ga4gh/vrsatile/tree/main) |
| metaschema-update | 0.1.X | [metaschema-update](https://github.com/ga4gh/vrs/tree/metaschema-update) | [metaschema-update](https://github.com/ga4gh/vrsatile/tree/metaschema-update)
| main | 0.2.X | [1.3.X](https://github.com/ga4gh/vrs/tree/1.3) | [main](https://github.com/ga4gh/vrsatile/tree/main)

## Installation

vrsatile-pydantic is available via PyPI:
```commandline
pip install ga4gh.vrsatile.pydantic
```

## Developer Instructions

Following are sections include instructions specifically for developers.

### Installation

For a development install, we recommend using Pipenv. See the
[pipenv docs](https://pipenv-fork.readthedocs.io/en/latest/#install-pipenv-today)
for direction on installing pipenv in your compute environment.

Once installed, clone the repo and initialize the environment:

```commandline
git clone https://github.com/ga4gh/vrsatile-pydantic
cd vrsatile-pydantic
pipenv shell
pipenv update
pipenv install --dev
```

### Init coding style tests

Code style is managed by [flake8](https://github.com/PyCQA/flake8) and checked prior to commit.

We use [pre-commit](https://pre-commit.com/#usage) to run conformance tests.

This ensures:

* Check code style
* Check for added large files
* Detect AWS Credentials
* Detect Private Key

Before first commit run:

```commandline
pre-commit install
```

### Running unit tests

Running unit tests is as easy as pytest.

```commandline
pipenv run pytest tests
```
