# vrsatile-pydantic
Translation of the GA4GH [VRS](https://vrs.ga4gh.org/en/stable/) and [VRSATILE](https://vrsatile.readthedocs.io/en/latest/) schemas to a Pydantic data model

The ga4gh/vrsatile/pydantic repo depends on VRS and VRSATILE models, and therefore each ga4gh.vrsatile.pydantic package on PyPI uses a particular version of VRS and VRSATILE. The correspondences between the packages may be summarized as:

- **0.0.X ~ 1.2.X**: ga4gh.vrsatile.pydantic 0.0.X tracks [VRS 1.2.X](https://vrs.ga4gh.org/en/1.2.1/) and the [VRSATILE latest](https://vrsatile.readthedocs.io/en/latest/)
- **0.1.X ~ metaschema-updates**: ga4gh.vrsatile.pydantic 0.1.X tracks VRS and VRSATILE metaschema-updates

# Developer instructions

To install vrstaile-pydantic:
```commandline
pip install ga4gh.vrsatile.pydantic
```

Following are sections include instructions specifically for developers.

For a development install, we recommend using Pipenv. See the
[pipenv docs](https://pipenv-fork.readthedocs.io/en/latest/#install-pipenv-today)
for direction on installing pipenv in your compute environment.

Once installed, from the project root dir, just run:

```commandline
pipenv lock
pipenv sync
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
pipenv run pytest
```
