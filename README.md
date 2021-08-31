# vrsatile-pydantic
Translation of the GA4GH [VRS v1.2.0](https://vrs.ga4gh.org/en/1.2.0) and [VRSATILE](https://vrsatile.readthedocs.io/en/latest/) schemas to a Pydantic data model

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
