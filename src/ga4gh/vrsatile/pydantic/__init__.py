"""Initialize GA4GH VRSATILE Pydantic."""
import logging
from abc import ABC

from pydantic import BaseModel, Extra


logger = logging.getLogger("ga4gh.vrsatile.pydantic")


class BaseModelForbidExtra(BaseModel, ABC):
    """Base Pydantic model class with extra values forbidden."""

    class Config:
        """Class configs."""

        extra = Extra.forbid


def return_value(cls, v):
    """Return value from object.

    :param ModelMetaclass cls: Pydantic Model ModelMetaclass
    :param v: Model from vrs or vrsatile
    :return: Value
    """
    if v:
        try:
            if hasattr(v, "__root__"):
                v = return_value(cls, v.__root__)
            elif isinstance(v, list):
                tmp = list()
                for item in v:
                    while True:
                        try:
                            item = item.__root__
                        except AttributeError:
                            break
                    tmp.append(item)
                v = tmp
        except AttributeError:
            pass
    return v
