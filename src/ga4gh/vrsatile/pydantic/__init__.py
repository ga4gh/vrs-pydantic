"""Initialize GA4GH VRSATILE Pydantic."""
import logging
from abc import ABC

from pydantic import BaseModel, Extra, root_validator


logger = logging.getLogger("FUSOR")


class BaseModelForbidExtra(BaseModel, ABC):
    """Base Pydantic model class with extra values forbidden."""

    class Config:
        """Class configs."""

        extra = Extra.forbid


class BaseModelDeprecated(BaseModel, ABC):
    """Base Pydantic model class to use for deprecated classes."""

    @root_validator(pre=True)
    def log_deprecated_warning(cls, values):
        """Log warning that object class is deprecated."""
        logger.warning(f"Using deprecated object: {cls.__name__}")
        return values


def return_value(cls, v):
    """Return value from object.

    :param ModelMetaclass cls: Pydantic Model ModelMetaclass
    :param v: Model from vrs or vrsatile
    :return: Value
    """
    if v is not None:
        try:
            if isinstance(v, list):
                tmp = list()
                for item in v:
                    while True:
                        try:
                            item = item.__root__
                        except AttributeError:
                            break
                    tmp.append(item)
                v = tmp
            else:
                v = v.__root__
        except AttributeError:
            pass
    return v
