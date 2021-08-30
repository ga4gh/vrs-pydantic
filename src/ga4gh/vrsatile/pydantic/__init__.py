"""Initialize GA4GH VRSATILE Pydantic."""


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
                    tmp.append(item.__root__)
                v = tmp
            else:
                v = v.__root__
        except AttributeError:
            pass
    return v
