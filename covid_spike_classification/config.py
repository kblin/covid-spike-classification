"""Manage configuration options."""


class CSCConfig:
    __slots__ = (
        'debug',
        '_failed',
        'input_format',
        'outdir',
        'quiet',
        'reads',
        'reference',
        'show_unexpected',
        'silence_warnings',
        'stdout',
        'zip_results',
    )

    def __init__(self, **kwargs):
        for attr in self.__slots__:
            # skip internal slots
            if attr.startswith("_"):
                continue
            setattr(self, attr, kwargs[attr])

        # set up internal slots
        self._failed = set()

    @classmethod
    def from_args(cls, namespace):
        kwargs = {}
        for attr in cls.__slots__:
            if hasattr(namespace, attr):
                kwargs[attr] = getattr(namespace, attr)

        return cls(**kwargs)
