"""Theano compatibility shim using aesara"""
try:
    import aesara as _T
    from aesara import *
    from aesara import tensor, scan
    __version__ = _T.__version__
except ImportError:
    raise ImportError("Neither theano nor aesara is available")
