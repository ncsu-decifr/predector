"""
Prediction of protein subcellular localization
"""

try:
    import theano
except ImportError:
    raise ImportError("""Could not import Theano.
Please make sure you install a recent enough version of Theano.
""")
else:
    del theano

try:
    import lasagne
except ImportError:
    raise ImportError("""Could not import Lasagne.
Please make sure you install a recent enough version of Lasagne.
""")
else:
    del lasagne

try:
    import numpy
except ImportError:
    raise ImportError("""Could not import Numpy.
Please make sure you install a recent enough version of Numpy.
""")
else:
    del numpy


from . import models
from . import utils