"""
I am sorry you are seeing this, but yes naughty monkeypatching is happening here.
Colab still runs on 3.7...
"""
import sys, typing_extensions, typing, functools
from singledispatchmethod import singledispatchmethod  # noqa it's okay, PyCharm, I am not a technoluddite


# hack to enable the backport:
def backport():
    if sys.version_info.major != 3:
        raise NotImplementedError('This is only for Python 3')
    if sys.version_info.minor > 8:
        return
    functools.singledispatchmethod = singledispatchmethod
    for key, fun in typing_extensions.__dict__.items():
      if key not in typing.__dict__:
        setattr(typing, key, fun)
