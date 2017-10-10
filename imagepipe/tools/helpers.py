import collections
import os
import types

import numpy as np
from PIL import Image


def safe_dir_create(path):
    if not os.path.isdir(path):
        os.makedirs(path)


class PipeArgError(ValueError):
    pass


def list_not_string(argument):
    """
    function that checks if a list is not a string

    credits: http://stackoverflow.com/questions/1055360/how-to-tell-a-variable-is-iterable-but-not-a-string

    :param argument:
    :return:
    """
    if isinstance(argument, collections.Iterable):
        if isinstance(argument, types.StringTypes):
            return False
        else:
            return True
    else:
        raise PipeArgError("Expected a name of channel or list of names. Found: '%s' " % argument)


def tiff_stack_2_np_arr(tiff_location):
    """
    Loads the image from the tiff stack to a 3D numpy array

    :param tiff_location:
    :return:
    """
    tiff_stack = Image.open(tiff_location)
    stack = [np.array(tiff_stack)]
    try:
        while 1:
            tiff_stack.seek(tiff_stack.tell() + 1)
            stack.append(np.array(tiff_stack))
    except EOFError:
        pass

    return np.array(stack)

