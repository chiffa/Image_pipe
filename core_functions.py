import numpy as np
import os
from PIL import Image
from debugger import CustomDebugger

debugger = CustomDebugger()


def tiff_stack_2_np_arr(tiff_location):
    """
    Loads the image from the tiff stack to a 3D numpy array

    :param tiff_stack:
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


def split_and_trim(prefix, main_root):
    trim_length = len(main_root)
    if main_root[-1] != os.sep:
        trim_length += 1
    return prefix[trim_length:].split(os.sep)

