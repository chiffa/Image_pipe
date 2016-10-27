from decorator import decorator
import numpy as np
import os
from PIL import Image
from debugger import CustomDebugger
from scipy.ndimage.filters import gaussian_filter


debugger = CustomDebugger()


dtype2bits = {'uint8': 8,
              'uint16': 16,
              'uint32': 32}


class PipeArgError(ValueError):
    pass


@decorator
def generator_wrapper(f, *args, **kwargs):
    """
    converts a function to accepting a generator of named dicts and adds channel selection logic
    """
    iterator = args[0]
    args = args[1:]

    for named_dict in iterator:
        local_args = (named_dict,) + args
        yield f(*local_args, **kwargs)


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


def split_and_trim(prefix, main_root):
    trim_length = len(main_root)
    if main_root[-1] != os.sep:
        trim_length += 1

    return prefix[trim_length:].split(os.sep)


@generator_wrapper
def gamma_stabilize(current_image, alpha_clean=5):
    bits = dtype2bits[current_image.dtype.name]
    stabilized = (current_image - np.min(current_image))/(float(2**bits) - np.min(current_image))
    stabilized[stabilized < alpha_clean*np.median(stabilized)] = 0

    return stabilized


@generator_wrapper
def smooth(current_image, smoothing_px=1.5):
    for i in range(0, current_image.shape[0]):
            current_image[i, :, :] = gaussian_filter(current_image[i, :, :],
                                                     smoothing_px, mode='constant')
            current_image[current_image < 5*np.mean(current_image)] = 0

    return current_image


@generator_wrapper
def extract(name_dict, channel):
    return name_dict[channel]

