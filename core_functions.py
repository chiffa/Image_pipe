import numpy as np
import os
from PIL import Image
from debugger import CustomDebugger
import traceback
from matplotlib import pyplot as plt
from scipy.ndimage.filters import gaussian_filter
from scipy import stats
from collections import defaultdict
from csv import writer
from skimage.segmentation import random_walker
from skimage.morphology import closing
from scipy import ndimage as ndi
from skimage.morphology import watershed
from skimage.feature import peak_local_max
from skimage.morphology import disk
from skimage.morphology import skeletonize, medial_axis
from skimage.filters import threshold_otsu
from scipy.stats import t
from scipy.stats import ttest_ind
from itertools import combinations
from functools import wraps
import types
import collections

debugger = CustomDebugger()


dtype2bits = {'uint8': 8,
              'uint16': 16,
              'uint32': 32}

class PipeArgError(ValueError):
    pass


def doublewrap(f):
    """
    a decorator decorator, allowing the decorator to be used as:
    @decorator(with, arguments, and=kwargs)
    or
    @decorator

    credits: http://stackoverflow.com/questions/653368/how-to-create-a-python-decorator-that-can-be-used-either-with-or-without-paramet

    """
    @wraps(f)
    def new_dec(*args, **kwargs):
        if len(args) == 1 and len(kwargs) == 0 and callable(args[0]):
            # actual decorated function
            return f(args[0])
        else:
            # decorator arguments
            return lambda realf: f(realf, *args, **kwargs)

    return new_dec


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


@doublewrap
def generator_wrapper(f, expected_dims=(3,)):

    @wraps(f)
    def inner_wrapper(*args, **kwargs):
        """
        converts a function to accepting a generator of named dicts and adds channel selection logic
        """

        iterator = args[0]
        args = args[1:]

        if 'in_channel' in kwargs:
            in_chan = kwargs['in_channel']  # Multiple arguments
            del kwargs['in_channel']

            out_chan = ''
            log_chan = ''

            if 'out_channel' in kwargs:
                out_chan = kwargs['out_channel']  # Single output
                del kwargs['out_channel']

                if list_not_string(out_chan) and len(out_chan) != 1:
                    raise PipeArgError('Only one output channel can be bound')

            if 'log_channel' in kwargs:
                log_chan = kwargs['log_channel']
                del kwargs['log_channel']

            if not list_not_string(in_chan):  # convert string to list of strings
                if not out_chan:
                    out_chan = in_chan
                in_chan = [in_chan]

            if len(in_chan) != len(expected_dims):
                print in_chan, expected_dims
                print len(in_chan), len(expected_dims)
                raise PipeArgError('More channels are piped than function allows')

            if not out_chan:
                raise PipeArgError('Please provide out_channel argument')

            for named_dict in iterator:
                args_puck = []
                for i, chan in enumerate(in_chan):
                    if len(named_dict[chan].shape) != expected_dims[i]:
                        raise PipeArgError('Mismatched channel dimension for channel. %s is of dim %s, expected %s' %
                                           chan, len(chan.shape), expected_dims[i])
                    args_puck.append(named_dict[chan])
                local_args = tuple(args_puck) + args
                if log_chan:
                    kwargs['log_to'] = (named_dict, log_chan)
                named_dict[out_chan] = f(*local_args, **kwargs)
                yield named_dict

        else:
            for named_dict in iterator:
                local_args = (named_dict,) + args
                yield f(*local_args, **kwargs)

    return inner_wrapper


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


def _3d_stack_2d_filter(_3d_stack, _2d_filter):
    new_stack = np.zeros_like(_3d_stack)
    new_stack[:, _2d_filter] = _3d_stack[:, _2d_filter]
    return new_stack


def _2d_stack_2d_filter(_2d_stack, _2d_filter):
    new_stack = np.zeros_like(_2d_stack)
    new_stack[_2d_filter] = _2d_stack[_2d_filter]
    return new_stack


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
def sum_projection(current_image):
    return np.sum(current_image, axis=0)


@generator_wrapper(expected_dims=(2,))
def segment_out_cells(current_image):

    # TODO: parameters that might need to be eventually factored in:
    # upper/ lower bounds
    # min distance between cores
    # abs threshold
    # closing_size, fusing size

    gfp_clustering_markers = np.zeros(current_image.shape, dtype=np.uint8)
    # random walker segment
    gfp_clustering_markers[current_image > np.mean(current_image) * 2] = 2
    gfp_clustering_markers[current_image < np.mean(current_image) * 0.20] = 1
    labels = random_walker(current_image, gfp_clustering_markers, beta=10, mode='bf')
    # round up the labels and set the background to 0 from 1
    sel_elem = disk(2)
    labels = closing(labels, sel_elem)
    labels -= 1
    # prepare distances for the watershed
    distance = ndi.distance_transform_edt(labels)
    local_maxi = peak_local_max(distance,
                                indices=False,  # we want the image mask, not peak position
                                min_distance=10,  # about half of a bud with our size
                                threshold_abs=10,  # allows to clear the noise
                                labels=labels)
    # we fuse the labels that are close together that escaped the min distance in local_maxi
    local_maxi = ndi.convolve(local_maxi, np.ones((5, 5)), mode='constant', cval=0.0)
    # finish the watershed
    expanded_maxi_markers = ndi.label(local_maxi, structure=np.ones((3, 3)))[0]
    segmented_cells_labels = watershed(-distance, expanded_maxi_markers, mask=labels)

    # # log debugging data
    # running_debug_frame.current_image = current_image
    # running_debug_frame.gfp_clustering_markers = gfp_clustering_markers
    # running_debug_frame.labels = labels
    # running_debug_frame.segmented_cells_labels = segmented_cells_labels
    #
    return segmented_cells_labels


@generator_wrapper(expected_dims=(2,))
def qualifying_gfp(max_sum_projection):
    return max_sum_projection > np.median(max_sum_projection[max_sum_projection > 0])


@generator_wrapper(expected_dims=(2, 2, 2))
def aq_gfp_per_region(cell_labels, max_sum_projection, qualifying_gfp_mask):

    cells_average_gfp_list = []

    for i in range(1, np.max(cell_labels) + 1):

        current_mask = cell_labels == i
        current_cell_gfp = max_sum_projection[np.logical_and(current_mask, qualifying_gfp_mask)]

        if len(current_cell_gfp) == 0:
            continue

        gfp_percentile = np.percentile(current_cell_gfp, 50)
        gfp_average = np.average(max_sum_projection[np.logical_and(current_mask, max_sum_projection > gfp_percentile)])
        cells_average_gfp_list.append(gfp_average)

    return np.array(cells_average_gfp_list)


@generator_wrapper(expected_dims=(1,))
def detect_upper_outliers(cells_average_gfp_list, log_to=None):
    arg_sort = np.argsort(np.array(cells_average_gfp_list))
    cells_average_gfp_list = sorted(cells_average_gfp_list)
    cell_no = range(0, len(cells_average_gfp_list))
    # Non-trivial logic selecting the regression basis
    regression_base = min(len(cells_average_gfp_list) - 5, 10)
    slope, intercept, _, _, _ = stats.linregress(np.array(cell_no)[1:regression_base],
                                                 np.array(cells_average_gfp_list)[1:regression_base])
    std_err = (np.max(np.array(cells_average_gfp_list)[1:regression_base]) -
               np.min(np.array(cells_average_gfp_list)[1:regression_base])) / 2
    std_err *= 8
    predicted_average_gfp = intercept + slope * np.array(cell_no)

    upper_outliers = arg_sort[np.array(cell_no)[np.array(predicted_average_gfp + std_err) <
                                                 np.array(cells_average_gfp_list)]]

    # TODO: redirect the injection to a logging value in the pipe
    # => side_store argument
    embedded_dict = {'average area value': cells_average_gfp_list,
                     'predicted area value': predicted_average_gfp,
                     'classfication bounds': std_err
                     }
    if log_to:
        log_to[0][log_to[1]] = embedded_dict

    return upper_outliers


@generator_wrapper(expected_dims=(2, 1))
def paint_mask(label_masks, labels_to_paint):
    mask_to_paint = np.zeros_like(label_masks).astype(np.uint8)

    if labels_to_paint.tolist() != np.array([]):
        for idx in labels_to_paint.tolist():
            mask_to_paint[label_masks == idx + 1] = 1  # indexing starts from 1, not 0 for the labels

    return mask_to_paint


@generator_wrapper(expected_dims=(3, 2))
def clear_based_on_2d_mask(stack, mask):
    return _3d_stack_2d_filter(stack, np.logical_not(mask))


