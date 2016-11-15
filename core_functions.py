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
from skimage.filters import threshold_otsu, rank
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


# TODO: convert to multiple input/output as well
# allow to override the expected dims check
# allow to perform expected output dims check => Problem, the mapping rules are not very clear.
# map in to out if only one output is passed to the function and function signature is n = n
# remove the logging channel as well

# can we pipe the in_dims into the out_dims if they are not explicitely provided?

@doublewrap
def generator_wrapper(f, in_dims=(3,), out_dims=None):

    if out_dims is None:
            out_dims = in_dims

    @wraps(f)
    def inner_wrapper(*args, **kwargs):
        """
        converts a function to accepting a generator of named dicts and adds channel selection logic
        """

        iterator = args[0]
        args = args[1:]

        if 'in_channel' in kwargs:
            #Start in/out channel logic

            in_chan = kwargs['in_channel']  # Multiple arguments
            del kwargs['in_channel']

            if not list_not_string(in_chan):  # convert string to list of strings
                in_chan = [in_chan]

            if 'out_channel' in kwargs:
                out_chan = kwargs['out_channel']  # output explicitely provided
                del kwargs['out_channel']
                if not list_not_string(out_chan):  # convert string to list of strings
                    out_chan = [out_chan]

            else:  # implicit output, bound to in_channel only if a single input is provided
                if len(in_chan) == 1:
                    print 'Input %s will be overwritten by function %s' % (in_chan[0], f.__name__)
                    out_chan = in_chan
                else:
                    raise PipeArgError('Please provide out_channel argument')

            if len(in_chan) != len(in_dims):
                print in_chan, in_dims
                print len(in_chan), len(in_dims)
                raise PipeArgError('More inbound channels are piped than function allows')

            if len(out_chan) != len(out_dims):
                print out_chan, out_dims
                print len(out_chan), len(out_dims)
                raise PipeArgError('More outbound channels are piped than function allows')

            # end in/out channel logic

            for name_space in iterator:
                # start args prepare
                args_puck = []

                for i, chan in enumerate(in_chan):
                    if in_dims[i] and len(name_space[chan].shape) != in_dims[i]:
                        print f.__name__
                        print chan, len(name_space[chan].shape), in_dims[i]
                        raise PipeArgError('Mismatched inbound channel dimension for channel. %s is of dim %s, expected %s'%
                                           (chan, len(name_space[chan].shape), in_dims[i]))

                    args_puck.append(name_space[chan])

                local_args = tuple(args_puck) + args
                # end args prepare

                return_puck = f(*local_args, **kwargs)

                if return_puck is None and out_chan[0] == '_':
                    yield name_space

                # start output prepare
                if not isinstance(return_puck, tuple):
                    return_puck = (return_puck, )

                for i, chan in enumerate(out_chan):
                    if out_dims[i] and len(return_puck[i].shape) != out_dims[i]:
                        raise PipeArgError('Mismatched outgoing channel dimension for channel. %s is of dim %s, expected %s' %
                                           (chan, len(return_puck[i].shape), out_dims[i]))
                    name_space[chan] = return_puck[i]
                # end output prepare

                yield name_space

        else:
            for name_space in iterator:
                local_args = (name_space,) + args
                yield f(*local_args, **kwargs)

    return inner_wrapper


def splitter(outer_generator, to, sources, mask):
    """
    Creates a secondary namespace by using mask as a pad to conserve only certain segments in sources

    :param outer_generator:
    :param to:
    :param sources:
    :param mask:
    :return:
    """
    for primary_namespace in outer_generator:
        primary_namespace[to] = {}
        unique_vals = np.unique(primary_namespace[mask])
        unique_vals = unique_vals[unique_vals > 0]

        primary_namespace[to]['_pad'] = (unique_vals, mask)  # used to rebuild padded images

        for val in unique_vals:
            secondary_namespace = {}
            primary_namespace[to][val] = secondary_namespace

            for chan in sources:
                local_mask = primary_namespace[mask] == val
                if len(primary_namespace[chan].shape) == 2:
                    base_chan = _2d_stack_2d_filter(primary_namespace[chan], local_mask)
                elif len(primary_namespace[chan].shape) == 3:
                    base_chan = _3d_stack_2d_filter(primary_namespace[chan], local_mask)
                else:
                    raise PipeArgError('masking impossible: dims not match, base channel %s is of dim %s' %
                                       (chan, len(primary_namespace[chan].shape)))

                secondary_namespace[chan] = base_chan

        yield primary_namespace


def for_each(outer_generator, embedded_transformer, inside, **kwargs):

    for primary_namespace in outer_generator:
        embedded_transformer(primary_namespace[inside].itervalues(), **kwargs)
        yield primary_namespace



def paint_from_mask(outer_generator, based_on, in_anchor, out_channel=None):

    if out_channel is None:
        out_channel = in_anchor

    for primary_namespace in outer_generator:
        secondary_namespace = primary_namespace[based_on]
        mask = secondary_namespace['_pad'][1]
        mask_values = secondary_namespace['_pad'][0]
        accumulator = np.zeros_like(mask)

        for unique_value in mask_values:
            accumulator[mask == unique_value] = secondary_namespace[unique_value][in_anchor]

        primary_namespace[out_channel] = accumulator
        yield primary_namespace


def tile_from_mask(outer_generator, based_on, in_anchor, out_channel=None):

    if out_channel is None:
        out_channel = in_anchor

    for primary_namespace in outer_generator:
        secondary_namespace = primary_namespace[based_on]
        mask = secondary_namespace['_pad'][1]
        mask_values = secondary_namespace['_pad'][0]
        accumulator = np.zeros_like(mask)

        for unique_value in mask_values:
            accumulator[mask == unique_value] = secondary_namespace[unique_value][in_anchor][mask == unique_value]

        primary_namespace[out_channel] = accumulator
        yield primary_namespace


def summarize(outer_generator, based_on, in_anchor, out_channel=None):

    if out_channel is None:
        out_channel = in_anchor

    # summarizes the inner generator based_on based on the bound name.
    # actually, should be a function wrapper
    pass


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


@generator_wrapper(in_dims=(3,), out_dims=(2,))
def sum_projection(current_image):
    return np.sum(current_image, axis=0)


@generator_wrapper(in_dims=(2,))
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

    return segmented_cells_labels


@generator_wrapper(in_dims=(2,))
def simple_segment(base_channel, min_px_radius=3):
    numbered_skeleton, object_no = ndi.label(base_channel, structure=np.ones((3, 3)))

    for label in range(1, object_no+1):
        px_radius = np.sqrt(np.sum((numbered_skeleton == label).astype(np.int8)))
        if px_radius < min_px_radius:
            numbered_skeleton[numbered_skeleton == label] = 0

    return numbered_skeleton


@generator_wrapper(in_dims=(2,))
def qualifying_gfp(max_sum_projection):
    return max_sum_projection > np.median(max_sum_projection[max_sum_projection > 0])


@generator_wrapper(in_dims=(2, 2, 2), out_dims=(1, 2))
def aq_gfp_per_region(cell_labels, max_sum_projection, qualifying_gfp_mask):

    cells_average_gfp_list = []
    cells_average_gfp_pad = np.zeros_like(cell_labels)

    for i in range(1, np.max(cell_labels) + 1):

        current_mask = cell_labels == i
        current_cell_gfp = max_sum_projection[np.logical_and(current_mask, qualifying_gfp_mask)]

        if len(current_cell_gfp) == 0:
            continue

        gfp_percentile = np.percentile(current_cell_gfp, 50)
        gfp_average = np.average(max_sum_projection[np.logical_and(current_mask, max_sum_projection > gfp_percentile)])
        cells_average_gfp_list.append(gfp_average)
        cells_average_gfp_pad[current_mask] = gfp_average

    return np.array(cells_average_gfp_list), cells_average_gfp_pad


@generator_wrapper(in_dims=(1,), out_dims=(1, 1, None))
def detect_upper_outliers(cells_average_gfp_list):
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

    return upper_outliers, predicted_average_gfp, std_err


@generator_wrapper(in_dims=(2, 1), out_dims=(2,))
def paint_mask(label_masks, labels_to_paint):
    mask_to_paint = np.zeros_like(label_masks).astype(np.uint8)

    if labels_to_paint.tolist() != np.array([]):
        for idx in labels_to_paint.tolist():
            mask_to_paint[label_masks == idx + 1] = 1  # indexing starts from 1, not 0 for the labels

    return mask_to_paint


@generator_wrapper(in_dims=(3, 2), out_dims=(3,))
def clear_based_on_2d_mask(stack, mask):
    return _3d_stack_2d_filter(stack, np.logical_not(mask))


@generator_wrapper
def binarize_3d(float_volume, mcc_cutoff):
    binary_volume = np.zeros_like(float_volume)
    binary_volume[float_volume > mcc_cutoff] = 1
    return binary_volume.astype(np.bool)


@generator_wrapper(in_dims=(3, 3), out_dims=(None,))
def binary_inclusion_3d(float_volume, binary_volume):
    m_q_v_i = np.median(float_volume[binary_volume])
    return m_q_v_i


@generator_wrapper(in_dims=(2,), out_dims=(2,))
def binarize_2d(float_surface, cutoff_type='static', mcc_cutoff=None):
    if cutoff_type == 'otsu':
        mcc_cutoff = threshold_otsu(float_surface)

    elif cutoff_type == 'local otsu':
        selem = disk(5)
        mcc_cutoff = rank.otsu(float_surface, selem)

    elif cutoff_type == 'static':
        pass

    else:
        raise PipeArgError('unknown cutoff type')

    binary_stack = np.zeros_like(float_surface).astype(np.bool)
    binary_stack[float_surface > mcc_cutoff] = 1

    return binary_stack


@generator_wrapper(in_dims=(2, 2), out_dims=(2,))
def agreeing_skeletons(float_surface, mito_labels):
    topological_skeleton = agreeing_skeletons(mito_labels)
    medial_skeleton, distance = medial_axis(mito_labels, return_distance=True)

    # TODO: test without the active threshold surface
    active_threshold = np.mean(float_surface[mito_labels]) * 5  # todo: remove * 5?
    transform_filter = np.zeros(mito_labels.shape, dtype=np.uint8)
    transform_filter[np.logical_and(medial_skeleton > 0, float_surface > active_threshold)] = 1
    # transform filter is basically medial_skeleton on a field above threshold (mean*5 - wow, that's a lot)
    medial_skeleton = transform_filter * distance

    median_skeleton_masked = np.ma.masked_array(medial_skeleton, medial_skeleton > 0)
    skeleton_convolve = ndi.convolve(median_skeleton_masked, np.ones((3, 3)),
                                     mode='constant', cval=0.0)
    divider_convolve = ndi.convolve(transform_filter, np.ones((3, 3)),
                                    mode='constant', cval=0.0)
    skeleton_convolve[divider_convolve > 0] = skeleton_convolve[divider_convolve > 0] / \
                                              divider_convolve[divider_convolve > 0]

    skeletons = np.zeros_like(medial_skeleton)
    skeletons[topological_skeleton] = skeleton_convolve[topological_skeleton]

    return skeletons


@generator_wrapper(in_dims=(2, 2), out_dims=(None, 2, 2, 2))
def classify_fragmentation_for_mitochondria(label_mask, skeletons):
    # what if no mitochondria currently found?
    # what if we want to compare the surface of fragmented mitochondria v.s. non-fragmented ones?

    # well, one thing for sure, there is no way of escalating the skeleton/mito supression if they
    # are too small => we will need a filter on the label

    # maybe it actually is a good idea to get the mask manipulation for all areas in the skeleton

    mask_items = np.unique(label_mask)
    mask_items = mask_items[mask_items > 0].tolist()

    radius_mask = np.zeros_like(label_mask)
    support_mask = np.zeros_like(label_mask)
    classification_mask = np.zeros_like(label_mask)
    classification_roll = []
    weights = []

    for label in mask_items:
        px_radius = np.sqrt(np.sum((label_mask == label).astype(np.int8)))
        support = len(skeletons[label_mask == label])

        if px_radius < 5 or support < 20:
            classification = 1  # fragment of a broken mitochondria
        else:
            classification = 2  # mitochondria is intact

        radius_mask += (label_mask == label).astype(np.float)*px_radius
        support_mask += (label_mask == label).astype(np.float)*support
        classification_mask += (label_mask == label).astype(np.float)*classification

        classification_roll.append(classification)
        weights.append(radius_mask)

    classification_roll = np.array(classification_roll)
    weights = np.array(weights)
    final_classification = np.average(classification_roll, weights=weights)

    return final_classification, classification_mask, radius_mask, support_mask
