import recipy
import os
import traceback
import numpy as np
from PIL import Image
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





# Module-level user-changeable parameters parameters:
scaling_factor = (1.0, 1.0, 3.5)
mcc_cutoff = 0.05

# 0 is the protein of interest marker
# 1 is the organelle of interest marker
translator = {'w1488': 0,
              'w2561': 1,
              'C1': 1,
              'C2': 0}

# needs to be modified to properly group classes
classes = ['ry233-1', 'ry233-2', 'ry130-1', 'ry130-2']

# translation of file names to time-stamps relative to the heatshock moment
time_stamps = {
    'beforehs': -30,
    'hs30min': 0,
    'rec15min': 15,
    'rec30min': 30,
    'rec45min': 45,
    'rec60min': 60,
    'rec75min': 75,
    'rec90min': 90
}

time_stamp_coll = time_stamps.keys()

header = ['name pattern', 'GFP', 'mito marker', 'cross',
          'MCC mito in GFP %', 'MCC GFP in mito %',
          'AQVI GFP', 'AQVI mito', 'ill', 'cell_no',
          'mean width', 'mean length', 'cells with intact mitochondria %',
          'area of intact mitochondria %']


# Support functions and logic:
dtype2bits = {'uint8': 8,
              'uint16': 16,
              'uint32': 32}


def safe_dir_create(path):
    if not os.path.isdir(path):
        os.makedirs(path)


safe_dir_create('verification_bank')


def rm_nans(np_array):
    fltr = np.logical_not(np.isnan(np_array))
    return np_array[fltr]


# DEBUG frame with rendering:
class DebugFrame(object):

    def __init__(self):
        self.name_pattern = None
        self.gfp_collector = None
        self.gfp_clustering_markers = None
        self.labels = None
        self.segmented_cells_labels = None

        self.average_gfp_in_cell_mask = None
        self.cells_average_gfp_list = None
        self.non_dying_predicted = None
        self.std_err = None
        self.non_dying_cells_mask = None
        self.qualifying_gfp_mask = None

        self.gfp_collector_post = None
        self.mch_collector_post = None
        self.gfp_collector_pre = None
        self.mch_collector_pre = None

        self.mch_collector = None
        self.skeleton = None
        self.mean_width = None
        self.mean_length = None
        self.mito_binary_mask = None
        self.segmented_cells = None
        self.numbered_lables = None
        self.classification_pad = None
        self.paint_area = None
        self.numbered_skeleton_label = None

    def render_1(self):
        plt.figure(figsize=(20.0, 15.0))
        plt.title(self.name_pattern)

        plt.subplot(241)
        plt.imshow(self.gfp_collector, interpolation='nearest')

        plt.subplot(242)
        plt.imshow(self.gfp_clustering_markers, cmap='hot', interpolation='nearest')

        plt.subplot(243)
        plt.imshow(self.labels, cmap='gray', interpolation='nearest')

        plt.subplot(244)
        plt.imshow(self.segmented_cells_labels, cmap=plt.cm.spectral, interpolation='nearest')

        if self.std_err is not None:

            plt.subplot(245)
            plt.imshow(self.average_gfp_in_cell_mask, cmap='hot', interpolation='nearest')
            plt.colorbar()

            plt.subplot(246)
            plt.plot(self.cells_average_gfp_list, 'ko')
            plt.plot(self.non_dying_predicted, 'r')
            plt.plot(self.non_dying_predicted + self.std_err, 'g')
            plt.plot(self.non_dying_predicted - self.std_err, 'g')

            plt.subplot(247)
            plt.imshow(self.non_dying_cells_mask, cmap='gray', interpolation='nearest')

            plt.subplot(248)
            plt.imshow(self.qualifying_gfp_mask)

        plt.savefig('verification_bank/%s.png' % self.name_pattern)
        # plt.show()
        plt.close()

    def render_2(self):
        plt.figure(figsize=(20.0, 15.0))

        plt.subplot(221)
        plt.imshow(self.gfp_collector_pre, cmap='Greens')

        plt.subplot(222)
        plt.imshow(self.mch_collector_pre, cmap='Reds')

        plt.subplot(223)
        plt.imshow(self.gfp_collector_post, cmap='Greens')

        plt.subplot(224)
        plt.imshow(self.mch_collector_post, cmap='Reds')

        plt.savefig('verification_bank/core-%s.png' % self.name_pattern)
        # plt.show()
        plt.close()

    def render_3(self):
        plt.figure(figsize=(20.0, 15.0))

        ax1 = plt.subplot(231)
        plt.title(self.name_pattern)
        plt.imshow(self.mch_collector, cmap='Reds')
        plt.colorbar()
        plt.contour(self.labels, [0.5], colors='k')

        plt.subplot(232, sharex=ax1, sharey=ax1)
        plt.title('width ; length - av: %.2f ; %.2f' % (self.mean_width, self.mean_length))
        plt.imshow(self.skeleton, cmap=plt.cm.spectral, interpolation='nearest')
        plt.colorbar()
        plt.contour(self.mito_binary_mask, [0.5], colors='w')

        plt.subplot(233, sharex=ax1, sharey=ax1)
        plt.imshow(self.segmented_cells, cmap=plt.cm.spectral, interpolation='nearest')
        plt.colorbar()
        plt.contour(self.mito_binary_mask, [0.5], colors='w')

        plt.subplot(234, sharex=ax1, sharey=ax1)
        # numbered_skeleton_label ?
        plt.imshow(self.numbered_lables, cmap=plt.cm.spectral, interpolation='nearest')
        plt.colorbar()
        plt.contour(self.mito_binary_mask, [0.5], colors='w')

        plt.subplot(235, sharex=ax1, sharey=ax1)
        plt.imshow(self.paint_area, cmap='hot', interpolation='nearest')
        plt.colorbar()
        plt.contour(self.mito_binary_mask, [0.5], colors='w')

        plt.subplot(236, sharex=ax1, sharey=ax1)
        plt.imshow(self.classification_pad, cmap=plt.cm.spectral, interpolation='nearest')
        plt.colorbar()
        plt.contour(self.mito_binary_mask, [0.5], colors='w')

        plt.savefig('verification_bank/mitochondria-%s.png' % self.name_pattern)
        # plt.show()
        plt.close()

    def dump_debug_frames(self):
        self.render_2()
        self.render_1()
        self.render_3()


running_debug_frame = DebugFrame()

# debug of failings list:
sucker_list = []


def tiff_stack_2_np_arr(tiff_stack):
    """
    Loads the image from the tiff stack to a 3D numpy array

    :param tiff_stack:
    :return:
    """
    stack = [np.array(tiff_stack)]
    try:
        while 1:
            tiff_stack.seek(tiff_stack.tell() + 1)
            stack.append(np.array(tiff_stack))
    except EOFError:
        pass

    return np.array(stack)


def gamma_stabilize_and_smooth(tiff_stack,
                               alpha_clean=5, smoothing_px=1.5, debug=False):
    """
    Performs the initial conversion and de-noising of the tiff stack

    :param tiff_stack:
    :param alpha_clean:
    :param smoothing_px:
    :param debug:
    :return:
    """
    current_image = tiff_stack_2_np_arr(tiff_stack)
    bits = dtype2bits[current_image.dtype.name]

    if debug:
        print np.max(current_image), np.min(current_image), np.median(current_image)
        plt.hist(current_image.flatten(), 100)
        plt.show()

    stabilized = (current_image - np.min(current_image))/(float(2**bits) - np.min(current_image))
    stabilized[stabilized < alpha_clean*np.median(stabilized)] = 0

    if debug:
        print np.max(current_image), np.min(current_image), np.median(current_image)
        plt.hist(current_image.flatten(), 100)
        plt.show()

    if smoothing_px:
        for i in range(0, stabilized.shape[0]):
            stabilized[i, :, :] = gaussian_filter(stabilized[i, :, :],
                                                  smoothing_px, mode='constant')
            stabilized[stabilized < 5*np.mean(stabilized)] = 0

    if debug:
        print np.max(current_image), np.min(current_image), np.median(current_image)
        plt.hist(current_image.flatten(), 100)
        plt.show()

        for i in range(0, stabilized.shape[0]):
            plt.imshow(stabilized[i, :, :] > mcc_cutoff, cmap='gray', vmin=0., vmax=1.)
            plt.show()

    return stabilized


# can be edited to a debug wrapper
def detect_ill_cells(cell_labels, gfp_collector):
    """
    Logic that determines outliers that look like dead cells in the gfp channel projection.
    Requires at least 5 non-dead cells in the image.

    :param cell_labels:
    :param gfp_collector:
    :return:
    """
    cells_average_gfp_list = extract_average_qualifying_gfp_for_cell_regions(cell_labels,
                                                                             gfp_collector)
    non_dying_cells = detect_upper_outliers(cells_average_gfp_list)
    non_dying_cells_mask = paint_mask(cell_labels, non_dying_cells)

    running_debug_frame.non_dying_cells_mask = non_dying_cells_mask

    return non_dying_cells_mask


def extract_average_qualifying_gfp_for_cell_regions(cell_labels, gfp_collector):

    cells_average_gfp_list = []
    average_gfp_in_cell_mask = np.zeros_like(cell_labels).astype(np.float64)
    qualifying_gfp_mask = gfp_collector > np.median(gfp_collector[gfp_collector > 0])

    for i in range(1, np.max(cell_labels) + 1):

        current_mask = cell_labels == i
        current_cell_gfp = gfp_collector[np.logical_and(current_mask, qualifying_gfp_mask)]

        if len(current_cell_gfp) == 0:
            continue

        gfp_percentile = np.percentile(current_cell_gfp, 50)
        gfp_average = np.average(gfp_collector[np.logical_and(current_mask,
                                                              gfp_collector > gfp_percentile)])
        cells_average_gfp_list.append(gfp_average)
        average_gfp_in_cell_mask[current_mask] = gfp_average

    running_debug_frame.average_gfp_in_cell_mask = average_gfp_in_cell_mask
    running_debug_frame.qualifying_gfp_mask = qualifying_gfp_mask

    return cells_average_gfp_list


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
    non_dying_predicted = intercept + slope * np.array(cell_no)
    non_dying_cells = arg_sort[np.array(cell_no)[np.array(non_dying_predicted + std_err) <
                                                 np.array(cells_average_gfp_list)]]

    running_debug_frame.cells_average_gfp_list = cells_average_gfp_list
    running_debug_frame.non_dying_predicted = non_dying_predicted
    running_debug_frame.std_err = std_err

    return non_dying_cells


def paint_mask(cell_labels, non_dying_cells):
    non_dying_cells_mask = np.zeros_like(cell_labels).astype(np.uint8)
    if non_dying_cells.tolist() != []:
        # print 'enter non_dying_cells correction'
        for idx in non_dying_cells.tolist():
            non_dying_cells_mask[
                cell_labels == idx + 1] = 1  # indexing starts from 1, not 0 for the labels
            # print 'updated %s' % (idx + 1)
    return non_dying_cells_mask


def segment_out_cells(base):
    # TODO: try using OTSU for GFP thresholding

    sel_elem = disk(2)
    gfp_collector = np.sum(base, axis=0)
    gfp_clustering_markers = np.zeros(gfp_collector.shape, dtype=np.uint8)
    # random walker segment
    gfp_clustering_markers[gfp_collector > np.mean(gfp_collector) * 2] = 2
    gfp_clustering_markers[gfp_collector < np.mean(gfp_collector) * 0.20] = 1
    labels = random_walker(gfp_collector, gfp_clustering_markers, beta=10, mode='bf')
    # round up the labels and set the background to 0 from 1
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

    # log debugging data
    running_debug_frame.gfp_collector = gfp_collector
    running_debug_frame.gfp_clustering_markers = gfp_clustering_markers
    running_debug_frame.labels = labels
    running_debug_frame.segmented_cells_labels = segmented_cells_labels

    return gfp_collector, segmented_cells_labels


def skeletonize_mitochondria(mch_channel):
    mch_collector = np.max(mch_channel, axis=0)  # TODO: check max projection v.s. sum
    skeleton_labels = np.zeros(mch_collector.shape, dtype=np.uint8)

    # thresh = np.max(mch_collector)/2.
    thresh = threshold_otsu(mch_collector)
    # use adaptative threshold? => otsu seems to be sufficient in this case

    skeleton_labels[mch_collector > thresh] = 1
    skeleton2 = skeletonize(skeleton_labels)
    skeleton, distance = medial_axis(skeleton_labels, return_distance=True)
    active_threshold = np.mean(mch_collector[skeleton_labels]) * 5

    # print active_threshold
    transform_filter = np.zeros(mch_collector.shape, dtype=np.uint8)
    transform_filter[np.logical_and(skeleton > 0, mch_collector > active_threshold)] = 1
    skeleton = transform_filter * distance

    skeleton_ma = np.ma.masked_array(skeleton, skeleton > 0)
    skeleton_convolve = ndi.convolve(skeleton_ma, np.ones((3, 3)), mode='constant', cval=0.0)
    divider_convolve = ndi.convolve(transform_filter, np.ones((3, 3)), mode='constant', cval=0.0)
    skeleton_convolve[divider_convolve > 0] = skeleton_convolve[divider_convolve > 0] / \
                                              divider_convolve[divider_convolve > 0]
    new_skeleton = np.zeros_like(skeleton)
    new_skeleton[skeleton2] = skeleton_convolve[skeleton2]
    skeleton = new_skeleton

    return skeleton_labels, mch_collector, skeleton, transform_filter


def measure_skeleton_stats(numbered_labels, skeleton, transform_filter):

    numbered_skeleton, object_no = ndi.label(transform_filter, structure=np.ones((3, 3)))

    # print numbered_skeleton.shape, np.min(numbered_skeleton), np.max(numbered_skeleton)
    collector = []
    paint_area = np.zeros_like(numbered_labels)
    paint_length = np.zeros_like(numbered_labels)

    for contig_no in range(1, object_no + 1):
        vals = skeleton[numbered_skeleton == contig_no]
        current_label = np.max(numbered_labels[numbered_skeleton == contig_no])
        area, support = (np.sqrt(np.sum((numbered_labels == current_label).astype(np.int8))),
                         len(vals))

        if area < 3:
            skeleton[numbered_skeleton == contig_no] = 0
            transform_filter[numbered_skeleton == contig_no] = 0

        else:
            paint_area[numbered_labels == current_label] = area
            paint_length[numbered_labels == current_label] = support
            collector.append([area, support])

    collector = np.array(collector)
    running_debug_frame.numbered_skeleton_label = numbered_labels

    return collector, paint_length, paint_area


def compute_mito_fragmentation(name_pattern, skeleton_labels, mch_collector, skeleton,
                               transform_filter, segmented_cells):

    numbered_lables, lables_no = ndi.label(skeleton_labels, structure=np.ones((3, 3)))
    collector, paint_length, paint_area = measure_skeleton_stats(numbered_lables,
                                                                  skeleton, transform_filter)
    classification_pad = np.zeros_like(segmented_cells)
    classification_roll = []

    for i in range(1, np.max(segmented_cells)+1):
        pre_mask = segmented_cells == i
        current_mask = np.logical_and(pre_mask, skeleton_labels > 0)
        if len(paint_length[current_mask]) == 0:
            classification_roll.append(-1)
            classification_pad[pre_mask] = -1
        else:
            length = np.mean(np.unique(paint_length[current_mask]))
            area = np.mean(np.unique(paint_area[current_mask]))
            if length < 20 or area < 5:
                classification_pad[pre_mask] = 1
                classification_roll.append(1)
            else:
                classification_pad[pre_mask] = 2
                classification_roll.append(2)

    intact = np.logical_and(paint_length > 20, paint_area > 5)
    broken = np.logical_and(np.logical_or(paint_length < 20, paint_area < 5), paint_area > 1)

    if np.any(intact) or np.any(intact):
        mito_summary = float(np.sum(intact.astype(np.int8))) / \
                       float(np.sum(intact.astype(np.int8)) + np.sum(broken.astype(np.int8)))
    else:
        mito_summary = np.nan

    if len(collector) == 0:
        mean_width, mean_length = [np.NaN, np.NaN]
    else:
        mean_width, mean_length = np.mean(collector, axis=0).tolist()

    running_debug_frame.mch_collector = mch_collector
    running_debug_frame.skeleton = skeleton
    running_debug_frame.mean_width = mean_width
    running_debug_frame.mean_length = mean_length
    running_debug_frame.mito_binary_mask = skeleton_labels
    running_debug_frame.segmented_cells = segmented_cells
    running_debug_frame.numbered_lables = numbered_lables
    running_debug_frame.classification_pad = classification_pad
    running_debug_frame.paint_area = paint_area

    classification_array = np.array(classification_roll)
    classification_array = classification_array[classification_array > 0] - 1

    return [mean_width, mean_length, np.mean(classification_array), mito_summary]


def _3d_stack_2d_filter(_3d_stack, _2d_filter):
    new_stack = np.zeros_like(_3d_stack)
    new_stack[:, _2d_filter] = _3d_stack[:, _2d_filter]
    return new_stack


def _2d_stack_2d_filter(_2d_stack, _2d_filter):
    new_stack = np.zeros_like(_2d_stack)
    new_stack[_2d_filter] = _2d_stack[_2d_filter]
    return new_stack


def extract_statistics(_labels, _marked_prot, _mch_collector, _organelle_marker, _skeleton,
                       _transform_filter, name_pattern, segmented_cells_labels,
                       cell_no='NaN', ill='NaN'):

    seg4 = compute_mito_fragmentation(name_pattern, _labels, _mch_collector,
                                      _skeleton, _transform_filter, segmented_cells_labels)

    seg0 = [name_pattern]
    seg1 = [np.sum(_marked_prot * _marked_prot),
            np.sum(_organelle_marker * _organelle_marker),
            np.sum(_marked_prot * _organelle_marker)]

    if np.sum(_marked_prot) == 0.0:
        pre_seg2_2 = np.nan
    else:
        pre_seg2_2 = np.sum(_marked_prot[_organelle_marker > mcc_cutoff]) / np.sum(_marked_prot)

    if np.sum(_organelle_marker) == 0.0:
        pre_seg2_1 = np.nan
    else:
        pre_seg2_1 = np.sum(_organelle_marker[_marked_prot > mcc_cutoff]) / np.sum(_organelle_marker)

    seg2 = [pre_seg2_1, pre_seg2_2]

    seg3 = [np.median(_marked_prot[_organelle_marker > mcc_cutoff]),
            np.median(_organelle_marker[_organelle_marker > mcc_cutoff]),
            ill,
            cell_no]

    return [seg0 + seg1 + seg2 + seg3 + seg4]


def analyze_gfp_mch_paired_stacks(name_pattern, gfp_marked_prot, mch_organelle_marker,
                                  segment_out_ill=True, per_cell=False):
    """
    Stitches the analysis pipeline together

    :param name_pattern:
    :param gfp_marked_prot:
    :param mch_organelle_marker:
    :param segment_out_ill:
    :param per_cell:
    :return:
    """
    running_debug_frame.name_pattern = name_pattern
    gfp_collector, segmented_cells_labels = segment_out_cells(gfp_marked_prot)

    running_debug_frame.gfp_collector_pre = np.sum(gfp_marked_prot, axis=0)
    running_debug_frame.mch_collector_pre = np.sum(mch_organelle_marker, axis=0)

    if segment_out_ill:
        non_dying_cells_mask = detect_ill_cells(segmented_cells_labels, gfp_collector)
        gfp_marked_prot = _3d_stack_2d_filter(gfp_marked_prot,
                                              np.logical_not(non_dying_cells_mask))
        mch_organelle_marker = _3d_stack_2d_filter(mch_organelle_marker,
                                                   np.logical_not(non_dying_cells_mask))

    running_debug_frame.gfp_collector_post = np.sum(gfp_marked_prot, axis=0)
    running_debug_frame.mch_collector_post = np.sum(mch_organelle_marker, axis=0)

    skeleton_labels, mch_collector, skeleton, transform_filter = \
        skeletonize_mitochondria(mch_organelle_marker)

    # required to properly load the debug state.
    compute_mito_fragmentation(name_pattern, skeleton_labels, mch_collector, skeleton,
                               transform_filter, segmented_cells_labels)

    running_debug_frame.dump_debug_frames()

    if per_cell:

        seg_stack = []

        for cell_no in range(1, np.max(segmented_cells_labels)+1):
            current_mask = segmented_cells_labels == cell_no
            current_mask = current_mask[:, :]

            if segment_out_ill:
                ill = np.median(non_dying_cells_mask[current_mask])
            else:
                ill = 'NaN'

            _organelle_marker = _3d_stack_2d_filter(mch_organelle_marker, current_mask)
            _marked_prot = _3d_stack_2d_filter(gfp_marked_prot, current_mask)
            _skeleton_labels = _2d_stack_2d_filter(skeleton_labels, current_mask)
            _mch_collector = _2d_stack_2d_filter(mch_collector, current_mask)
            _skeleton = _2d_stack_2d_filter(skeleton, current_mask)
            _transform_filter = _2d_stack_2d_filter(transform_filter, current_mask)

            seg_stack += extract_statistics(_skeleton_labels, _marked_prot, _mch_collector,
                                            _organelle_marker, _skeleton, _transform_filter,
                                            name_pattern, segmented_cells_labels, cell_no, ill)
        return seg_stack

    else:

        return extract_statistics(skeleton_labels, gfp_marked_prot, mch_collector, mch_organelle_marker,
                                  skeleton, transform_filter, name_pattern, segmented_cells_labels)


def folder_structure_traversal(main_root, per_cell=False, mammalian=False):
    replicas = defaultdict(lambda: [0, 0])
    results_collector = []

    for current_location, sub_directories, files in os.walk(main_root):
        print current_location
        print '\t', files

        if files:
            pair_channels_in_folder(current_location, files, replicas, mammalian=mammalian)
            results_collector += analyze_matched_stack(per_cell, replicas, results_collector)
            replicas = defaultdict(lambda: [0, 0])

    # table_post_processing(results_collector)

    name_mod = 'yeast'
    per_cell_name = 'full_image'

    if mammalian:
        name_mod = 'mammalian'

    if per_cell:
        per_cell_name = 'per_cell'

    dump_name = 'results-nn-%s-%s.csv' % (per_cell_name, name_mod)

    dump_results_table(results_collector, dump_name)


def pair_channels_in_folder(current_location, files, replicas, mammalian=False):

    for img in files:
        if ('.TIF' in img or '.tif' in img) and '_thumb_' not in img:
            prefix = current_location.split('\\')[4:]

            if mammalian:
                img_codename = img.split('-')
                color = translator[img_codename[0]]
                name_pattern = ' - '.join(prefix + img_codename[1:])

            else:
                img_codename = img.split(' ')[0].split('_')
                color = translator[img_codename[-1]]
                name_pattern = ' - '.join(prefix + img_codename[:-1])

            current_image = Image.open(os.path.join(current_location, img))
            print '%s image was parsed, code: %s %s' % (img, name_pattern, color)
            replicas[name_pattern][color] = gamma_stabilize_and_smooth(current_image)


def analyze_matched_stack(per_cell, replicas, results_collector):
    results_subcollector = []
    for name_pattern, (w1448, w2561) in replicas.iteritems():
        print name_pattern
        try:
            results_subcollector += analyze_gfp_mch_paired_stacks(name_pattern, w1448, w2561,
                                                               segment_out_ill=True,
                                                               per_cell=per_cell)
        except Exception as my_exception:
            print traceback.print_exc(my_exception)
            sucker_list.append(name_pattern)
    return results_subcollector


def dump_results_table(results_collector, fname):
    with open(fname, 'wb') as output:
        csv_writer = writer(output, )
        csv_writer.writerow(header)
        for item in results_collector:
            csv_writer.writerow(item)


def table_post_processing(results_collector):
    stats_collector = []

    print '\n Summary of the analysis:'

    for line in results_collector:
        class_name, _time_stamps = classify(line[0])
        stats_collector.append([class_name, _time_stamps, line[6], line[7]])

    stats_collector = np.array(stats_collector)

    for class_name in classes:
        class_set_filter = stats_collector[:, 0] == class_name
        if any(class_set_filter):
            class_set = stats_collector[class_set_filter, :]
            print class_name
            final_stats_collector_x = []
            final_stats_collector_y = []
            final_stats_collector_e = []
            raw_times = {}
            for time_stamp in time_stamp_coll:
                time_stamp_filter = class_set[:, 1] == time_stamp
                if any(time_stamp_filter):
                    time_stamp_set = class_set[time_stamp_filter, :]
                    mean = np.nanmean(time_stamp_set[:, 2].astype(np.float64))
                    err = np.nanstd(time_stamp_set[:, 2].astype(np.float64)) / \
                          np.sqrt(len(time_stamp_set[:, 2]))*1.96
                    raw_times[time_stamp] = rm_nans(time_stamp_set[:, 2].astype(np.float64))
                    print '\t time: %s, mean: %s, err: %s' % (time_stamp, mean, err)
                    final_stats_collector_x.append(time_stamps[time_stamp])
                    final_stats_collector_y.append(mean)
                    final_stats_collector_e.append(err)

            time_translator = dict([(item, i) for i, item in enumerate(sorted(raw_times.keys()))])

            samples_n = len(raw_times.keys())
            print samples_n
            p_val_array = np.array((samples_n, samples_n))

            for time1, time2 in combinations(sorted(raw_times.keys()), 2):
                print time1, time2
                print time_translator[time1], time_translator[time2]
                print ttest_ind(raw_times[time1], raw_times[time2])
                _, p_val = ttest_ind(raw_times[time1], raw_times[time2])
                p_val_array[time_translator[time1], time_translator[time2]] = p_val

            print p_val_array

            plt.errorbar(final_stats_collector_x, final_stats_collector_y, final_stats_collector_e,
                         label=class_name)

    plt.legend()
    plt.show()


def classify(naming_code):
    naming_code = str.lower(naming_code)
    time_stamp = None
    class_name = 'unclassified'

    for _time_stamp in time_stamps.keys():
        if _time_stamp in naming_code:
            time_stamp = _time_stamp
            break

    for _class_name in classes:
        if _class_name in naming_code:
            class_name = _class_name

    return class_name, time_stamp


if __name__ == "__main__":

    # folder_structure_traversal("L:\\Users\linghao\\Data for quantification\\Yeast\\NEW data for analysis",
    #                            per_cell=True)

    "L:\\Users\\linghao\\Data for quantification\\Yeast\\NEW data for analysis"

    folder_structure_traversal("L:\\Users\\linghao\\Data for quantification\\Yeast\\NEW data for analysis",
                               per_cell=False)



