from matplotlib import pyplot as plt
import numpy as np
import core_functions as cf
from collections import defaultdict
import os
from debugger import CustomDebugger

debugger = CustomDebugger()
cf.debugger = debugger


def traverse_and_match(main_root,
                       matching_rule='c', matching_map=None):
    """
    Traverses the main_root directory, looking for all the '.tif/.TIF' files, performs name matching
    then iterates through the resulting matched dictironary.

    Matching assumption is that except for the matching keys, the names are identical

    :param main_root: folder from which will be traversed in depth
    :param matching_rule: name modification to type mapping. Currently '' for no matching, 'color' for colors
    :param matching_map: {'pattern in the file name': color channel number}
    :return:
    """
    matched_images = defaultdict(lambda: ['']*len(matching_map.keys()))

    if matching_rule:
        assert(matching_map is not None)

    for current_location, sub_directories, files in os.walk(main_root):
        if files:
            for img in files:
                if ('.TIF' in img or '.tif' in img) and '_thumb_' not in img:
                    prefix = cf.split_and_trim(current_location, main_root)

                    img_codename = img.split(' ')[0].split('_')
                    color = matching_map[img_codename[-1]]
                    name_pattern = ' - '.join(prefix + img_codename[:-1])
                    matched_images[name_pattern][color] = os.path.join(current_location, img)

    delset = []

    for name_pattern, (color_set) in matched_images.iteritems():
        # debugger.logger.debug(color_set)
        if any([color == '' for color in color_set]):
            debugger.logger.info('in %s, colorset is broken:', name_pattern)
            for color_name in color_set:
                debugger.logger.info('\t %s', color_name)
            debugger.logger.info('name_pattern will be deleted')
            delset.append(name_pattern)

    for name_pattern in delset:
        del matched_images[name_pattern]

    for name_pattern, color_set in matched_images.iteritems():
        channels = []
        for color in color_set:
            channels.append(cf.tiff_stack_2_np_arr(color))

        yield name_pattern, channels


def stack_splitter(stack_group_generator):
    """
    Used when the tiff stack encodes a different information from the z-stack. assumes a single
    stack in the input generator

    :param stack_group_generator:
    :return:
    """
    for name_pattern, stack in stack_group_generator:
        channels = np.split(stack, stack.shape)
        yield name_pattern, channels


def name_channels(stack_group_generator, channel_names):
    """
    Assigns names to the channel for the future processing and bundles them together

    :param stack_group_generator:
    :param channel_names:
    :return:
    """

    for name_pattern, channels in stack_group_generator:
        group_dict = {'name pattern': name_pattern}
        group_dict['channel list'] = channel_names
        for chan_name, chan in zip(channel_names, channels):
            group_dict[chan_name] = chan

        yield group_dict

