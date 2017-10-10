import logging
import os
import csv
from collections import defaultdict
import numpy as np
from imagepipe.raw_functions import split_and_trim
from imagepipe.tools.helpers import tiff_stack_2_np_arr

logger = logging.getLogger('Default Debug Logger')
logger.setLevel(logging.DEBUG)

fh = logging.FileHandler('debug_log.log')
fh.setLevel(logging.DEBUG)

# create console handler with a higher log level
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)

# create formatter and add it to the handlers
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
fh.setFormatter(formatter)
ch.setFormatter(formatter)

# add the handlers to the logger
logger.addHandler(fh)
logger.addHandler(ch)


def name_channels(stack_group_generator, channel_names):
    """
    Assigns names to the channel for the future processing and bundles them together

    :param stack_group_generator: generator returning image stack,
    :param channel_names:
    :return:
    """
    for name_pattern, group_ids, channels in stack_group_generator:
        print 'name pattern', name_pattern
        print 'group ids', group_ids
        print 'channels', [chan.shape for chan in channels]
        group_dict = {'name pattern': name_pattern,
                      'group id': group_ids}
        group_dict['channel list'] = channel_names
        for chan_name, chan in zip(channel_names, channels):
            group_dict[chan_name] = chan

        yield group_dict


def color_based_traversal(main_root):
    """
    Traverses the main root directory pulling the data from the images. different layers are assumed
    to be different colors, images are assumed to be 2D.

    :param main_root:
    :return:
    """

    matched_images = []

    for current_location, sub_directories, files in os.walk(main_root):
        if files:
            for img in files:
                if ('.TIF' in img or '.tif' in img) and '_thumb_' not in img:
                    prefix = split_and_trim(current_location, main_root)
                    img_codename = [img.split('.')[0]]
                    name_pattern = ' - '.join(prefix + img_codename)
                    group_by = img_codename[0].split('rpe')[1].split('dapi')[0].strip()
                    matched_images.append((name_pattern, group_by, os.path.join(current_location, img)))

    for name_pattern, group_by, image_location in matched_images:
        stack = tiff_stack_2_np_arr(image_location)
        stack = np.rollaxis(stack[0, :, :, :], 2)  # on the data where channels have not been split into z stacks
        channels = np.split(stack, stack.shape[0])
        channels = [chan[0, :, :] for chan in channels]
        yield name_pattern, group_by, channels


def z_stack_based_traversal(main_root, matching_rule='c', matching_map=None):
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
    tags_dict = defaultdict(lambda: [])
    if matching_rule:
        assert(matching_map is not None)

    for current_location, sub_directories, files in os.walk(main_root):

        if files:

            for img in files:
                if ('.TIF' in img or '.tif' in img) and '_thumb_' not in img and 'w' in img:

                    prefix = split_and_trim(current_location, main_root)

                    img_codename = img.split(' ')[0].split('_')
                    color = matching_map[img_codename[-1]]
                    name_pattern = ' - '.join(prefix + img_codename[:-1])
                    matched_images[name_pattern][color] = os.path.join(current_location, img)
                    time_stamp = prefix[-1]

                    if time_stamp == 'HS':
                        time = 0
                    elif 'HS' in time_stamp:
                        time = -30
                    else:
                        time = time_stamp[3:-3]  # int(time_stamp[3:-3])
                    tags_dict[name_pattern] = []
                    tags_dict[name_pattern].append(time)  # time in the times series
                    _date = prefix[0][:8]
                    tags_dict[name_pattern].append("%s-%s-%s" % (_date[:2], _date[2:4], _date[4:]))
                    tags_dict[name_pattern].append(prefix[-2])

    delset = []

    for name_pattern, (color_set) in matched_images.iteritems():
        # debugger.logger.debug(color_set)
        if any([color == '' for color in color_set]):
            logger.info('in %s, colorset is broken:', name_pattern)
            for color_name in color_set:
                logger.info('\t %s', color_name)
            logger.info('name_pattern will be deleted')
            delset.append(name_pattern)

    for name_pattern in delset:
        del matched_images[name_pattern]

    user_input_about_new_csv_file = raw_input("To continue, enter 2. To start the process from the beginning, enter 1")

    if user_input_about_new_csv_file == '1':
        print "Preparing a new CSV file"
        initial_open = open("matched_images.csv", 'wb')
        # this is the file we need to save unless user provides input saying we can override it
        writer = csv.writer(initial_open, delimiter='\t')
        for key in matched_images:
            writer.writerow([key] + matched_images[key] + [0])
        initial_open.close()

    else:
        print "Continuing where the process last left off"
        file_exists = os.path.isfile("matched_images.tmp")

        if file_exists:

            open_tmp = open('matched_images.tmp', 'r')
            read_preexisting_tmp = csv.reader(open_tmp, delimiter = '\t')
            tmp_list = []
            for row in read_preexisting_tmp:
                tmp_list.append(row)
            open_tmp.close()

            open_csv = open('matched_images.csv', 'r')
            read_preexisting_csv = csv.reader(open_csv, delimiter = '\t')
            csv_list = []
            for row in read_preexisting_csv:
                csv_list.append(row)
            open_csv.close()

            for csv_row in csv_list:
                for tmp_row in tmp_list:
                    if csv_row[0] == tmp_row[0]:
                        csv_row[3] = tmp_row[3]

            open_csv_write = open('matched_images.csv', 'wb')
            override_csv = csv.writer(open_csv_write, delimiter='\t')
            for new_csv_row in csv_list:
                override_csv.writerow(new_csv_row)

    open_updated_csv_to_read = open('matched_images.csv', 'rb')
    csv_reader = csv.reader(open_updated_csv_to_read, delimiter='\t')

    open_tmp_to_write = open("matched_images.tmp", 'wb')
    writer_check_tmp = csv.writer(open_tmp_to_write, delimiter='\t')

    for row in csv_reader:
        name_pattern = row[0]
        color_set = [row[1], row[2]]
        if row[3] == 1:
            writer_check_tmp.writerow(row)
            continue
        channels = []
        for color in color_set:
            channels.append(tiff_stack_2_np_arr(color))
        print "name pattern", name_pattern
        print 'channels', channels
        yield name_pattern, tags_dict[name_pattern], channels
        row[3] = 1
        writer_check_tmp.writerow(row)