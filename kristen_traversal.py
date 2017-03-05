import numpy as np
import core_functions as cf
from collections import defaultdict
import logging
import os


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

def Kristen_traverse(main_root,
                    matching_rule='c', matching_map=None):

    # main root = /run/user/1000/gvfs/smb-share:server=10.17.0.219,share=common/Users/kristen/Split GFP quant_Andrei/20170209

    matched_images = defaultdict(lambda: [''] * len(matching_map.keys()))
    tags_dict = defaultdict(lambda: [])

    if matching_rule:
        assert (matching_map is not None)

    for current_location, sub_directories, files in os.walk(main_root):
        if files:
            for img in files:
                # if ('.TIF' in img or '.tif' in img) and '_thumb_' not in img:
                if ('.TIF' in img or '.tif' in img) and '_thumb_' not in img:
                    prefix = cf.split_and_trim(current_location, main_root)
                    img_codename = img.split(' ')[0].split('_')
                    color = matching_map[img_codename[-1]]
                    name_pattern = ' - '.join(prefix + img_codename[:-1])
                    matched_images[name_pattern][color] = os.path.join(current_location, img)
                    time_stamp = prefix[-1]

                    # if time_stamp == 'HS':
                    #     time = 0
                    # elif 'HS' in time_stamp:
                    #     time = -30
                    # else:
                    #     time = time_stamp[3:-3]  # int(time_stamp[3:-3])
                    # tags_dict[name_pattern] = []
                    # tags_dict[name_pattern].append(time)  # i dont think we need this commented out section above?
                    _date = prefix[0][3:11]
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

    for name_pattern, color_set in matched_images.iteritems():
        channels = []
        for color in color_set:
            channels.append(cf.tiff_stack_2_np_arr(color))

        print 'starting to analyze', name_pattern

        yield name_pattern, tags_dict[name_pattern], channels
