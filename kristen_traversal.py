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


    matched_images = []

    for current_location, sub_directories, files in os.walk(main_root):
        if files:
            print "Files: True"
            for img in files:
                if ('.TIF' in img or '.tif' in img) and '_thumb_' not in img:
                    prefix = cf.split_and_trim(current_location, main_root)

                    img_codename = [img.split('.')[0]]
                    name_pattern = ' - '.join(prefix + img_codename)
                    print "image codename is", img_codename
                    print "prefix is", prefix
                    # group_by = img_codename[0].split('rpe')[1].split('dapi')[0].strip()
                    group_by = pass
                    matched_images.append((name_pattern, group_by, os.path.join(current_location, img)))

    for name_pattern, group_by, image_location in matched_images:
        stack = cf.tiff_stack_2_np_arr(image_location)
        stack = np.rollaxis(stack[0, :, :, :], 2)  # on the data where channels have not been split into z stacks
        channels = np.split(stack, stack.shape[0])
        channels = [chan[0, :, :] for chan in channels]

        yield name_pattern, group_by, channels

Kristen_traverse("/run/user/1000/gvfs/smb-share:server=10.17.0.219,share=common/Users/kristen/Split GFP quant_Andrei/")