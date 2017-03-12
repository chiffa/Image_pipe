import numpy as np
import core_functions as cf
from collections import defaultdict
import logging
import os



def Kristen_traverse(main_root,
                    matching_rule='c', matching_map=None):

    matched_images = defaultdict(lambda: ['']*len(matching_map.keys()))
    tags_dict = defaultdict(lambda: [])
#     provide matching map based on translator (a dictionary with color codes and associated colors/dyes)

    for current_location, sub_directories, files in os.walk(main_root):
        if files:
            for img in files:
                if ('.TIF' in img or '.tif' in img) and '_thumb_' not in img:
                    # prefix = cf.split_and_trim(current_location, main_root)
                    # img_codename = img.split('.')[0].split(' ')
                    # img_codename = [img.split('.')[0]]   #linhao
                    prefix = cf.split_and_trim(current_location, main_root)

                    img_codename = [img.split('.')[0]]
                    name_pattern = ' - '.join(prefix + img_codename)
                    print "prefix", prefix
                    print "img codename", img_codename
                    # name_pattern = ' - '.join([prefix[1]] + img_codename)
                    print "name pattern", name_pattern
                    group_by = img_codename[0][:2]
                    # grouping by color
                    print "group by", group_by, type(group_by)
                    print
                    print
                    color = group_by
                    print type(matched_images)
                    matched_images[name_pattern][color] = os.path.join(current_location, img)
                    # group_by = img_codename[0].split('rpe')[1].split('dapi')[0].strip()
                    # matched_images.append((name_pattern, group_by, os.path.join(current_location, img)))

    for name_pattern, group_by, image_location in matched_images:
        print "reached for loop"
        print name_pattern
        stack = cf.tiff_stack_2_np_arr(image_location)
        stack = np.rollaxis(stack[0, :, :, :], 2)  # on the data where channels have not been split into z stacks
        channels = np.split(stack, stack.shape[0])
        channels = [chan[0, :, :] for chan in channels]


translator = {'C1':'DAPI', 'C3':'GFP', 'C4':'mCherry'}
Kristen_traverse("/run/user/1000/gvfs/smb-share:server=10.17.0.219,share=common/Users/kristen/Split GFP quant_Andrei/", matching_map=translator)
