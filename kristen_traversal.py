import numpy as np
import core_functions as cf
from collections import defaultdict
import logging
import os



def Kristen_traverse(main_root,
                    matching_rule='c', matching_map=None):

    matched_images = defaultdict(lambda: ['']*len(matching_map.keys()))
    tags_dict = defaultdict(lambda: [])
#     using a dictionary instead of list since grouped by mutliple colors and may need to superimpose them
#     provide matching map based on translator (a dictionary with color codes and associated colors/dyes)

    for current_location, sub_directories, files in os.walk(main_root):
        if files:
            for img in files:
                if ('.TIF' in img or '.tif' in img) and '_thumb_' not in img:
                    # img_codename = [img.split('.')[0]]   #linhao
                    prefix = cf.split_and_trim(current_location, main_root)

                    img_codename = [img.split('.')[0]]
                    name_pattern = ' - '.join(prefix + img_codename)
                    print "prefix", prefix
                    print "img codename", img_codename
                    print "name pattern", name_pattern
                    group_by = img_codename[0][:2]
                    # grouping by color
                    print "group by", group_by, type(group_by)
                    print
                    print
                    color = group_by
                    # matched_images[name_pattern][color] = os.path.join(current_location, img)

    for name_pattern, group_by, image_location in matched_images:
        print "reached for loop"
        print name_pattern
        stack = cf.tiff_stack_2_np_arr(image_location)
        stack = np.rollaxis(stack[0, :, :, :], 2)  # on the data where channels have not been split into z stacks
        channels = np.split(stack, stack.shape[0])
        channels = [chan[0, :, :] for chan in channels]
#     dont forget to add the csv/user input portion after this part is running smoothly

translator = {'C1':3, 'C3':0, 'C4':1}
# C1 = DAPI, C2 = GFP, C3 = mCherry
Kristen_traverse("/run/user/1000/gvfs/smb-share:server=10.17.0.219,share=common/Users/kristen/Split GFP quant_Andrei/", matching_map=translator)
