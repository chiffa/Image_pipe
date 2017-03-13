import numpy as np
import core_functions as cf
from collections import defaultdict
import logging as logger
import os
from pympler import  muppy, summary
import csv

translator = {'C1':3,
              'C3':0,
              'C4':1}


def Kristen_traverse(main_root,
                    matching_rule='c', matching_map=None):
    matched_images = defaultdict(lambda: [''] * len(matching_map.keys()))
    print matched_images
    tags_dict = defaultdict(lambda: [])

    if matching_rule:
        assert (matching_map is not None)

    for current_location, sub_directories, files in os.walk(main_root):
        print files

        if files:

            for img in files:
                if ('.TIF' in img or '.tif' in img) and '_thumb_' not in img and 'w' in img:

                    prefix = cf.split_and_trim(current_location, main_root)
                    print "prefix", prefix, type(prefix)

                    img_codename = img.split(' ')[0].split('_')
                    print "img_codename", img_codename, type(img_codename)
                    color = matching_map[img_codename[-1]]
                    name_pattern = ' - '.join(prefix + img_codename[:-1])
                    print "name pattern", name_pattern, type(name_pattern)
                    print
                    print
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
            print key
            writer.writerow([key] + matched_images[key] + [0])
        initial_open.close()

    else:
        print "Continuing where the process last left off"
        file_exists = os.path.isfile("matched_images.tmp")

        if file_exists:

            open_tmp = open('matched_images.tmp', 'r')
            read_preexisting_tmp = csv.reader(open_tmp, delimiter='\t')
            tmp_list = []
            for row in read_preexisting_tmp:
                tmp_list.append(row)
            open_tmp.close()

            open_csv = open('matched_images.csv', 'r')
            read_preexisting_csv = csv.reader(open_csv, delimiter='\t')
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
            channels.append(cf.tiff_stack_2_np_arr(color))
        yield name_pattern, tags_dict[name_pattern], channels
        print row[3]
        row[3] = 1
        writer_check_tmp.writerow(row)



source = Kristen_traverse("/run/user/1000/gvfs/smb-share:server=10.17.0.219,share=common/Users/kristen/Split GFP quant_Andrei/", matching_map=translator)
print source
#
#     matched_images = defaultdict(lambda: ['']*len(matching_map.keys()))
#     tags_dict = defaultdict(lambda: [])
# #     using a dictionary instead of list since grouped by mutliple colors and may need to superimpose them
# #     provide matching map based on translator (a dictionary with color codes and associated colors/dyes)
#     image_list = {'B_1': [], 'B_2':[], 'B_3':[], 'B_4':[], 'B_5':[], 'B_6':[], 'B_7':[], 'C_2':[], 'C_3':[], 'C_4':[]}
#
#     for current_location, sub_directories, files in os.walk(main_root):
#         if files:
#             for img in files:
#                 if ('.TIF' in img or '.tif' in img) and '_thumb_' not in img:
#                     # img_codename = [img.split('.')[0]]   #linhao
#                     prefix = cf.split_and_trim(current_location, main_root)
#
#                     img_codename = [img.split('.')[0]]
#                     name_pattern = ' - '.join(prefix + img_codename)
#                     print "prefix", prefix
#                     print "img codename", img_codename
#                     print "name pattern", name_pattern
#                     group_by = img_codename[0][:2]
#                     # grouping by color
#                     print "group by", group_by
#                     if name_pattern[-3::] in image_list:
#                         image_list[name_pattern[-3::]].append(img_codename)
#                     # project the three channels
#                     color = group_by
#                     # matched_images[name_pattern][color] = os.path.join(current_location, img)
#                     print type(img)
#             print image_list
#
#     for name_pattern, group_by, image_location in matched_images:
#         print "reached for loop"
#         print name_pattern
#         stack = cf.tiff_stack_2_np_arr(image_location)
#         stack = np.rollaxis(stack[0, :, :, :], 2)  # on the data where channels have not been split into z stacks
#         channels = np.split(stack, stack.shape[0])
#         channels = [chan[0, :, :] for chan in channels]
#         pass
# #     dont forget to add the csv/user input portion after this part is running smoothly
# # need to add yield, but first fix matched images so it works properly
# #

# # C1 = DAPI, C2 = GFP, C3 = mCherry

