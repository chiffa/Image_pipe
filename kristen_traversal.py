import numpy as np
import core_functions as cf
from collections import defaultdict
import logging as logger
import os
from pympler import  muppy, summary
import csv
from matplotlib import pyplot as plt


translator = {'C1':0,
              'C3':1,
              'C4':2}


def Kristen_traverse(main_root, matching_rule='c', matching_map=None):
    print "starting kristen's traversal"
    matched_images = defaultdict(lambda: [''] * len(matching_map))
    # reference: {name_pattern:[location of DAPI, location of GFP, location of mCherry]}
    print len(matching_map)
    print matched_images
    tags_dict = defaultdict(lambda: [])
    # do we even need this? For linhao's case this was used to keep track of HS time
    name_pattern_list = []
    if matching_rule:
        assert (matching_map is not None)
    for current_location, sub_directories, files in os.walk(main_root):
            if files:
                for img in files:
                    if ('.TIF' in img or '.tif' in img) and '_thumb_' not in img:
                        prefix = cf.split_and_trim(current_location, main_root)
                        img_codename = [img.split('.')[0]]
                        print img_codename[0].split(' ')[1:]
                        name_pattern = ' - '.join(prefix + img_codename[0].split(' ')[1:])
                        print "prefix", prefix
                        print "img codename", img_codename
                        print "name pattern", name_pattern, type(name_pattern)
                        group_by = img_codename[0][:2]
                        color = matching_map[img_codename[0].split('-')[0]]

                        print "group by", group_by
                        print "color", color, type(color)
                        matched_images[name_pattern][color] = os.path.join(current_location, img)

                        print
                        print
    for name in matched_images:
        name_pattern_list.append(name_pattern)
        plt.figure(figsize=(20.0, 15.0))
        plt.suptitle('Projected DAPI. GFP, mCherry')
        main_ax = plt.subplot(121)
        plt.title('DAPI')
        plt.imshow(name[0], interpolation='nearest', cmap=plt.cm.spectral)
        cbar = plt.colorbar()
        plt.title('GFP')
        plt.imshow(name[1], interpolation='nearest', cmap='gray')
        plt.imshow(name[2], interpolation='nearest', alpha=0.3)
    print "here"
    print set(name_pattern_list)

    for key in matched_images:
        print key
        print matched_images[key]

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

    # for row in csv_reader:
    #     name_pattern = row[0]
    #     color_set = [row[1], row[2]]
    #     if row[3] == 1:
    #         writer_check_tmp.writerow(row)
    #         continue
    #     channels = []
    #     for color in color_set:
    #         channels.append(cf.tiff_stack_2_np_arr(color))
    #     yield name_pattern, tags_dict[name_pattern], channels
    #     print row[3]
    #     row[3] = 1
    #     writer_check_tmp.writerow(row)


Kristen_traverse("/run/user/1000/gvfs/smb-share:server=10.17.0.219,share=common/Users/kristen/Split GFP quant_Andrei/", matching_map=translator)


# # C1 = DAPI, C3 = GFP, C4 = mCherry

