import csv
import logging as logger
import os
from _csv import writer as csv_writer
from collections import defaultdict

import numpy as np
import scipy
from matplotlib import pyplot as plt

import imagepipe.wrapped_functions
from imagepipe import core_functions as cf
from imagepipe.core_functions import generator_wrapper


@generator_wrapper(in_dims=(None, 2, 2, 2, 2, 1, 1, None, None, 2, 2), out_dims=(None,))
def linhao_gfp_render(name_pattern, proj_gfp, qual_gfp, cell_labels, average_gfp_pad, average_gfp_list,
                      predicted_average, std_err, upper_outliers, gfp_outliers, proj_mch,
                      save=False, directory_to_save_to='verification'):

    # To Consider: bind all the x/y axes together.

    plt.figure(figsize=(20, 15))

    plt.suptitle(name_pattern)

    main_ax = plt.subplot(241)
    plt.imshow(proj_gfp, interpolation='nearest')
    plt.contour(cell_labels > 0, [0.5], colors='w')

    plt.subplot(242, sharex=main_ax, sharey=main_ax)
    plt.imshow(np.log(proj_gfp + np.min(proj_gfp[proj_gfp > 0])), cmap='hot', interpolation='nearest')
    plt.contour(cell_labels > 0, [0.5], colors='w')

    plt.subplot(243, sharex=main_ax, sharey=main_ax)
    plt.imshow(qual_gfp, cmap='gray', interpolation='nearest')
    plt.contour(cell_labels > 0, [0.5], colors='w')

    ax = plt.subplot(244, sharex=main_ax, sharey=main_ax)
    plt.imshow(cell_labels, cmap=plt.cm.spectral, interpolation='nearest')
    unique = np.sort(np.unique(cell_labels))[1:]
    for i in unique:
        mask = cell_labels == i
        x, y = scipy.ndimage.measurements.center_of_mass(mask)
        ax.text(y-8, x+8, '%s' % i, fontsize=10)

    plt.subplot(245, sharex=main_ax, sharey=main_ax)
    plt.imshow(average_gfp_pad, cmap='hot', interpolation='nearest')
    plt.colorbar()
    plt.contour(cell_labels > 0, [0.5], colors='w')

    plt.subplot(246)
    plt.plot(np.sort(average_gfp_list), 'ko')
    plt.plot(predicted_average, 'r')
    plt.plot(predicted_average + std_err, 'g')
    plt.plot(predicted_average - std_err, 'g')

    ax = plt.subplot(247, sharex=main_ax, sharey=main_ax)
    plt.imshow(gfp_outliers, cmap=plt.cm.spectral, interpolation='nearest')
    unique = np.unique(gfp_outliers)
    for i in unique:
        mask = gfp_outliers == i
        x, y = scipy.ndimage.measurements.center_of_mass(mask)
        ax.text(y-8, x+8, '%s' % i, fontsize=10)

    plt.subplot(248, sharex=main_ax, sharey=main_ax)
    plt.imshow(proj_gfp, interpolation='nearest')
    plt.contour(cell_labels > 0, [0.5], colors='w')

    if not save:
        plt.show()

    else:
        name_puck = directory_to_save_to+'/'+'gfp_base-'+name_pattern+'.png'
        plt.savefig(name_puck)
        plt.close()


@generator_wrapper(in_dims=(None, 2, 2, 2, 2, 2, 2, 2, 2), out_dims=(None,))
def linhao_mch_render(name_pattern, proj_mCh, mitochondria, skeleton,
                      mito_classes, final_classes, cell_labels,
                      radius_mask, support_mask,
                      save=False, directory_to_save_to='verification'):

    # To Consider: bind all the x/y axes together.

    cell_binary = cell_labels > 0
    plt.figure(figsize=(20.0, 15.0))

    plt.suptitle(name_pattern)

    plt.subplot(241)
    plt.imshow(proj_mCh, interpolation='nearest')
    plt.contour(cell_binary, [0.5], colors='k')

    plt.subplot(242)
    plt.imshow(np.log(proj_mCh + np.min(proj_mCh[proj_mCh > 0])), cmap='hot', interpolation='nearest')
    plt.contour(cell_binary, [0.5], colors='w')

    plt.subplot(243)
    plt.imshow(mitochondria, interpolation='nearest', cmap='gray')

    plt.subplot(244)
    plt.imshow(skeleton, interpolation='nearest', cmap=plt.cm.spectral)
    plt.contour(mitochondria, [0.5], colors='w')

    plt.subplot(245)
    plt.imshow(radius_mask, cmap='hot', interpolation='nearest')
    plt.contour(cell_binary, [0.5], colors='w')

    plt.subplot(246)
    plt.imshow(support_mask, cmap='hot', interpolation='nearest')
    plt.contour(cell_binary, [0.5], colors='w')

    plt.subplot(247)
    plt.imshow(mito_classes, interpolation='nearest', cmap='coolwarm')
    plt.contour(cell_binary, [0.5], colors='k')

    plt.subplot(248)
    plt.imshow(final_classes, interpolation='nearest', cmap='coolwarm')
    plt.contour(cell_binary, [0.5], colors='k')

    if not save:
        plt.show()

    else:
        name_puck = directory_to_save_to+'/'+'mCh_base-'+name_pattern+'.png'
        plt.savefig(name_puck)
        plt.close()


@generator_wrapper(in_dims=(None, 2, 2, 2, 2, 2, 2), out_dims=(None,))
def linhao_mqvi_render(name_pattern, mito_outlines, cell_labels,
                       projected_gfp, projected_mch,
                       gfp_mqvi, mch_mqvi,
                       save=False, directory_to_save_to='verification'):
    plt.figure(figsize=(26.0, 15.0))
    plt.suptitle(name_pattern)

    main_ax = plt.subplot(221)
    plt.title('mCherry')
    plt.imshow(np.log(projected_mch+np.min(projected_mch[projected_mch > 0])),
               interpolation='nearest', cmap='hot')
    plt.contour(mito_outlines, [0.5], colors='b')
    plt.contour(cell_labels, [0.5], colors='g')

    plt.subplot(222, sharex=main_ax, sharey=main_ax)
    plt.title('GFP')
    plt.imshow(np.log(projected_gfp+np.min(projected_gfp[projected_gfp > 0])),
               interpolation='nearest', cmap='hot')
    plt.contour(mito_outlines, [0.5], colors='b')
    plt.contour(cell_labels, [0.5], colors='g')

    plt.subplot(224, sharex=main_ax, sharey=main_ax)
    plt.title('GFP MQVI')
    plt.imshow(gfp_mqvi, interpolation='nearest', cmap='hot')
    plt.colorbar()

    plt.subplot(223, sharex=main_ax, sharey=main_ax)
    plt.title('mCherry MQVI')
    plt.imshow(mch_mqvi, interpolation='nearest', cmap='hot')
    plt.colorbar()


    if not save:
        plt.show()

    else:
        name_puck = directory_to_save_to+'/'+'mqvi-'+name_pattern+'.png'
        plt.savefig(name_puck)
        plt.close()


@generator_wrapper
def linhao_summarize(primary_namespace, output):
    with open(output, 'ab') as output_file:
        writer = csv_writer(output_file)
        namespace = primary_namespace['name pattern']
        tag_group = primary_namespace['group id']
        secondary_namespace = primary_namespace['per_cell']
        pre_puck = [namespace] + tag_group
        for key, value in secondary_namespace.iteritems():
            if key != '_pad':
                proper_puck = pre_puck+[key, value['gfp_mqvi'], value['mch_mqvi'], value['final_classification']]
                writer.writerow(proper_puck)

    return primary_namespace


@generator_wrapper
def linhao_secondary_summarize(primary_namespace, output):
    with open(output, 'ab') as output_file:
        writer = csv_writer(output_file)
        namespace = primary_namespace['name pattern']
        tag_group = primary_namespace['group id']
        total_cells = len(np.unique(primary_namespace['pre_cell_labels'])) - 1
        analyzed_cells = len(np.unique(primary_namespace['cell_labels'])) - 1
        writer.writerow([namespace] + tag_group + [total_cells, analyzed_cells])

    return primary_namespace


def Linhao_traverse(main_root,
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
    tags_dict = defaultdict(lambda: [])
    if matching_rule:
        assert(matching_map is not None)

    for current_location, sub_directories, files in os.walk(main_root):

        if files:

            for img in files:
                if ('.TIF' in img or '.tif' in img) and '_thumb_' not in img and 'w' in img:

                    prefix = imagepipe.wrapped_functions.split_and_trim(current_location, main_root)

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
            channels.append(imagepipe.wrapped_functions.tiff_stack_2_np_arr(color))
        print "name pattern", name_pattern
        print 'channels', channels
        yield name_pattern, tags_dict[name_pattern], channels
        row[3] = 1
        writer_check_tmp.writerow(row)