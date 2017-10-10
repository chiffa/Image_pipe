import csv
import logging as logger
import os
from _csv import writer as csv_writer
from collections import defaultdict

import numpy as np
from matplotlib import pyplot as plt
from scipy import ndimage as ndi, stats

import imagepipe.raw_functions
import imagepipe.tools.helpers
import imagepipe.wrapped_functions
import imagepipe.core_functions
from imagepipe import core_functions as cf, density_plot
from imagepipe.core_functions import generator_wrapper


@generator_wrapper(in_dims=2,out_dims=(None,))
def Kristen_render_single_image(dapi, gfp, mcherry):
    plt.figure(figsize=(26.0, 15.0))
    plt.title('Max Projection')

    plt.subplot(221)
    plt.title('DAPI')
    plt.imshow(dapi,interpolation='nearest')

    plt.subplot(222)
    plt.title('GFP')
    plt.imshow(gfp, interpolation='nearest')

    plt.subplot(221)
    plt.title('mCherry')
    plt.imshow(mcherry, interpolation='nearest')


@generator_wrapper(in_dims=(None,None, 2, 2, 3, 3), out_dims=(None,))
def Kristen_render(name_pattern,
                   group_id,
                   mCherry,
                   extranuclear_mCherry_pad,
                   GFP_orig,
                   mCherry_orig, output,
                   save=False, directory_to_save_to='verification'):
    labels, _ = ndi.label(extranuclear_mCherry_pad)
    unique_segmented_cells_labels = np.unique(labels)[1:]
    mCherry_cutoff = np.zeros_like(mCherry)
    qualifying_cell_label = []
    qualifying_regression_stats = []

    for cell_label in unique_segmented_cells_labels:
        mCherry_2 = np.zeros_like(mCherry)
        my_mask = labels == cell_label
        average_apply_mask = np.mean(mCherry[my_mask])
        intensity = np.sum(mCherry[my_mask])
        binary_pad = np.zeros_like(mCherry)
        binary_pad[my_mask] = 1
        pixel = np.sum(binary_pad[my_mask])

        if (average_apply_mask > .05 or intensity > 300) and pixel > 4000:

            GFP_limited_to_cell_mask = imagepipe.raw_functions.f_3d_stack_2d_filter(GFP_orig, my_mask)
            mCherry_limited_to_cell_mask = imagepipe.raw_functions.f_3d_stack_2d_filter(mCherry_orig, my_mask)

            qualifying_3d_GFP = GFP_limited_to_cell_mask[mCherry_limited_to_cell_mask>50]
            average_3d_GFP = np.mean(qualifying_3d_GFP)
            median_3d_GFP = np.median(qualifying_3d_GFP)
            std_3d_GFP = np.std(qualifying_3d_GFP)
            sum_qualifying_GFP = np.sum(qualifying_3d_GFP)

            nonqualifying_3d_GFP = GFP_limited_to_cell_mask[mCherry_limited_to_cell_mask<=50]
            average_nonqualifying_3d_GFP = np.mean(nonqualifying_3d_GFP)
            median_nonqualifying_3d_GFP = np.median(nonqualifying_3d_GFP)
            std_nonqualifying_3d_GFP = np.std(nonqualifying_3d_GFP)
            sum_nonqualifying_GFP = np.sum(nonqualifying_3d_GFP)

            sum_total_GFP = sum_qualifying_GFP + sum_nonqualifying_GFP
            percent_qualifying_over_total_GFP = sum_qualifying_GFP/sum_total_GFP
            # report the percentage too or sums are sufficient?

            GFP_orig_qualifying = imagepipe.raw_functions.f_3d_stack_2d_filter(GFP_orig, my_mask)
            mCherry_orig_qualifying = imagepipe.raw_functions.f_3d_stack_2d_filter(mCherry_orig, my_mask)
            mCherry_1d = mCherry_orig_qualifying[mCherry_orig_qualifying > 50]
            GFP_1d = GFP_orig_qualifying[mCherry_orig_qualifying>50]
            regression_results = stats.linregress(GFP_1d, mCherry_1d)

            mCherry_2[my_mask] = mCherry[my_mask]
            mCherry_cutoff[my_mask] = mCherry[my_mask]
            qualifying_cell_label.append(cell_label)
            qualifying_regression_stats.append((regression_results[0], regression_results[2], regression_results[3]))

            name_pattern_split = name_pattern.split(' - ')
            transfection_label = name_pattern_split[0]
            cell_type = name_pattern_split[1]
            exp_time = name_pattern_split[2]
            image_number = name_pattern_split[4]

            with open(output, 'ab') as output_file:
                writer = csv_writer(output_file, delimiter='\t')
                writer.writerow([transfection_label, cell_type, exp_time, image_number, cell_label, sum_qualifying_GFP, sum_total_GFP, average_3d_GFP, median_3d_GFP, std_3d_GFP, average_nonqualifying_3d_GFP, median_nonqualifying_3d_GFP, std_nonqualifying_3d_GFP, regression_results[0], regression_results[2], regression_results[3]])

            plt.figure(figsize=(26.0, 15.0))
            plt.title('Kristen\'s Data')
            plt.suptitle(name_pattern)

            main_ax = plt.subplot(221)
            plt.subplot(221, sharex=main_ax, sharey=main_ax)
            plt.title('mCherry Binary')
            im = plt.imshow(extranuclear_mCherry_pad, interpolation='nearest', cmap = 'hot')
            plt.colorbar(im)
            plt.subplot(222, sharex=main_ax, sharey=main_ax)
            plt.title('mCherry')
            plt.imshow(mCherry, interpolation='nearest')
            plt.contour(extranuclear_mCherry_pad, [0.5], colors='k')
            plt.subplot(223)
            plt.better2D_desisty_plot(GFP_1d, mCherry_1d)
            plt.title('mCherry Intensity as a Function of GFP Voxel')
            plt.xlabel('GFP Voxel')
            plt.ylabel('mCherry Intensity')
            plt.subplot(224, sharex=main_ax, sharey=main_ax)
            plt.title('mCherry-cutoff applied')
            plt.imshow(mCherry_2, interpolation='nearest')

            if not save:
                plt.show()

            else:
                name_puck = directory_to_save_to + '/' + 'Kristen-' + name_pattern+ '_cell' + str(cell_label)+ '.png'
                plt.savefig(name_puck)
                plt.close()
    plt.figure(figsize=(26.0, 15.0))
    main_ax = plt.subplot(121)
    plt.subplot(121, sharex=main_ax, sharey=main_ax)
    plt.suptitle('mCherry Before and After Qualifying Cell Cutoff is Applied')
    plt.title('mCherry')
    im = plt.imshow(mCherry, interpolation='nearest')
    plt.colorbar(im)
    plt.subplot(122, sharex=main_ax, sharey=main_ax)
    plt.title('mCherry')
    plt.imshow(mCherry_cutoff, interpolation='nearest')
    if not save:
        plt.show()

    else:
        name_puck = directory_to_save_to + '/' + 'Kristen-' + name_pattern + 'cutoff_app' + '.png'
        plt.savefig(name_puck)
        plt.close()

    return qualifying_regression_stats


@generator_wrapper(in_dims=(None, 0, 0, 0, 0, 0, 0, 0, 0, 0 ), out_dims=(None,))
def Kristen_summarize_a(name_pattern, q_mean,q_median, q_std, nq_mean, nq_median, nq_std, slope, r2, p, output):

    with open(output, 'ab') as output_file:
        # csv_read = csv.res
        writer = csv_writer(output_file, delimiter = '\t')
        for i in name_pattern:
            writer.writerow([i, q_mean, q_median, q_std, nq_mean, nq_median, nq_std, slope, r2, p])


def Kristen_traverse(main_root, matching_rule='c', matching_map=None):
    print "starting kristen's traversal"
    matched_images = defaultdict(lambda: [''] * len(matching_map))

    # name_pattern_list = []
    if matching_rule:
        assert (matching_map is not None)
    for current_location, sub_directories, files in os.walk(main_root):
            if files:
                for img in files:
                    if ('.TIF' in img or '.tif' in img) and '_thumb_' not in img:
                        prefix = imagepipe.raw_functions.split_and_trim(current_location, main_root)
                        img_codename = [img.split('.')[0]]



                        # # choosing one image to work with
                        # c = img_codename[0].split('-')[0]
                        # b = img_codename[0].split(' ')[-1][0]
                        # print c
                        # print b
                        # if c == 'C1' and b[0] == 'B':
                        #     # change these conditions back to original to test all images
                        #     print "found image"
                        #
                        #     name_pattern = ' - '.join(prefix + img_codename[0].split(' ')[1:])
                        #     group_by = img_codename[0][:2]
                        #     color = matching_map[img_codename[0].split('-')[0]]
                        #     # print matched_images[name_pattern][color]
                        #     # print os.path.join(current_location, img)
                        #     matched_images[name_pattern][color] = os.path.join(current_location, img)
                    name_pattern = ' - '.join(prefix + img_codename[0].split(' ')[1:])
                    group_by = img_codename[0][:2]
                    color = matching_map[img_codename[0].split('-')[0]]
                    # print matched_images[name_pattern][color]
                    # print os.path.join(current_location, img)
                    matched_images[name_pattern][color] = os.path.join(current_location, img)

    # shift tab upper portion out/ placed inside for loop to study a single image but originally only inside the if(.TIF...)

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
        color_set = [row[1], row[2], row[3]]
        if row[3] == 1:
            writer_check_tmp.writerow(row)
            continue
        channels = []
        plot_list = []
        for color in color_set:
            channels.append(imagepipe.tools.helpers.tiff_stack_2_np_arr(color))
            plot_list.append(imagepipe.tools.helpers.tiff_stack_2_np_arr(color))

#         plt.figure(figsize=(20.0, 15.0))
#         plt.suptitle('Projected DAPI. GFP, mCherry')
#         plt.title('DAPI')
#         dapi = np.max(plot_list[0], axis=0)
#         plt.imshow(dapi, interpolation='nearest', cmap='gray')
#         plt.title('GFP')
#         gfp = np.max(plot_list[1], axis=0)
#         plt.imshow(gfp, interpolation='nearest', cmap='gray')
#         mcherry = np.max(plot_list[2], axis=0)
#         plt.imshow(mcherry, interpolation='nearest', alpha=0.3)
#         plt.show()

        yield name_pattern, matched_images, channels
        row[3] = 1
        writer_check_tmp.writerow(row)