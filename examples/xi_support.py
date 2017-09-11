import os
from _csv import writer as csv_writer
from collections import defaultdict

import numpy as np
import scipy
from chiffatools.dataviz import better2D_desisty_plot
from matplotlib import pyplot as plt
from scipy.stats import linregress

from imagepipe import core_functions as cf
from imagepipe.core_functions import generator_wrapper


def xi_traverse(main_root, matching_map=None):
    """
    Traverses the main_root directory, looking for all the '.tif/.TIF' files, performs name matching
    then iterates through the resulting matched dictironary.

    Matching assumption is that except for the matching keys, the names are identical

    :param main_root: folder from which will be traversed in depth
    :param matching_rule: name modification to type mapping. Currently '' for no matching, 'color' for colors
    :param matching_map: {'pattern in the file name': color channel number}
    :return:
    """
    matched_images = defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))
    tags_dict = defaultdict(lambda: [])

    assert(matching_map is not None)

    for current_location, sub_directories, files in os.walk(main_root):
        for img in files:
            if ('.TIF' in img or '.tif' in img) and '_thumb_' not in img:
                prefix = cf.split_and_trim(current_location, main_root)

                pre_name = '-'.join(img.split('-')[1:])[:-4]
                # print pre_name[-10:]
                _time, _z = pre_name[-9:].split('_')
                time_stamp = int(_time[1:])
                z_position = int(_z[1:])
                color = matching_map[img.split('-')[0]]
                name_pattern = ' - '.join(prefix + [pre_name[:-10]])
                matched_images[name_pattern][time_stamp][color][z_position] = os.path.join(current_location, img)
                # print name_pattern
                # print time_stamp, color, z_position

    for name_pattern, time_dict in matched_images.iteritems():
        for time_stamp, color_dict in time_dict.iteritems():
            channels = ['', '']
            for color, z_position_dict in color_dict.iteritems():
                z_collector = []
                for z_position, file_name in sorted(z_position_dict.items()):
                    z_collector.append(cf.tiff_stack_2_np_arr(file_name)[0, :, :])
                channels[color] = np.array(z_collector)

            yield name_pattern, str(time_stamp), channels


@generator_wrapper(in_dims=(None, 2, 2, 2, 2, 2, 3, 3, None), out_dims=(None,))
def xi_pre_render(name_pattern, proj_gfp, qual_gfp, cell_labels, average_gfp_pad, proj_mch,
                  mch, gfp, timestamp,
                  save=False, directory_to_save_to='verification', mch_cutoff=0.2, slector_cutoff=0.1):

    plt.figure(figsize=(20, 15))

    plt.suptitle(name_pattern)

    main_ax = plt.subplot(231)
    plt.title('GFP')
    plt.imshow(proj_gfp, interpolation='nearest')
    plt.contour(cell_labels > 0, [0.5], colors='w')

    plt.subplot(232, sharex=main_ax, sharey=main_ax)
    plt.title('log-GFP')
    plt.imshow(np.log(proj_gfp + np.min(proj_gfp[proj_gfp > 0])), cmap='hot', interpolation='nearest')
    plt.contour(cell_labels > 0, [0.5], colors='w')

    plt.subplot(233, sharex=main_ax, sharey=main_ax)
    plt.title('raw segmentation')
    plt.imshow(qual_gfp, cmap='gray', interpolation='nearest')
    plt.contour(cell_labels > 0, [0.5], colors='w')

    ax = plt.subplot(234, sharex=main_ax, sharey=main_ax)
    plt.title('labeled segmentation')
    plt.imshow(cell_labels, cmap=plt.cm.spectral, interpolation='nearest')
    unique = np.unique(cell_labels)
    for i in unique:
        mask = cell_labels == i
        x, y = scipy.ndimage.measurements.center_of_mass(mask)
        ax.text(y-8, x+8, '%s' % i, fontsize=10)

    plt.subplot(235)
    selector = np.logical_and(mch > slector_cutoff, gfp > slector_cutoff)
    plt.title('mCh-GFP correlation - %s, qual GFP intensity: %s' %
              (np.corrcoef(mch[selector], gfp[selector])[0, 1], np.median(gfp[mch > mch_cutoff])))
    slope, intercept, rvalue, pvalue, stderr = linregress(mch[selector], gfp[selector])
    better2D_desisty_plot(mch[selector], gfp[selector])
    linarray = np.arange(0.1, 0.5, 0.05)
    plt.plot(linarray, intercept+slope*linarray, 'r')
    plt.xlabel('mCherry')
    plt.ylabel('GFP')

    plt.subplot(236, sharex=main_ax, sharey=main_ax)
    plt.title('mCherry')
    plt.imshow(proj_mch, interpolation='nearest')
    plt.contour(cell_labels > 0, [0.5], colors='w')

    with open('xi_analys_results.csv', 'ab') as output_file:
        writer = csv_writer(output_file)

        puck = [name_pattern, timestamp,
                np.corrcoef(mch[selector], gfp[selector])[0, 1],
                np.median(gfp[mch > mch_cutoff]), np.average(gfp[mch > mch_cutoff]),
                slope, rvalue, pvalue]
        writer.writerow(puck)

    if not save:
        plt.show()

    else:
        name_puck = directory_to_save_to+'/'+'xi_pre_render-'+timestamp+'-'+name_pattern+'.png'
        plt.savefig(name_puck)
        plt.close()