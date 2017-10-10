import os
from _csv import writer as csv_writer

import numpy as np
import scipy
from matplotlib import pyplot as plt

from imagepipe.raw_functions import split_and_trim
from imagepipe.tools.helpers import tiff_stack_2_np_arr
import imagepipe.wrapped_functions
from imagepipe import core_functions as cf
from imagepipe.core_functions import generator_wrapper


@generator_wrapper(in_dims=(None, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2), out_dims=(None,))
def akshay_render(name_pattern, DAPI, p53, p21,
                  nuclei, vor_segment,
                  extra_nuclear_p53, nuclear_p53_pad, extranuclear_p53_pad,
                  extra_nuclear_p21, nuclear_p21_pad, extranuclear_p21_pad,
                  save=False, directory_to_save_to='verification'):

    plt.figure(figsize=(26.0, 15.0))
    plt.suptitle(name_pattern)

    main_ax = plt.subplot(231)
    plt.title('DAPI')
    plt.imshow(DAPI, interpolation='nearest')
    plt.contour(nuclei, [0.5], colors='k')

    plt.subplot(232, sharex=main_ax, sharey=main_ax)
    plt.title('p53')
    plt.imshow(p53, interpolation='nearest')
    plt.contour(nuclei, [0.5], colors='k')
    plt.contour(extra_nuclear_p53, [0.5], colors='w')

    plt.subplot(233, sharex=main_ax, sharey=main_ax)
    plt.title('p21')
    plt.imshow(p21, interpolation='nearest')
    plt.contour(nuclei, [0.5], colors='k')
    plt.contour(extra_nuclear_p21, [0.5], colors='w')

    ax = plt.subplot(234, sharex=main_ax, sharey=main_ax)
    plt.title('nuclei & Voronoi segmentation')
    plt.imshow(vor_segment, interpolation='nearest', cmap='spectral', vmin=0)
    plt.contour(nuclei, [0.5], colors='k')
    unique = np.unique(vor_segment)
    for i in unique:
        mask = nuclei == i
        x, y = scipy.ndimage.measurements.center_of_mass(mask)
        ax.text(y-8, x+8, '%s' % i, fontsize=10)

    plt.subplot(235, sharex=main_ax, sharey=main_ax)
    plt.title('p53 nucleus/cell intensity')
    p_53_summmary = np.zeros_like(nuclear_p53_pad)
    p_53_summmary[extranuclear_p53_pad > 0] = extranuclear_p53_pad[extranuclear_p53_pad > 0]
    p_53_summmary[nuclear_p53_pad > 0] = nuclear_p53_pad[nuclear_p53_pad > 0]
    im = plt.imshow(p_53_summmary, interpolation='nearest', cmap='hot')
    plt.colorbar(im)
    plt.contour(nuclei, [0.5], colors='b')
    plt.contour(extra_nuclear_p53, [0.5], colors='g')

    plt.subplot(236, sharex=main_ax, sharey=main_ax)
    plt.title('p21 nucleus/cell intensity')
    p_21_summmary = np.zeros_like(nuclear_p21_pad)
    p_21_summmary[extranuclear_p21_pad > 0] = extranuclear_p21_pad[extranuclear_p21_pad > 0]
    p_21_summmary[nuclear_p21_pad > 0] = nuclear_p21_pad[nuclear_p21_pad > 0]
    im = plt.imshow(p_21_summmary, interpolation='nearest', cmap='hot')
    plt.colorbar(im)
    plt.contour(nuclei, [0.5], colors='b')
    plt.contour(extra_nuclear_p21, [0.5], colors='g')

    if not save:
        plt.show()

    else:
        name_puck = directory_to_save_to+'/'+'akshay-'+name_pattern+'.png'
        plt.savefig(name_puck)
        plt.close()


@generator_wrapper(in_dims=(None, None, 1, 1, 1, 1), out_dims=(None,))
def akshay_summarize(name_pattern, group_by, av_nuc_p53, av_en_p53, av_nuc_p21, av_en_p21, output):
    with open(output, 'ab') as output_file:
        writer = csv_writer(output_file)
        for i, nuc_pac in enumerate(zip(av_nuc_p53, av_en_p53, av_nuc_p21, av_en_p21)):
            writer.writerow([name_pattern, group_by, i, nuc_pac[0], nuc_pac[1], nuc_pac[2], nuc_pac[3]])


def Akshay_traverse(main_root):

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
