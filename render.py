import os
import numpy as np
from matplotlib import pyplot as plt
from core_functions import generator_wrapper
from csv import writer as csv_writer
import scipy


@generator_wrapper(in_dims=(None, 2, 2, 2, 2, 1, 1, None, None, 2), out_dims=(None,))
def linhao_gfp_render(name_pattern, proj_gfp, qual_gfp,                      cell_labels, average_gfp_pad, average_gfp_list,
                      predicted_average, std_err, upper_outliers, gfp_outliers,
                      save=False, directory_to_save_to='verification'):

    # To Consider: bind all the x/y axes together.

    plt.figure(figsize=(20, 15))

    plt.suptitle(name_pattern)

    plt.subplot(241)
    plt.imshow(proj_gfp, interpolation='nearest')

    plt.subplot(242)
    plt.imshow(np.log(proj_gfp + np.min(proj_gfp[proj_gfp > 0])), cmap='hot', interpolation='nearest')

    plt.subplot(243)
    plt.imshow(qual_gfp, cmap='gray', interpolation='nearest')

    plt.subplot(244)
    plt.imshow(cell_labels, cmap=plt.cm.spectral, interpolation='nearest')

    plt.subplot(245)
    plt.imshow(average_gfp_pad, cmap='hot', interpolation='nearest')
    plt.colorbar()

    plt.subplot(246)
    plt.plot(np.sort(average_gfp_list), 'ko')
    plt.plot(predicted_average, 'r')
    plt.plot(predicted_average + std_err, 'g')
    plt.plot(predicted_average - std_err, 'g')

    plt.subplot(247)
    plt.imshow(gfp_outliers, cmap='gray', interpolation='nearest')

    # plt.subplot(248)
    # plt.imshow(self.qualifying_gfp_mask)

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