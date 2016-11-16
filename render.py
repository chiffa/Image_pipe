import os
import numpy as np
from matplotlib import pyplot as plt
from core_functions import generator_wrapper


@generator_wrapper(in_dims=(None, 2, 2, 2, 2, 1, 1, None, None, 2), out_dims=(None,))
def linhao_gfp_render(name_pattern, proj_gfp, qual_gfp,
                      cell_labels, average_gfp_pad, average_gfp_list,
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


@generator_wrapper(in_dims=(None, 2, 2, 2, 2, 2, 2, 2), out_dims=(None,))
def akshay_render(name_pattern, DAPI, p53, p21, nuclei, p53_pad, p53_o_n, p53_o_n_seg):

    plt.figure(figsize=(20.0, 15.0))
    plt.suptitle(name_pattern)

    main_ax = plt.subplot(241)
    plt.title('DAPI')
    plt.imshow(DAPI, interpolation='nearest')
    plt.contour(nuclei, [0.5], colors='k')

    plt.subplot(242, sharex=main_ax, sharey=main_ax)
    plt.title('p53')
    plt.imshow(p53, interpolation='nearest')
    plt.contour(nuclei, [0.5], colors='k')

    plt.subplot(243, sharex=main_ax, sharey=main_ax)
    plt.title('p21')
    plt.imshow(p21, interpolation='nearest')
    plt.contour(nuclei, [0.5], colors='k')

    plt.subplot(244, sharex=main_ax, sharey=main_ax)
    plt.title('p53 no nucleus')
    plt.imshow(p53_o_n, cmap='hot', interpolation='nearest')
    plt.contour(nuclei, [0.5], colors='k')

    plt.subplot(245, sharex=main_ax, sharey=main_ax)
    plt.title('nuclei')
    plt.imshow(nuclei, interpolation='nearest', cmap='spectral')

    plt.subplot(246, sharex=main_ax, sharey=main_ax)
    plt.title('p53 in nucleus average intensity')
    plt.imshow(p53_pad, interpolation='nearest', cmap='hot')
    plt.contour(nuclei, [0.5], colors='w')
    plt.contour(p53_o_n_seg, [0.5], colors='w')

    plt.subplot(248, sharex=main_ax, sharey=main_ax)
    plt.title('p53 no nucleus segmented')
    plt.imshow(p53_o_n_seg, interpolation='nearest', cmap='gray')
    plt.contour(nuclei, [0.5], colors='k')

    # plt.subplot(245, sharex=main_ax, sharey=main_ax)
    # plt.title('log-DAPI')
    # plt.imshow(np.log(DAPI + np.min(DAPI[DAPI > 0])), cmap='hot', interpolation='nearest')
    # plt.contour(nuclei, [0.5], colors='w')

    # plt.subplot(245)
    # plt.hist(DAPI.flatten(), bins=100)
    #
    # plt.subplot(246)
    # plt.hist(p53.flatten(), bins=100)
    #
    # plt.subplot(247)
    # plt.hist(p21.flatten(), bins=100)

    plt.show()
