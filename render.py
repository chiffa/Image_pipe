import os
import numpy as np
from matplotlib import pyplot as plt
from core_functions import generator_wrapper, safe_dir_create
from csv import writer as csv_writer
import scipy
from scipy import ndimage as ndi
from scipy import stats
import density_plot as dplt
import core_functions as cf
from chiffatools.dataviz import better2D_desisty_plot
from scipy.stats import linregress


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

            GFP_limited_to_cell_mask = cf._3d_stack_2d_filter(GFP_orig, my_mask)
            mCherry_limited_to_cell_mask = cf._3d_stack_2d_filter(mCherry_orig, my_mask)

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

            GFP_orig_qualifying = cf._3d_stack_2d_filter(GFP_orig, my_mask)
            mCherry_orig_qualifying = cf._3d_stack_2d_filter(mCherry_orig, my_mask)
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
            dplt.better2D_desisty_plot(GFP_1d, mCherry_1d)
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



@generator_wrapper(in_dims=(None, None, 1, 1, 1, 1), out_dims=(None,))
def akshay_summarize(name_pattern, group_by, av_nuc_p53, av_en_p53, av_nuc_p21, av_en_p21, output):
    with open(output, 'ab') as output_file:
        writer = csv_writer(output_file)
        for i, nuc_pac in enumerate(zip(av_nuc_p53, av_en_p53, av_nuc_p21, av_en_p21)):
            writer.writerow([name_pattern, group_by, i, nuc_pac[0], nuc_pac[1], nuc_pac[2], nuc_pac[3]])


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

# @generator_wrapper(in_dims=(None, 0, 0, 0, 0, 0, 0, 0, 0, 0 ), out_dims=(None,))
# def Kristen_summarize_a(name_pattern, q_mean,q_median, q_std, nq_mean, nq_median, nq_std, slope, r2, p, output):
#
#     with open(output, 'ab') as output_file:
#         # csv_read = csv.res
#         writer = csv_writer(output_file, delimiter = '\t')
#         for i in name_pattern:
#             writer.writerow([i, q_mean, q_median, q_std, nq_mean, nq_median, nq_std, slope, r2, p])

safe_dir_create('verification')
