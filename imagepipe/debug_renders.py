from matplotlib import pyplot as plt
import numpy as np
from scipy import histogram2d
from core_functions import watershed
from scipy import ndimage as ndi
# my addition to imports: watershed, ndi

def robust_binarize_debug(base, smooth, median, otsu, labels, binary, u_median, u_otsu):
    plt.figure(figsize=(20.0, 15.0))
    plt.suptitle('Robust Binarize_debug')

    main_ax = plt.subplot(241)
    plt.title('base')
    plt.imshow(base, interpolation='nearest')
    plt.contour(binary, [0.5], colors='k')

    plt.subplot(242, sharex=main_ax, sharey=main_ax)
    plt.title('smooth')
    plt.imshow(smooth, interpolation='nearest')
    plt.contour(binary, [0.5], colors='k')

    plt.subplot(243, sharex=main_ax, sharey=main_ax)
    plt.title('median')
    plt.imshow(median, interpolation='nearest')
    plt.contour(binary, [0.5], colors='k')

    plt.subplot(244, sharex=main_ax, sharey=main_ax)
    plt.title('otsu')
    plt.imshow(otsu, interpolation='nearest')
    plt.contour(binary, [0.5], colors='k')

    plt.subplot(245, sharex=main_ax, sharey=main_ax)
    plt.title('otsu - median')
    plt.imshow(otsu-median, interpolation='nearest', cmap='coolwarm')
    plt.colorbar()
    plt.contour(binary, [0.5], colors='k')

    plt.subplot(246, sharex=main_ax, sharey=main_ax)
    plt.title('otsu - unsmoothed median')
    plt.imshow(u_otsu - u_median, interpolation='nearest')
    plt.contour(binary, [0.5], colors='k')

    plt.subplot(247, sharex=main_ax, sharey=main_ax)
    plt.title('labels')
    plt.imshow(labels, interpolation='nearest')
    plt.contour(binary, [0.5], colors='k')

    plt.subplot(248, sharex=main_ax, sharey=main_ax)
    plt.title('binary')
    plt.imshow(binary, interpolation='nearest')


def voronoi_debug(base, points, labels, watershed):
    plt.figure(figsize=(20.0, 15.0))
    plt.suptitle('Voronoi debug')

    main_ax = plt.subplot(221)
    plt.title('base')
    plt.imshow(base, interpolation='nearest')

    plt.subplot(222, sharex=main_ax, sharey=main_ax)
    plt.title('points')
    plt.imshow(points, interpolation='nearest', cmap='hot')

    plt.subplot(223, sharex=main_ax, sharey=main_ax)
    plt.title('labels')
    plt.imshow(labels, interpolation='nearest', cmap='hot')

    plt.subplot(224, sharex=main_ax, sharey=main_ax)
    plt.title('watershed')
    plt.imshow(watershed, interpolation='nearest', cmap=plt.cm.spectral, vmin=0)


def filter_labels_debug(labels, binary, result):
    plt.figure(figsize=(20.0, 15.0))
    plt.suptitle('filter labels')

    main_ax = plt.subplot(221)
    plt.title('labels')
    plt.imshow(labels, interpolation='nearest', cmap=plt.cm.spectral)

    plt.subplot(222, sharex=main_ax, sharey=main_ax)
    plt.title('binary')
    plt.imshow(binary, interpolation='nearest', cmap='gray')

    plt.subplot(223, sharex=main_ax, sharey=main_ax)
    plt.title('result')
    plt.imshow(result, interpolation='nearest', cmap=plt.cm.spectral, vmin=0)


def weight_sum_zero_debug(label_mask, img_codename):
    plt.figure(figsize=(20.0, 15.0))
    plt.suptitle('weight sum zero debugger')

    main_ax = plt.subplot(121)
    plt.title('label_mask')
    plt.imshow(label_mask, interpolation='nearest', cmap=plt.cm.spectral)

    plt.subplot(122, sharex=main_ax, sharey=main_ax)
    plt.title('image code')
    plt.imshow(img_codename, interpolation='nearest', cmap='gray')


def skeleton_debug(label_mask, skeleton, float_surface):
    plt.figure(figsize=(20.0, 15.0))
    plt.suptitle('skeleton debugger')

    main_ax = plt.subplot(221)
    plt.title('label_mask')
    plt.imshow(label_mask, interpolation='nearest', cmap=plt.cm.spectral)

    plt.subplot(222, sharex=main_ax, sharey=main_ax)
    plt.title('skeleton')
    plt.imshow(skeleton, interpolation='nearest', cmap='gray')

    plt.subplot(223, sharex=main_ax, sharey=main_ax)
    plt.title('float surface')
    plt.imshow(float_surface, interpolation='nearest', cmap='gray')


def max_projection_debug(current_image):
    plt.figure(figsize=(20.0, 15.0))
    plt.suptitle('current image')
    main_ax = plt.subplot(121)
    plt.title('current')
    plt.imshow(current_image, interpolation='nearest', cmap=plt.cm.spectral)
    plt.show()


def weight_sum_debug_see_full_image(mitochondria, proj_mCh, skeleton, radius_mask, support_mask, mito_classes, final_classes, cell_labels, name_pattern):

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


# def robust_binarize_debug(base_image, binary_labels):
#     plt.figure(figsize=(20.0, 15.0))
#     plt.suptitle('robust binarize')
#
#     main_ax = plt.subplot(121)
#     plt.title('base image')
#     plt.imshow(base_image, interpolation='nearest', cmap=plt.cm.spectral)
#
#     plt.subplot(122, sharex=main_ax, sharey=main_ax)
#     plt.title('binary labels')
#     plt.imshow(binary_labels, interpolation='nearest', cmap='gray')


def improved_watershed_debug(segmented_cells_labels, mcherry):
    plt.figure(figsize=(20.0, 15.0))
    plt.suptitle('improved watershed')

    main_ax = plt.subplot(121)
    plt.title('segmented labels')
    plt.imshow(segmented_cells_labels, interpolation='nearest', cmap=plt.cm.spectral)
    cbar = plt.colorbar()
    plt.subplot(122, sharex=main_ax, sharey=main_ax)
    plt.title('projected mCherry')
    plt.imshow(mcherry, interpolation='nearest', cmap='gray')
    plt.imshow(segmented_cells_labels, interpolation='nearest', alpha= 0.3)


    plt.savefig('watershed_debug_image')


def improved_watershed_plot_intensities(unique_segmented_cell_labels, average_apply_mask):
    plt.figure()
    plt.title("Average Intensity as a Function of Cell Number")
    plt.xlabel('Cell Number')
    plt.ylabel('Average Intensity (pixels)')
    ax = plt.axes()
    x = [0,1,2,3,4,5,6,7,8,9,10,11,12,13]
    y = average_apply_mask
    ax.plot(y)
    my_xticks = unique_segmented_cell_labels
    plt.xticks(x, my_xticks)


def plot_GFP_as_a_function_of_cell_number(average_apply_mask):
    plt.figure()
    plt.title("GFP as a Function of Cell Number")
    plt.xlabel('Cell Number')
    plt.ylabel('GFP')
    ax = plt.axes()
    y = average_apply_mask
    ax.plot(y)
    plt.show()


def label_based_aq(ar):
    plt.figure()
    plt.title("GFP as a Function of Cell Number")
    plt.xlabel('Cell Number')
    plt.ylabel('GFP')
    ax = plt.axes()
    x = np.array[1:33]
    y = ar
    ax.plot(y)
    plt.show()


def better2D_desisty_plot(xdat, ydat, thresh=3, bins=(100, 100)):
    xdat = np.max(xdat, axis=0)
    ydat = np.max(ydat, axis=0)
    xyrange = [[min(xdat), max(xdat)], [min(ydat), max(ydat)]]
    distortion = (xyrange[1][1] - xyrange[1][0]) / \
        (xyrange[0][1] - xyrange[0][0])
    xdat = xdat * distortion

    xyrange = [[min(xdat), max(xdat)], [min(ydat), max(ydat)]]
    hh, locx, locy = histogram2d(xdat, ydat, range=xyrange, bins=bins)
    posx = np.digitize(xdat, locx)
    posy = np.digitize(ydat, locy)

    ind = (posx > 0) & (posx <= bins[0]) & (posy > 0) & (posy <= bins[1])
    # values of the histogram where the points are
    hhsub = hh[posx[ind] - 1, posy[ind] - 1]
    xdat1 = xdat[ind][hhsub < thresh]  # low density points
    ydat1 = ydat[ind][hhsub < thresh]
    hh[hh < thresh] = np.nan  # fill the areas with low density by NaNs

    plt.imshow(
        np.flipud(
            hh.T),
        cmap='jet',
        extent=np.array(xyrange).flatten(),
        interpolation='none')
    plt.plot(xdat1, ydat1, '.')
    plt.show()


def Kristen_robust_binarize(binary_labels, base_image):
    print "reached Kristen robust binarize debug"
    plt.figure()
    plt.subplot(121)
    plt.imshow(binary_labels)
    plt.title("Binary labels")
    plt.subplot(122)
    plt.imshow(base_image)
    plt.figure("Base Image")
    labels, _ = ndi.label(binary_labels)
    # dist = ndi.morphology.distance_transform_edt(np.logical_not(labels))
    # segmented_cells_labels = watershed(dist, labels)
    unique_segmented_cells_labels = np.unique(labels)
    average_intensity_list = []
    total_intensity_list = []
    pixel_list = []
    for cell_label in unique_segmented_cells_labels:
        my_mask = labels == cell_label
        average_apply_mask = np.mean(base_image[my_mask])

        average_intensity_list.append(average_apply_mask)
        intensity = np.sum(base_image[my_mask])
        total_intensity_list.append(intensity)
        binary_pad = np.zeros_like(base_image)
        binary_pad[my_mask] = 1
        pixel = np.sum(binary_pad[my_mask])
        pixel_list.append(pixel)

    print
    print
    print "avg intensity",average_intensity_list
    print 'total intensity', total_intensity_list
    print 'pixels', pixel_list
    plt.subplot(131)
    plt.plot(average_intensity_list)
    plt.title('average intensity')
    plt.subplot(132)
    plt.plot(total_intensity_list)
    plt.title('total intensity')
    plt.subplot(133)
    plt.plot(pixel_list)
    plt.title('pixel values')
    # plt.show


def in_contact_debug(ic1, ic2):
    plt.figure(figsize=(20.0, 15.0))
    plt.suptitle('In Contact')

    main_ax = plt.subplot(121)
    plt.title('mCherry')
    plt.imshow(ic1, interpolation='nearest', cmap=plt.cm.spectral)

    plt.subplot(122, sharex=main_ax, sharey=main_ax)
    plt.title('GFP')
    plt.imshow(ic2, interpolation='nearest', cmap='gray')
    plt.show()


def label_and_correct_debug(segmented):
    plt.figure(figsize=(20.0, 15.0))
    plt.title('Label and Correct')
    plt.imshow(segmented)


def random_walker_debug(p1, p2):
    plt.figure()
    plt.subplot(211)
    plt.plot(p1)
    plt.subplot(212)
    plt.plot(p2)
    plt.show()


def sum_proj_debug(image):
    plt.figure(figsize=(20.0, 15.0))
    plt.title('Sum Projection')
    plt.imshow(image)
    plt.show()


def DAPI_debug(stabilized, smoothed):
    plt.figure(figsize=(20.0, 15.0))
    plt.suptitle('DAPI Debug')

    #
    # main_ax = plt.subplot(221)
    # plt.title('stabilized')
    # stabilized_2 = np.max(stabilized, axis=0)
    # plt.imshow(stabilized_2, interpolation='nearest', cmap='gray')
    #
    # plt.subplot(222, sharex=main_ax, sharey=main_ax)
    # plt.title('smoothed')
    # smoothed_2 = np.max(smoothed, axis=0)
    # plt.imshow(smoothed_2, interpolation='nearest', cmap='gray')
    #
    # plt.show()
    #


def nuclei_debug(binarized, segmented):
    plt.figure(figsize=(20.0, 15.0))
    plt.suptitle('DAPI Debug')
    print binarized

    main_ax = plt.subplot(221)
    plt.title('binzarized nuclei')
    binarized_2 = np.max(binarized, axis=0)
    plt.imshow(binarized_2, interpolation='nearest', cmap='gray')

    plt.subplot(222, sharex=main_ax, sharey=main_ax)
    plt.title('segmented nuclei')
    segmented_2 = np.max(segmented, axis=0)
    plt.imshow(segmented, interpolation='nearest', cmap='gray')

    plt.show()
