from matplotlib import pyplot as plt


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

    plt.show()


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

    plt.show()


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
    plt.imshow(result, interpolation='nearest', cmap=plt.cm.spectral, vmin=0),

    # plt.subplot(224, sharex=main_ax, sharey=main_ax)
    # plt.title('watershed')
    # plt.imshow(watershed, interpolation='nearest', cmap='hot')

    plt.show()


def weight_sum_zero_debug(label_mask, img_codename):
    plt.figure(figsize=(20.0, 15.0))
    plt.suptitle('weight sum zero debugger')

    main_ax = plt.subplot(121)
    plt.title('label_mask')
    plt.imshow(label_mask, interpolation='nearest', cmap=plt.cm.spectral)

    plt.subplot(122, sharex=main_ax, sharey=main_ax)
    plt.title('image code')
    plt.imshow(img_codename, interpolation='nearest', cmap='gray')

    plt.show()



