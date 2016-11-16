from matplotlib import pyplot as plt


def robust_binarize_debug(base, smooth, median, otsu, labels, binary):
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
    plt.title('median - otsu')
    plt.imshow(median-otsu, interpolation='nearest')
    plt.contour(binary, [0.5], colors='k')


    plt.subplot(247, sharex=main_ax, sharey=main_ax)
    plt.title('labels')
    plt.imshow(labels, interpolation='nearest')
    plt.contour(binary, [0.5], colors='k')

    plt.subplot(248, sharex=main_ax, sharey=main_ax)
    plt.title('binary')
    plt.imshow(binary, interpolation='nearest')

    plt.show()