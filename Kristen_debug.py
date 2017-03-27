from matplotlib import pyplot as plt
import numpy as np


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

# problem: you're trying to graph a generator, not a 3-D numpy array