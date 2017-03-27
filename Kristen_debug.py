from matplotlib import pyplot as plt
import numpy as np


def DAPI_debug(stabilized, smoothed):
    plt.figure(figsize=(20.0, 15.0))
    plt.suptitle('DAPI Debug')

    main_ax = plt.subplot(221)
    plt.title('stabilized')
    plt.imshow(stabilized, interpolation='nearest', cmap=plt.cm.spectral)

    plt.subplot(222, sharex=main_ax, sharey=main_ax)
    plt.title('smoothed')
    plt.imshow(smoothed, interpolation='nearest', cmap='gray')

    plt.show()

