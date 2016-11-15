import os
import numpy as np
from matplotlib import pyplot as plt
from core_functions import generator_wrapper


@generator_wrapper(in_dims=(None, 2, 2, 2, 2, 1, 1, None, None, 2), out_dims=(None,))
def gfp_render(name_pattern, proj_gfp, qual_gfp,
               cell_labels, average_gfp_pad, average_gfp_list,
               predicted_average, std_err, upper_outliers, gfp_outliers,
               save=False, directory_to_save_to='verification'):
    plt.figure(figsize=(20, 15))

    plt.title(name_pattern)

    plt.subplot(241)
    plt.imshow(proj_gfp, interpolation='nearest')

    # plt.subplot(242)
    # plt.imshow(qual_gfp, cmap='hot', interpolation='nearest')

    plt.subplot(243)
    plt.imshow(qual_gfp, cmap='gray', interpolation='nearest')

    plt.subplot(244)
    plt.imshow(cell_labels, cmap=plt.cm.spectral, interpolation='nearest')

    plt.subplot(245)
    plt.imshow(average_gfp_pad, cmap='hot', interpolation='nearest')
    plt.colorbar()

    plt.subplot(246)
    plt.plot(average_gfp_list, 'ko')
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
        name_puck = directory_to_save_to+'/'+name_pattern+'.png'
        plt.savefig(name_puck)
        plt.close()
