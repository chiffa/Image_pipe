from matplotlib import pyplot as plt
import numpy as np
from debugger_skeleton import Debugger, DebugFrame


class DebugFrame1(DebugFrame):

    _fields = ['name_pattern', 'gfp_collector', 'gfp_clustering_markers', 'labels',
               'segmented_cells_labels', 'std_err', 'average_gfp_in_cell_mask',
               'cells_average_gfp_list', 'non_dying_predicted', 'non_dying_cells_mask',
               'qualifying_gfp_mask']

    def render(self, show=False):
        plt.figure(figsize=(20.0, 15.0))

        plt.title(self.name_pattern)

        plt.subplot(241)
        plt.imshow(self.gfp_collector, interpolation='nearest')

        plt.subplot(242)
        plt.imshow(self.gfp_clustering_markers, cmap='hot', interpolation='nearest')

        plt.subplot(243)
        plt.imshow(self.labels, cmap='gray', interpolation='nearest')

        plt.subplot(244)
        plt.imshow(self.segmented_cells_labels, cmap=plt.cm.spectral, interpolation='nearest')

        if self.std_err is not None:

            plt.subplot(245)
            plt.imshow(self.average_gfp_in_cell_mask, cmap='hot', interpolation='nearest')
            plt.colorbar()

            plt.subplot(246)
            plt.plot(self.cells_average_gfp_list, 'ko')
            plt.plot(self.non_dying_predicted, 'r')
            plt.plot(self.non_dying_predicted + self.std_err, 'g')
            plt.plot(self.non_dying_predicted - self.std_err, 'g')

            plt.subplot(247)
            plt.imshow(self.non_dying_cells_mask, cmap='gray', interpolation='nearest')

            plt.subplot(248)
            plt.imshow(self.qualifying_gfp_mask)

        if show:
            plt.show()

        else:
            plt.savefig('verification_bank/%s.png' % self.name_pattern)
            plt.close()


class DebugFrame2(DebugFrame):

    _fields = ['gfp_collector_pre', 'mch_collector_pre',
               'gfp_collector_post', 'mch_collector_post',
               'name_pattern']

    def render(self, show=False):

        plt.figure(figsize=(20.0, 15.0))

        plt.subplot(221)
        plt.imshow(self.gfp_collector_pre, cmap='Greens')

        plt.subplot(222)
        plt.imshow(self.mch_collector_pre, cmap='Reds')

        plt.subplot(223)
        plt.imshow(self.gfp_collector_post, cmap='Greens')

        plt.subplot(224)
        plt.imshow(self.mch_collector_post, cmap='Reds')

        if show:
            plt.show()

        else:
            plt.savefig('verification_bank/core-%s.png' % self.name_pattern)
            plt.close()


class DebugFrame3(DebugFrame):

    _fields = ['name_pattern', 'mch_collector', 'labels', 'mean_width', 'mean_length', 'skeleton'
               'mito_binary_mask', 'segmented_cells', 'numbered_lables', 'paint_area',
               'classification_pad']

    def render(self, show=False):
        plt.figure(figsize=(20.0, 15.0))

        ax1 = plt.subplot(231)
        plt.title(self.name_pattern)
        plt.imshow(self.mch_collector, cmap='Reds')
        plt.colorbar()
        plt.contour(self.labels, [0.5], colors='k')

        plt.subplot(232, sharex=ax1, sharey=ax1)
        plt.title('width ; length - av: %.2f ; %.2f' % (self.mean_width, self.mean_length))
        plt.imshow(self.skeleton, cmap=plt.cm.spectral, interpolation='nearest')
        plt.colorbar()
        plt.contour(self.mito_binary_mask, [0.5], colors='w')

        plt.subplot(233, sharex=ax1, sharey=ax1)
        plt.imshow(self.segmented_cells, cmap=plt.cm.spectral, interpolation='nearest')
        plt.colorbar()
        plt.contour(self.mito_binary_mask, [0.5], colors='w')

        plt.subplot(234, sharex=ax1, sharey=ax1)
        # numbered_skeleton_label ?
        plt.imshow(self.numbered_lables, cmap=plt.cm.spectral, interpolation='nearest')
        plt.colorbar()
        plt.contour(self.mito_binary_mask, [0.5], colors='w')

        plt.subplot(235, sharex=ax1, sharey=ax1)
        plt.imshow(self.paint_area, cmap='hot', interpolation='nearest')
        plt.colorbar()
        plt.contour(self.mito_binary_mask, [0.5], colors='w')

        plt.subplot(236, sharex=ax1, sharey=ax1)
        plt.imshow(self.classification_pad, cmap=plt.cm.spectral, interpolation='nearest')
        plt.colorbar()
        plt.contour(self.mito_binary_mask, [0.5], colors='w')

        if show:
            plt.show()

        else:
            plt.savefig('verification_bank/mitochondria-%s.png' % self.name_pattern)
            plt.close()


class CustomDebugger(Debugger):

    frame1 = DebugFrame1()
    frame2 = DebugFrame2()
    frame3 = DebugFrame3()


# TODO: Redo it in terms of the rendering of the final_namespace namespace:
#   - if item exists, bind it to namespace and render according to the procedure
#   - if the item doesn't exist, leave empty space in the render frame.





if __name__ == "__main__":
    debugger = CustomDebugger()
