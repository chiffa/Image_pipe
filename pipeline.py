import debugger_skeleton
import filters as uf
import core_functions as cf
from matplotlib import pyplot as plt

# different debugger injection can be performed by doing a
# import module_of_interest
# module_of_interest.debugger = new_debugger

translator = {'w1488': 0,
              'w2561': 1
              }

source = uf.traverse_and_match("L:\\Users\\linghao\\Data for quantification\\Yeast\\NEW data for analysis",
                               matching_map=translator)

# that architecture actually removes the need for the debug line

named_source = uf.name_channels(source, ['GFP', 'mCherry'])
stabilized_GFP = cf.gamma_stabilize(named_source, in_channel='GFP')
smoothed_GFP = cf.smooth(stabilized_GFP, in_channel='GFP')
projected_GFP = cf.sum_projection(smoothed_GFP,
                                  in_channel='GFP',
                                  out_channel='projected_GFP')
segmented_GFP = cf.segment_out_cells(projected_GFP,
                                     in_channel='projected_GFP',
                                     out_channel='cell_labels')
qualifying_GFP = cf.qualifying_gfp(segmented_GFP,
                                   in_channel='projected_GFP',
                                   out_channel='qualifying_GFP')
average_GFP = cf.aq_gfp_per_region(qualifying_GFP,
                                   in_channel=['cell_labels', 'projected_GFP', 'qualifying_GFP'],
                                   out_channel='average_GFP', )
GFP_upper_outlier_cells = cf.detect_upper_outliers(average_GFP,
                                                   in_channel='average_GFP',
                                                   out_channel='upper_outliers',
                                                   log_channel='outlier_log')
non_GFP_outliers = cf.paint_mask(GFP_upper_outlier_cells,
                                 in_channel=['cell_labels', 'upper_outliers'],
                                 out_channel='non_GFP_outliers')


for payload in non_GFP_outliers:
    plt.imshow(payload['cell_labels'])
    plt.show()
