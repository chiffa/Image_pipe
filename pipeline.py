import debugger_skeleton
import filters as uf
import core_functions as cf
from matplotlib import pyplot as plt
import render as rdr

# different debugger injection can be performed by doing a
# import module_of_interest
# module_of_interest.debugger = new_debugger

translator = {'w1488': 0,
              'w2561': 1
              }

source = uf.traverse_and_match("L:\\Users\\linghao\\Data for quantification\\Yeast\\NEW data for analysis06152016",
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
                                   out_channel=['average_GFP', 'average_GFP_pad'])

GFP_upper_outlier_cells = cf.detect_upper_outliers(average_GFP,
                                                   in_channel='average_GFP',
                                                   out_channel=['upper_outliers', 'pred_gpf_av', 'gfp_std'])

GFP_outliers = cf.paint_mask(GFP_upper_outlier_cells,
                                 in_channel=['cell_labels', 'upper_outliers'],
                                 out_channel='GFP_outliers')

no_outliers = cf.clear_based_on_2d_mask(GFP_outliers,
                                        in_channel=['GFP', 'GFP_outliers'],
                                        out_channel='GFP')

rendered = rdr.gfp_render(no_outliers,
                          in_channel=['name pattern', 'projected_GFP', 'qualifying_GFP',
                                      'cell_labels', 'average_GFP_pad', 'average_GFP',
                                      'pred_gpf_av', 'gfp_std', 'upper_outliers', 'GFP_outliers'],
                          out_channel='_')

cleared = cf.clear_based_on_2d_mask(rendered,
                                        in_channel=['mCherry', 'GFP_outliers'],
                                        out_channel='mCherry')

per_cell_split = cf.splitter(cleared, 'per_cell',
                             sources=['GFP', 'mCherry'],
                             mask='cell_labels')

per_cell_mito = cf.for_each(per_cell_split, cf.binarize_2d, 'per_cell',
                            in_channel='mCherry',
                            out_channel='mito_binary')

skeletonized = cf.for_each(per_cell_mito, cf.agreeing_skeletons, 'per_cell',
                           in_channel=['mCherry', 'mito_binary'],
                           out_channel='mCh_skeleton')

classified = cf.for_each(skeletonized, cf.classify_fragmentation_for_mitochondria, 'per_cell',
                         in_channel=['mito_binary', 'mCh_skeleton'],
                         out_channel=['final_classification', 'classification_mask',
                                      'radius_mask', 'support_mask'])


# TODO: tiling skeletons, mitochondria and classifications

for payload in classified:
    for key, value in payload['per_cell'].iteritems():
        print key
        if not (isinstance(key, basestring) and key[0] == '_'):
            print value.keys()
        else:
            print 'Not a namespace'
    plt.imshow(payload['GFP'][6, :, :])
    plt.show()
