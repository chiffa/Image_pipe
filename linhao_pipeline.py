import traversals as uf
import core_functions as cf
from matplotlib import pyplot as plt
import render as rdr
from csv import writer as csv_writer
from time import time, strftime

# different debugger injection can be performed by doing a
# import module_of_interest
# module_of_interest.debugger = new_debugger

translator = {'w1488': 0,
              'w2561': 1
              }


# source = uf.Linhao_traverse("L:\\Users\\jerry\\Image\\ForAndrei\\Ry233282285",
#                             matching_map=translator)

source = uf.Linhao_traverse("/run/user/1000/gvfs/smb-share:server=10.17.0.219,share=common/Users/jerry/Image/ForAndrei/SSAmutant/",
                            matching_map=translator)

# that architecture actually removes the need for the debug line
# that can be incorporated into the traversal

named_source = uf.name_channels(source, ['GFP', 'mCherry'])

stabilized_GFP = cf.gamma_stabilize(named_source, in_channel='GFP')
smoothed_GFP = cf.smooth(stabilized_GFP, in_channel='GFP')

stabilized_mCh = cf.gamma_stabilize(smoothed_GFP, in_channel='mCherry')

projected_GFP = cf.sum_projection(stabilized_mCh,
                                  in_channel='GFP',
                                  out_channel='projected_GFP')

projected_mCh = cf.max_projection(projected_GFP,
                                  in_channel='mCherry',
                                  out_channel='projected_mCh')


binarized_GFP = cf.robust_binarize(projected_mCh,
                                   in_channel='projected_mCh',
                                   out_channel='cell_tags')

segmented_GFP = cf.improved_watershed(binarized_GFP,
                                      in_channel='cell_tags',
                                      out_channel='pre_cell_labels')

qualifying_GFP = cf.qualifying_gfp(segmented_GFP,
                                   in_channel='projected_GFP',
                                   out_channel='qualifying_GFP')

average_GFP = cf.aq_gfp_per_region(qualifying_GFP,
                                   in_channel=['pre_cell_labels', 'projected_GFP', 'qualifying_GFP'],
                                   out_channel=['average_GFP', 'average_GFP_pad'])

GFP_upper_outlier_cells = cf.detect_upper_outliers(average_GFP,
                                                   in_channel='average_GFP',
                                                   out_channel=['non_outliers', 'pred_gpf_av',
                                                                'gfp_std'])

GFP_outliers = cf.paint_mask(GFP_upper_outlier_cells,
                             in_channel=['pre_cell_labels', 'non_outliers'],
                             out_channel='kept_cells')

GFP_filtered = cf.mask_filter_2d(GFP_outliers,
                                 in_channel=['pre_cell_labels', 'kept_cells'],
                                 out_channel='cell_labels')


per_cell_split = cf.splitter(GFP_filtered, 'per_cell',
                             sources=['mCherry', 'GFP',
                                      'projected_mCh', 'projected_GFP'],
                             mask='cell_labels')

per_cell_mito = cf.for_each(per_cell_split, cf.binarize_2d, 'per_cell', cutoff_type='otsu',
                            in_channel='projected_mCh',
                            out_channel='mito_binary')

segmented_mito = cf.for_each(per_cell_mito, cf.label_and_correct, 'per_cell',
                             in_channel=['mito_binary', 'projected_mCh'],
                             out_channel='mito_labels')

mito_3d_from_2d_mask = cf.for_each(segmented_mito, cf._3d_mask_from_2d_mask, 'per_cell',
                                    in_channel=['mCherry', 'mito_labels'],
                                    out_channel='mito_labels_3d')

# problem - mqvi does not seem to be working on an indifidual basis

GFP_AEQVI = cf.for_each(mito_3d_from_2d_mask, cf.volume_aqvi, 'per_cell',
                        in_channel=['GFP', 'mito_labels_3d'],
                        out_channel='gfp_mqvi')

MCH_AEQVI = cf.for_each(GFP_AEQVI, cf.volume_aqvi, 'per_cell',
                        in_channel=['mCherry', 'mito_labels_3d'],
                        out_channel='mch_mqvi')

skeletonized = cf.for_each(MCH_AEQVI, cf.agreeing_skeletons, 'per_cell',
                           in_channel=['projected_mCh', 'mito_binary'],
                           out_channel='mCh_skeleton')



classified = cf.for_each(skeletonized, cf.classify_fragmentation_for_mitochondria, 'per_cell',
                         in_channel=['mito_labels', 'mCh_skeleton'],
                         out_channel=['final_classification', 'classification_mask',
                                      'radius_mask', 'support_mask'])

mito_tiled = cf.tile_from_mask(classified, 'per_cell', 'mito_binary')

skeleton_tiled = cf.tile_from_mask(mito_tiled, 'per_cell', 'mCh_skeleton')

classification_tiled = cf.tile_from_mask(skeleton_tiled, 'per_cell', 'classification_mask')

cell_class_tiled = cf.paint_from_mask(classification_tiled, 'per_cell', 'final_classification')

rad_mask_tiled = cf.tile_from_mask(cell_class_tiled, 'per_cell', 'radius_mask')

supp_mask_tiled = cf.tile_from_mask(rad_mask_tiled, 'per_cell', 'support_mask')

gfp_mqvi_tiled = cf.paint_from_mask(supp_mask_tiled, 'per_cell', 'gfp_mqvi')

mch_mqvi_tiled = cf.paint_from_mask(gfp_mqvi_tiled, 'per_cell', 'mch_mqvi')

gfp_rendered = rdr.linhao_gfp_render(mch_mqvi_tiled,
                                     in_channel=['name pattern',
                                                 'projected_GFP', 'qualifying_GFP',
                                                 'pre_cell_labels',
                                                 'average_GFP_pad', 'average_GFP',
                                                 'pred_gpf_av', 'gfp_std', 'non_outliers',
                                                 'cell_labels', 'projected_mCh'],
                                     out_channel='_',
                                     save=True)

mqvi_render = rdr.linhao_mqvi_render(gfp_rendered,
                                    in_channel=['name pattern', 'mito_binary', 'cell_labels',
                                                'projected_GFP', 'projected_mCh',
                                                'gfp_mqvi', 'mch_mqvi'],
                                    out_channel='_',
                                    save=True)

mch_render = rdr.linhao_mch_render(mqvi_render,
                                        in_channel=['name pattern', 'projected_mCh', 'mito_binary',
                                             'mCh_skeleton', 'classification_mask', 'final_classification',
                                             'cell_labels', 'radius_mask', 'support_mask'],
                                        out_channel='_',
                                        save=True)

per_cell_render = rdr.linhao_summarize(mch_render, output='linhao_analys_results.csv')

cell_count = rdr.linhao_secondary_summarize(per_cell_render, output='linhao_raw counts.csv')

with open('linhao_analys_results.csv', 'wb') as output_file:
        writer = csv_writer(output_file)
        writer.writerow(['file', 'time in curve', 'date', 'cell type',
                         'cell no', 'gfp_mqvi', 'mch_mqvi', 'mito fragmentation'])


with open('linhao_raw counts.csv', 'wb') as output_file:
        writer = csv_writer(output_file)
        writer.writerow(['file', 'time in curve', 'date', 'cell type', 'cells detected', 'cells analyzed'])

prev_time = time()

for primary_namespace in cell_count:

    #my addition
    print
    print '****************************************************'
    print "primary namespace", primary_namespace


    print '%s - analyzed %s in %s' % (strftime('%X %x'), primary_namespace['name pattern'], time() - prev_time)
    prev_time = time()
