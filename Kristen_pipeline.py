import kristen_traversal as uf
import core_functions as cf
from matplotlib import pyplot as plt
import render as rdr
from csv import writer as csv_writer
import debug_renders as dbg
import numpy as np
from time import time

# Goal of this pipeline
#     1. Detect the number of cells that were properly stained
#     2. For the successfully stained cells, determine how much GFP is located inside the mitochondria


translator = {'C1':0,
              'C3':1,
              'C4':2}

source = uf.Kristen_traverse('/run/user/1000/gvfs/smb-share:server=10.17.0.219,share=common/Users/kristen/Split GFP quant_Andrei/20170209/Transfection B', matching_map=translator)
named_source = uf.name_channels(source, ['DAPI','GFP', 'mCherry'])

max_mCherry = cf.max_projection(named_source, in_channel = 'mCherry', out_channel = 'max_mCherry')
stabilized_mCherry = cf.gamma_stabilize(max_mCherry, in_channel = 'max_mCherry', min='min', alpha_clean=.5)
smoothed_mCherry = cf.smooth_2d(stabilized_mCherry, in_channel = 'max_mCherry', smoothing_px=.5)




mCherry_o_n_segmented = cf.robust_binarize(smoothed_mCherry,
                                       in_channel='max_mCherry',
                                       out_channel='max_mCherry_binary',
                                       heterogeity_size=10, feature_size=250)
# original: 10, 100

running_render = rdr.Kristen_render(mCherry_o_n_segmented, in_channel=['name pattern',
                                                               'max_mCherry',
                                                               'max_mCherry_binary'],
                                   out_channel='_',
                                   save=False)


Kristen_summary = rdr.Kristen_summarize_a(running_render, in_channel = ['name pattern', 'group id', 'max_mCherry', 'max_mCherry_binary'],
                                          # in_channel=['name pattern', 'group id', 'av_GFP', 'av_en_GFP',
                                          #  'nuc_mCherry', 'av_en_mCherry'],
                               out_channel='_',
                               output='Kristen_analysis_results_1.csv')















# Derivation from Linhao's Pipeline
# named_source = uf.name_channels(source, ['GFP', 'mCherry'])
# stabilized_GFP = cf.gamma_stabilize(named_source, in_channel = 'GFP')
# smoothed_GFP = cf.smooth(stabilized_GFP, in_channel = 'GFP')
# stabilized_mCh = cf.gamma_stabilize(smoothed_GFP, in_channel='mCherry')
#already stabilized and smoothed (top of akshay's portion)




# COMMENTING OUT LINHAO PORTION FOR NOW
# projected_GFP = cf.sum_projection(Kristen_summary,
#                                   in_channel='GFP',
#                                   out_channel='projected_GFP')
#
# projected_mCh = cf.max_projection(projected_GFP,
#                                   in_channel='mCherry',
#                                   out_channel='projected_mCh')
#
#
# binarized_GFP = cf.robust_binarize(projected_mCh,
#                                    in_channel='projected_mCh',
#                                    out_channel='cell_tags')
#
# segmented_GFP = cf.improved_watershed(binarized_GFP,
#                                       in_channel=['cell_tags', 'projected_mCh'],
#                                       out_channel='pre_cell_labels')
#
# qualifying_GFP = cf.qualifying_gfp(segmented_GFP,
#                                    in_channel='projected_GFP',
#                                    out_channel='qualifying_GFP')
#
# average_GFP = cf.aq_gfp_per_region(qualifying_GFP,
#                                    in_channel=['pre_cell_labels', 'projected_GFP', 'qualifying_GFP'],
#                                    out_channel=['average_GFP', 'average_GFP_pad'])
#
# GFP_upper_outlier_cells = cf.detect_upper_outliers(average_GFP,
#                                                    in_channel='average_GFP',
#                                                    out_channel=['non_outliers', 'pred_gpf_av',
#                                                                 'gfp_std'])
#
# GFP_outliers = cf.paint_mask(GFP_upper_outlier_cells,
#                              in_channel=['pre_cell_labels', 'non_outliers'],
#                              out_channel='kept_cells')
# GFP_filtered = cf.mask_filter_2d(GFP_outliers,
#                                  in_channel=['pre_cell_labels', 'kept_cells'],
#                                  out_channel='cell_labels')
# mito_3d_from_2d_mask = cf.for_each(GFP_filtered, cf._3d_mask_from_2d_mask, 'per_cell',
#                                     in_channel=['mCherry', 'mito_labels'],
#                                     out_channel='mito_labels_3d')
#
# # problem - mqvi does not seem to be working on an individual basis
#
# GFP_AEQVI = cf.for_each(mito_3d_from_2d_mask, cf.volume_aqvi, 'per_cell',
#                         in_channel=['GFP', 'mito_labels_3d'],
#                         out_channel='gfp_mqvi')
#
# MCH_AEQVI = cf.for_each(GFP_AEQVI, cf.volume_aqvi, 'per_cell',
#                         in_channel=['mCherry', 'mito_labels_3d'],
#                         out_channel='mch_mqvi')
#
# skeletonized = cf.for_each(MCH_AEQVI, cf.agreeing_skeletons, 'per_cell',
#                            in_channel=['projected_mCh', 'mito_binary'],
#                            out_channel='mCh_skeleton')
#
# classified = cf.for_each(skeletonized, cf.classify_fragmentation_for_mitochondria, 'per_cell',
#                          in_channel=['mito_labels', 'mCh_skeleton'],
#                          out_channel=['final_classification', 'classification_mask',
#                                       'radius_mask', 'support_mask'])
#
# mito_tiled = cf.tile_from_mask(classified, 'per_cell', 'mito_binary')
#
# skeleton_tiled = cf.tile_from_mask(mito_tiled, 'per_cell', 'mCh_skeleton')
#
# classification_tiled = cf.tile_from_mask(skeleton_tiled, 'per_cell', 'classification_mask')
#
# cell_class_tiled = cf.paint_from_mask(classification_tiled, 'per_cell', 'final_classification')
#
# rad_mask_tiled = cf.tile_from_mask(cell_class_tiled, 'per_cell', 'radius_mask')
#
# supp_mask_tiled = cf.tile_from_mask(rad_mask_tiled, 'per_cell', 'support_mask')
#
# gfp_mqvi_tiled = cf.paint_from_mask(supp_mask_tiled, 'per_cell', 'gfp_mqvi')
#
# mch_mqvi_tiled = cf.paint_from_mask(gfp_mqvi_tiled, 'per_cell', 'mch_mqvi')
# print "gfp render"
# gfp_rendered = rdr.linhao_gfp_render(mch_mqvi_tiled,
#                                      in_channel=['name pattern',
#                                                  'projected_GFP', 'qualifying_GFP',
#                                                  'pre_cell_labels',
#                                                  'average_GFP_pad', 'average_GFP',
#                                                  'pred_gpf_av', 'gfp_std', 'non_outliers',
#                                                  'cell_labels', 'projected_mCh'],
#                                      out_channel='_',
#                                      save=False)
#
# mqvi_render = rdr.linhao_mqvi_render(gfp_rendered,
#                                     in_channel=['name pattern', 'mito_binary', 'cell_labels',
#                                                 'projected_GFP', 'projected_mCh',
#                                                 'gfp_mqvi', 'mch_mqvi'],
#                                     out_channel='_',
#                                     save=False)
#
#
#
# mch_render = rdr.linhao_mch_render(mqvi_render,
#                                         in_channel=['name pattern', 'projected_mCh', 'mito_binary',
#                                              'mCh_skeleton', 'classification_mask', 'final_classification',
#                                              'cell_labels', 'radius_mask', 'support_mask'],
#                                         out_channel='_',
#                                         save=False)
#
# per_cell_render = rdr.linhao_summarize(mch_render, output='Kristen_per_cell_render_results.csv')
#
# cell_count = rdr.linhao_secondary_summarize(per_cell_render, output='Kristen_raw counts.csv')
#
# with open('Kristen_analysis_results_1.csv', 'wb') as output_file:
#     writer = csv_writer(output_file, delimiter = '\t')
#     writer.writerow(['file', 'group id', 'cell number', '   nuclear GFP',
#                      '  cellular GFP', '  nuclear mCherry', '   cellular mCherry'])
#
# with open('Kristen_per_cell_render_results.csv', 'wb') as output_file:
#         writer = csv_writer(output_file, delimiter = '\t')
#         writer.writerow(['file', '  time in curve', '   date', '    cell type',
#                          '  cell no', ' gfp_mqvi', '    mch_mqvi', '    mito fragmentation'])
#
#
# with open('Kristen_raw counts.csv', 'wb') as output_file:
#         writer = csv_writer(output_file, delimiter = '\t')
#         writer.writerow(['file', '  time in curve', '   date', '    cell type', '   cells detected', '  cells analyzed'])


# which variable from Linhao's pipeline represents the quantification of GFP Included in the volume encompassed by mCherry???
# Next steps
#     once pipeline is perceived to be finished, choose one image (i guess 3 superimposed images-for DAPI, GFP, mCherry) and run pipeline on
#     add debug renders to see what exactly the pipeline is doing
#     repeat this with other images until pipeline functions properly
# use mqvi render to unwrap the generator so you can plot the array in plt.imshow (also use gfp render?)
# the portion of this pipeline needed for kristen's pipeline is basically everything except the mitochondria segmentation (up to per cell split)
# add one by one and check each time to make sure pipeline is working instead of it crashing at the end which makes tracing the problem very difficult
# akshay pipeline still not working: "Nonetype object is not iterable"


for i, elt in enumerate(Kristen_summary):
    print 'operation %s analyzed group %s - image %s' % (i, elt['group id'], elt['name pattern'])