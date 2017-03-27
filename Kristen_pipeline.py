import kristen_traversal as uf
import core_functions as cf
from matplotlib import pyplot as plt
import render as rdr
from csv import writer as csv_writer
import Kristen_debug as dbg

# Goal of this pipeline
#     1. Detect the number of cells that were properly stained
#     2. For the successfully stained cells, determine how much GFP is located inside the mitochondria


translator = {'C1':0,
              'C3':1,
              'C4':2}

source = uf.Kristen_traverse("/run/user/1000/gvfs/smb-share:server=10.17.0.219,share=common/Users/kristen/Split GFP quant_Andrei/", matching_map=translator)

named_source = uf.name_channels(source, ['DAPI','GFP', 'mCherry'])


stabilized_DAPI = cf.gamma_stabilize(named_source, in_channel='DAPI', min='min', alpha_clean=.5)
smoothed_DAPI = cf.smooth_2d(stabilized_DAPI, in_channel='DAPI', smoothing_px=.5)
# dbg.DAPI_debug(stabilized_DAPI, smoothed_DAPI)

stabilized_GFP = cf.gamma_stabilize(smoothed_DAPI, in_channel='GFP', min='min', alpha_clean=.0)
smoothed_GFP = cf.smooth_2d(stabilized_GFP, in_channel='GFP', smoothing_px=.5)

stabilized_mCherry = cf.gamma_stabilize(smoothed_GFP, in_channel='mCherry', min='5p', alpha_clean=.5)
smoothed_mCherry = cf.smooth_2d(stabilized_mCherry, in_channel='p21', smoothing_px=.5)

print "stabilization complete"

binarized_nuclei = cf.robust_binarize(smoothed_mCherry,
                                      in_channel='DAPI',
                                      out_channel=['nuclei'],
                                      _dilation=0,
                                      heterogeity_size=5, feature_size=50)
print "binarization of nucleus complete"
segmented_nuclei = cf.label_and_correct(binarized_nuclei,
                                        in_channel=['nuclei', 'DAPI'],
                                        out_channel='nuclei',
                                        min_px_radius=15, min_intensity=20)

print "segmentation of nuclei complete"
# need segmentation of DAPI?
# Segmentation of GFP
GFP_aq =  cf.label_based_aq(segmented_nuclei,
                           in_channel=['nuclei', 'GFP'],
                           out_channel=['av_GFP','av_GFP_pad'])
GFP_o_n = cf.exclude_region(GFP_aq,
                            in_channel=['nuclei', 'GFP'],
                            out_channel='GFP_o_n')

vor_seg = cf.voronoi_segment_labels(GFP_o_n,
                                    in_channel='nuclei',
                                    out_channel='vor_segment')

GFP_o_n_segmented = cf.robust_binarize(vor_seg,
                                       in_channel='GFP_o_n',
                                       out_channel='GFP_o_n_seg',
                                       heterogeity_size=10, feature_size=100)

GFP_seg_contacted = cf.in_contact(GFP_o_n_segmented,
                                  in_channel=['GFP_o_n_seg', 'nuclei'],
                                  out_channel=['GFP_o_n_seg', '_'])



GFP_o_n_filtered = cf.filter_labels(GFP_seg_contacted,
                                    in_channel=['vor_segment', 'GFP_o_n_seg'],
                                    out_channel=['extra_nuclear_GFP'],
                                    min_feature_size=30)

GFP_en_eq = cf.label_based_aq(GFP_o_n_filtered,
                              in_channel=['extra_nuclear_GFP','GFP_o_n'],
                              out_channel=['av_en_GFP', 'av_en_GFP_pad'])

# Segmentation of mCherry
mCherry_aq = cf.label_based_aq(GFP_en_eq,
                           in_channel=['nuclei', 'mCherry'],
                           out_channel=['nuc_mCherry', 'nuc_mCherry_pad'])

mCherry_o_n = cf.exclude_region(mCherry_aq,
                            in_channel=['nuclei', 'mCherry'],
                            out_channel='mCherry_o_n')

mCherry_o_n_segmented = cf.robust_binarize(mCherry_o_n,
                                       in_channel='mCherry_o_n',
                                       out_channel='mCherry_o_n_seg',
                                       heterogeity_size=10, feature_size=100)

mCherry_seg_contacted = cf.in_contact(mCherry_o_n_segmented,
                                  in_channel=['mCherry_o_n_seg', 'nuclei'],
                                  out_channel=['mCherry_o_n_seg', '_'])

mCherry_o_n_filtered = cf.filter_labels(mCherry_seg_contacted,
                                    in_channel=['vor_segment', 'mCherry_o_n_seg'],
                                    out_channel=['extra_nuclear_mCherry'],
                                    min_feature_size=30)
mCherry_en_eq = cf.label_based_aq(mCherry_o_n_filtered,
                              in_channel=['extra_nuclear_mCherry', 'mCherry_o_n'],
                              out_channel=['av_en_mCherry', 'av_en_mCherry_pad'])

# Derivation from Linhao's Pipeline
projected_GFP = cf.sum_projection(stabilized_mCherry,
                                  in_channel='GFP',
                                  out_channel='projected_GFP')

projected_mCherry = cf.max_projection(projected_GFP,
                                  in_channel='mCherry',
                                  out_channel='projected_mCh')

# which variable from Linhao's pipeline represents the quantification of GFP Included in the volume encompassed by mCherry???
# Next steps
#     once pipeline is perceived to be finished, choose one image (i guess 3 superimposed images-for DAPI, GFP, mCherry) and run pipeline on
#     add debug renders to see what exactly the pipeline is doing
#     repeat this with other images until pipeline functions properly








#
# rdr.Kristen_render(mCherry_en_eq, in_channel=['name pattern', 'DAPI', 'GFP', 'mCherry',
#                                                'nuclei', 'vor_segment',
#                                                'extra_nuclear_GFP', 'av_GFP_pad', 'av_en_GFP_pad',
#                                                'extra_nuclear_mCherry', 'nuc_mCherry_pad', 'av_en_mCherry_pad'], out_channel = '_', save=False)
#
# running_render = rdr.Kristen_render(mCherry_en_eq,
#                                    in_channel=['name pattern', 'DAPI', 'GFP', 'mCherry',
#                                                'nuclei', 'vor_segment',
#                                                'extra_nuclear_GFP', 'av_GFP_pad', 'av_en_GFP_pad',
#                                                'extra_nuclear_mCherry', 'nuc_mCherry_pad', 'av_en_mCherry_pad'],
#                                    out_channel='_',
#                                    save=True)
#
#
# # why is the warning above showing up? (same exact formatting as in Akshay's pipeline)
# Kristen_summary = rdr.Kristen_summarize(running_render, in_channel=['name pattern', 'group id', 'av_GFP', 'av_en_GFP',
#                                            'nuc_mCherry', 'av_en_mCherry'],
#                                out_channel='_',
#                                output='kristen_analysis_results.csv' )
#
# with open('Kristen_analysis_results.csv', 'wb') as output_file:
#     writer = csv_writer(output_file)
#     writer.writerow(['file', 'group id', 'cell no', 'nuclear GFP',
#                      'cellular GFP', 'nuclear mCherry', 'cellular mCherry'])
#
# for i, elt in enumerate(Kristen_summary):
#     print 'operation %s analyzed group %s - image %s' % (i, elt['group id'], elt['name pattern'])