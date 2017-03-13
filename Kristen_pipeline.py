import kristen_traversal as uf
import core_functions as cf
from matplotlib import pyplot as plt
import render as rdr
from csv import writer as csv_writer

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

stabilized_GFP = cf.gamma_stabilize(smoothed_DAPI, in_channel='GFP', min='min', alpha_clean=.0)
smoothed_GFP = cf.smooth_2d(stabilized_GFP, in_channel='GFP', smoothing_px=.5)

stabilized_mCherry = cf.gamma_stabilize(smoothed_GFP, in_channel='mCherry', min='5p', alpha_clean=.5)
smoothed_mCherry = cf.smooth_2d(stabilized_mCherry, in_channel='p21', smoothing_px=.5)

binarized_nuclei = cf.robust_binarize(smoothed_mCherry,
                                      in_channel='DAPI',
                                      out_channel=['nuclei'],
                                      _dilation=0,
                                      heterogeity_size=5, feature_size=50)
segmented_nuclei = cf.label_and_correct(binarized_nuclei,
                                        in_channel=['nuclei', 'DAPI'],
                                        out_channel='nuclei',
                                        min_px_radius=15, min_intensity=20)


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
