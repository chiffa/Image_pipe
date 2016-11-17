import debugger_skeleton
import filters as uf
import core_functions as cf
from matplotlib import pyplot as plt
import render as rdr

source = uf.Akshay_traverse("L:\\Common\\AKN\\IF test images nuc vs cyto\\tiff")
named_source = uf.name_channels(source, ['DAPI', 'p53', 'p21'])

stablilized_1 = cf.gamma_stabilize(named_source, in_channel='DAPI', min='min', alpha_clean=.5)
smoothed_1 = cf.smooth_2d(stablilized_1, in_channel='DAPI', smoothing_px=.5)  # corresponds to a pixel of 1500 nm

stablilized_2 = cf.gamma_stabilize(smoothed_1, in_channel='p53', min='min', alpha_clean=2)
smoothed_2 = cf.smooth_2d(stablilized_2, in_channel='p53', smoothing_px=.5)

stablilized_3 = cf.gamma_stabilize(smoothed_2, in_channel='p21', min='5p', alpha_clean=2)
smoothed_3 = cf.smooth_2d(stablilized_3, in_channel='p21', smoothing_px=.5)

binarized_nuclei = cf.robust_binarize(smoothed_3,
                                      in_channel='DAPI',
                                      out_channel=['nuclei'],
                                      _dilation=0)

segmented_nuclei = cf.label_and_correct(binarized_nuclei,
                                        in_channel=['nuclei', 'DAPI'],
                                        out_channel='nuclei',
                                        min_px_radius=15, min_intensity=5)

p53_aq = cf.label_based_aq(segmented_nuclei,
                           in_channel=['nuclei', 'p53'],
                           out_channel=['av_p53', 'av_p53_pad'])

p53_o_n = cf.exclude_region(p53_aq,
                            in_channel=['nuclei', 'p53'],
                            out_channel='p53_o_n')

# Current problem - what if the P53 level is so low, we can't segment anything?
# => in this case the answer is to use voronoi segmentation and use the levels reported from it
# as p53 outside levels.
p53_o_n_segmented = cf.robust_binarize(p53_o_n,
                                       in_channel='p53_o_n',
                                       out_channel='p53_o_n_seg',
                                       heterogeity_size=10, feature_size=70)

running_render = rdr.akshay_render(p53_o_n_segmented,
                                   in_channel=['name pattern', 'DAPI', 'p53', 'p21',
                                               'nuclei', 'av_p53_pad', 'p53_o_n', 'p53_o_n_seg'],
                                   out_channel='_')

for elt in running_render:
    pass
