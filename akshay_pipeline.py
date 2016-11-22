import debugger_skeleton
import filters as uf
import core_functions as cf
from matplotlib import pyplot as plt
import render as rdr
from csv import writer as csv_writer

# source = uf.Akshay_traverse("L:\\Common\\AKN\\IF test images nuc vs cyto\\tiff")
source = uf.Akshay_traverse("L:\\Common\\AKN\\IF test images nuc vs cyto\\11-21-16 RPE IF untreat mps nut images")
named_source = uf.name_channels(source, ['DAPI', 'p53', 'p21'])

stablilized_1 = cf.gamma_stabilize(named_source, in_channel='DAPI', min='min', alpha_clean=.5)
smoothed_1 = cf.smooth_2d(stablilized_1, in_channel='DAPI', smoothing_px=.5)  # corresponds to a pixel of 1500 nm

stablilized_2 = cf.gamma_stabilize(smoothed_1, in_channel='p53', min='min', alpha_clean=.0)
smoothed_2 = cf.smooth_2d(stablilized_2, in_channel='p53', smoothing_px=.5)

stablilized_3 = cf.gamma_stabilize(smoothed_2, in_channel='p21', min='5p', alpha_clean=.5)
smoothed_3 = cf.smooth_2d(stablilized_3, in_channel='p21', smoothing_px=.5)

binarized_nuclei = cf.robust_binarize(smoothed_3,
                                      in_channel='DAPI',
                                      out_channel=['nuclei'],
                                      _dilation=0,
                                      heterogeity_size=5, feature_size=50)

segmented_nuclei = cf.label_and_correct(binarized_nuclei,
                                        in_channel=['nuclei', 'DAPI'],
                                        out_channel='nuclei',
                                        min_px_radius=15, min_intensity=20)

# p53 segmentation

p53_aq = cf.label_based_aq(segmented_nuclei,
                           in_channel=['nuclei', 'p53'],
                           out_channel=['av_p53', 'av_p53_pad'])

p53_o_n = cf.exclude_region(p53_aq,
                            in_channel=['nuclei', 'p53'],
                            out_channel='p53_o_n')

vor_seg = cf.voronoi_segment_labels(p53_o_n,
                                    in_channel='nuclei',
                                    out_channel='vor_segment')

p53_o_n_segmented = cf.robust_binarize(vor_seg,
                                       in_channel='p53_o_n',
                                       out_channel='p53_o_n_seg',
                                       heterogeity_size=10, feature_size=100)

p53_seg_contacted = cf.in_contact(p53_o_n_segmented,
                                  in_channel=['p53_o_n_seg', 'nuclei'],
                                  out_channel=['p53_o_n_seg', '_'])

# TODO: test the elimination

p53_o_n_filtered = cf.filter_labels(p53_seg_contacted,
                                    in_channel=['vor_segment', 'p53_o_n_seg'],
                                    out_channel=['extra_nuclear_p53'],
                                    min_feature_size=30)

p53_en_eq = cf.label_based_aq(p53_o_n_filtered,
                              in_channel=['extra_nuclear_p53', 'p53_o_n'],
                              out_channel=['av_en_p53', 'av_en_p53_pad'])

# p21 now

p21_aq = cf.label_based_aq(p53_en_eq,
                           in_channel=['nuclei', 'p21'],
                           out_channel=['nuc_p21', 'nuc_p21_pad'])

p21_o_n = cf.exclude_region(p21_aq,
                            in_channel=['nuclei', 'p21'],
                            out_channel='p21_o_n')

p21_o_n_segmented = cf.robust_binarize(p21_o_n,
                                       in_channel='p21_o_n',
                                       out_channel='p21_o_n_seg',
                                       heterogeity_size=10, feature_size=100)

p21_seg_contacted = cf.in_contact(p21_o_n_segmented,
                                  in_channel=['p21_o_n_seg', 'nuclei'],
                                  out_channel=['p21_o_n_seg', '_'])

p21_o_n_filtered = cf.filter_labels(p21_seg_contacted,
                                    in_channel=['vor_segment', 'p21_o_n_seg'],
                                    out_channel=['extra_nuclear_p21'],
                                    min_feature_size=30)

p21_en_eq = cf.label_based_aq(p21_o_n_filtered,
                              in_channel=['extra_nuclear_p21', 'p21_o_n'],
                              out_channel=['av_en_p21', 'av_en_p21_pad'])

running_render = rdr.akshay_render(p21_en_eq,
                                   in_channel=['name pattern', 'DAPI', 'p53', 'p21',
                                               'nuclei', 'vor_segment',
                                               'extra_nuclear_p53', 'av_p53_pad', 'av_en_p53_pad',
                                               'extra_nuclear_p21', 'nuc_p21_pad', 'av_en_p21_pad'],
                                   out_channel='_',
                                   save=False)

summary = rdr.akshay_summarize(running_render,
                               in_channel=['name pattern', 'group id', 'av_p53', 'av_en_p53',
                                           'nuc_p21', 'av_en_p21'],
                               out_channel='_',
                               output='analys_results.csv')

with open('analys_results.csv', 'wb') as output_file:
        writer = csv_writer(output_file)
        writer.writerow(['file', 'group id', 'cell no', 'nuclear p53',
                         'cellular p53', 'nuclear p21', 'cellular p21'])

for elt in summary:
    print 'analyzed %s' % elt['name pattern']

# TODO: test multi-scale median as a way of detecting the DAPI nuclei as opposed to local OTSU
