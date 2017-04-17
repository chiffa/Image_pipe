import traversals as uf
import core_functions as cf
from matplotlib import pyplot as plt
import render as rdr
from csv import writer as csv_writer

# source = uf.Akshay_traverse("L:\\Common\\AKN\\IF test images nuc vs cyto\\tiff")
# source = uf.Akshay_traverse("L:\\Common\\AKN\\IF test images nuc vs cyto\\11-21-16 RPE IF untreat mps nut images")
import wrapped_image_proc_steps

source = uf.Akshay_traverse('/run/user/1000/gvfs/smb-share:server=10.17.0.219,share=common/Common/AKN/IF test images nuc vs cyto/11-21-16 RPE IF untreat mps nut - Analysis run')
print 'source', source
named_source = uf.name_channels(source, ['DAPI', 'p53', 'p21'])

stablilized_1 = wrapped_image_proc_steps.gamma_stabilize(named_source, in_channel='DAPI', min='min', alpha_clean=.5)
smoothed_1 = wrapped_image_proc_steps.smooth_2d(stablilized_1, in_channel='DAPI', smoothing_px=.5)  # corresponds to a pixel of 1500 nm

stablilized_2 = wrapped_image_proc_steps.gamma_stabilize(smoothed_1, in_channel='p53', min='min', alpha_clean=.0)
smoothed_2 = wrapped_image_proc_steps.smooth_2d(stablilized_2, in_channel='p53', smoothing_px=.5)

stablilized_3 = wrapped_image_proc_steps.gamma_stabilize(smoothed_2, in_channel='p21', min='5p', alpha_clean=.5)
smoothed_3 = wrapped_image_proc_steps.smooth_2d(stablilized_3, in_channel='p21', smoothing_px=.5)

binarized_nuclei = wrapped_image_proc_steps.robust_binarize(smoothed_3,
                                                            in_channel='DAPI',
                                                            out_channel=['nuclei'],
                                                            _dilation=0,
                                                            heterogeity_size=5, feature_size=50)

segmented_nuclei = wrapped_image_proc_steps.label_and_correct(binarized_nuclei,
                                                              in_channel=['nuclei', 'DAPI'],
                                                              out_channel='nuclei',
                                                              min_px_radius=15, min_intensity=20)

# p53 segmentation

p53_aq = wrapped_image_proc_steps.label_based_aq(segmented_nuclei,
                                                 in_channel=['nuclei', 'p53'],
                                                 out_channel=['av_p53', 'av_p53_pad'])

p53_o_n = wrapped_image_proc_steps.exclude_region(p53_aq,
                                                  in_channel=['nuclei', 'p53'],
                                                  out_channel='p53_o_n')

vor_seg = wrapped_image_proc_steps.voronoi_segment_labels(p53_o_n,
                                                          in_channel='nuclei',
                                                          out_channel='vor_segment')

p53_o_n_segmented = wrapped_image_proc_steps.robust_binarize(vor_seg,
                                                             in_channel='p53_o_n',
                                                             out_channel='p53_o_n_seg',
                                                             heterogeity_size=10, feature_size=100)

p53_seg_contacted = wrapped_image_proc_steps.in_contact(p53_o_n_segmented,
                                                        in_channel=['p53_o_n_seg', 'nuclei'],
                                                        out_channel=['p53_o_n_seg', '_'])


p53_o_n_filtered = wrapped_image_proc_steps.filter_labels(p53_seg_contacted,
                                                          in_channel=['vor_segment', 'p53_o_n_seg'],
                                                          out_channel=['extra_nuclear_p53'],
                                                          min_feature_size=30)

p53_en_eq = wrapped_image_proc_steps.label_based_aq(p53_o_n_filtered,
                                                    in_channel=['extra_nuclear_p53', 'p53_o_n'],
                                                    out_channel=['av_en_p53', 'av_en_p53_pad'])

# p21 now

p21_aq = wrapped_image_proc_steps.label_based_aq(p53_en_eq,
                                                 in_channel=['nuclei', 'p21'],
                                                 out_channel=['nuc_p21', 'nuc_p21_pad'])

p21_o_n = wrapped_image_proc_steps.exclude_region(p21_aq,
                                                  in_channel=['nuclei', 'p21'],
                                                  out_channel='p21_o_n')

p21_o_n_segmented = wrapped_image_proc_steps.robust_binarize(p21_o_n,
                                                             in_channel='p21_o_n',
                                                             out_channel='p21_o_n_seg',
                                                             heterogeity_size=10, feature_size=100)

p21_seg_contacted = wrapped_image_proc_steps.in_contact(p21_o_n_segmented,
                                                        in_channel=['p21_o_n_seg', 'nuclei'],
                                                        out_channel=['p21_o_n_seg', '_'])

p21_o_n_filtered = wrapped_image_proc_steps.filter_labels(p21_seg_contacted,
                                                          in_channel=['vor_segment', 'p21_o_n_seg'],
                                                          out_channel=['extra_nuclear_p21'],
                                                          min_feature_size=30)

p21_en_eq = wrapped_image_proc_steps.label_based_aq(p21_o_n_filtered,
                                                    in_channel=['extra_nuclear_p21', 'p21_o_n'],
                                                    out_channel=['av_en_p21', 'av_en_p21_pad'])

running_render = rdr.akshay_render(p21_en_eq,
                                   in_channel=['name pattern', 'DAPI', 'p53', 'p21',
                                               'nuclei', 'vor_segment',
                                               'extra_nuclear_p53', 'av_p53_pad', 'av_en_p53_pad',
                                               'extra_nuclear_p21', 'nuc_p21_pad', 'av_en_p21_pad'],
                                   out_channel='_',
                                   save=False)

# after checking this pipeline, replace False with True

summary = rdr.akshay_summarize(running_render,
                               in_channel=['name pattern', 'group id', 'av_p53', 'av_en_p53',
                                           'nuc_p21', 'av_en_p21'],
                               out_channel='_',
                               output='akshay_analysis_results.csv')

with open('akshay_analysis_results.csv', 'wb') as output_file:
        writer = csv_writer(output_file)
        writer.writerow(['file', 'group id', 'cell no', 'nuclear p53',
                         'cellular p53', 'nuclear p21', 'cellular p21'])

for i, elt in enumerate(summary):
    print 'operation %s analyzed group %s - image %s' % (i, elt['group id'], elt['name pattern'])

