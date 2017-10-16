from csv import writer as csv_writer

import examples.akshay_support
import imagepipe.traversals as uf
import imagepipe.wrapped_functions as wf
from imagepipe import core_functions as cf

# source = uf.Akshay_traverse("L:\\Common\\AKN\\IF test images nuc vs cyto\\tiff")
# source = uf.Akshay_traverse("L:\\Common\\AKN\\IF test images nuc vs cyto\\11-21-16 RPE IF untreat mps nut images")
# source_directory = '/run/user/1000/gvfs/smb-share:server=10.17.0.219,share=common/Common/AKN/IF test images nuc vs cyto/11-21-16 RPE IF untreat mps nut - Analysis run'  # because Linux
source_directory = "L:\\Common\\AKN\\IF test images nuc vs cyto\\11-21-16 RPE IF untreat mps nut images"  # because Windows
# source_directory = "C:\\Users\\Andrei\Desktop\\tmp\\11-21-16 RPE IF untreat mps nut images"  # because windows and smb are not collaborating today

source = examples.akshay_support.Akshay_traverse(source_directory)
named_source = uf.name_channels(source, ['DAPI', 'p53', 'p21'])

stablilized_1 = wf.gamma_stabilize(named_source, in_channel='DAPI', floor_method='min', alpha_clean=.5)
smoothed_1 = wf.smooth_2d(stablilized_1, in_channel='DAPI', smoothing_px=.5)  # corresponds to a pixel of 1500 nm

stablilized_2 = wf.gamma_stabilize(smoothed_1, in_channel='p53', floor_method='min', alpha_clean=.0)
smoothed_2 = wf.smooth_2d(stablilized_2, in_channel='p53', smoothing_px=.5)

stablilized_3 = wf.gamma_stabilize(smoothed_2, in_channel='p21', floor_method='5p', alpha_clean=.5)
smoothed_3 = wf.smooth_2d(stablilized_3, in_channel='p21', smoothing_px=.5)

binarized_nuclei = wf.robust_binarize(smoothed_3,
                                                           in_channel='DAPI',
                                                           out_channel=['nuclei'],
                                                           _dilation=0,
                                                           heterogeity_size=5, feature_size=50)

segmented_nuclei = wf.label_and_correct(binarized_nuclei,
                                                             in_channel=['nuclei', 'DAPI'],
                                                             out_channel='nuclei',
                                                             min_px_radius=15, min_intensity=20)

# p53 segmentation

p53_aq = wf.label_based_aq(segmented_nuclei,
                                                in_channel=['nuclei', 'p53'],
                                                out_channel=['av_p53', 'av_p53_pad'])

p53_o_n = wf.exclude_region(p53_aq,
                                                 in_channel=['nuclei', 'p53'],
                                                 out_channel='p53_o_n')

vor_seg = wf.voronoi_segment_labels(p53_o_n,
                                                         in_channel='nuclei',
                                                         out_channel='vor_segment')

p53_o_n_segmented = wf.robust_binarize(vor_seg,
                                                            in_channel='p53_o_n',
                                                            out_channel='p53_o_n_seg',
                                                            heterogeity_size=10, feature_size=100)

p53_seg_contacted = wf.in_contact(p53_o_n_segmented,
                                                       in_channel=['p53_o_n_seg', 'nuclei'],
                                                       out_channel=['p53_o_n_seg', '_'])

# TODO: test the elimination

p53_o_n_filtered = wf.filter_labels(p53_seg_contacted,
                                                         in_channel=['vor_segment', 'p53_o_n_seg'],
                                                         out_channel=['extra_nuclear_p53'],
                                                         min_feature_size=30)

p53_en_eq = wf.label_based_aq(p53_o_n_filtered,
                                                   in_channel=['extra_nuclear_p53', 'p53_o_n'],
                                                   out_channel=['av_en_p53', 'av_en_p53_pad'])

# p21 now

p21_aq = wf.label_based_aq(p53_en_eq,
                                                in_channel=['nuclei', 'p21'],
                                                out_channel=['nuc_p21', 'nuc_p21_pad'])

p21_o_n = wf.exclude_region(p21_aq,
                                                 in_channel=['nuclei', 'p21'],
                                                 out_channel='p21_o_n')

p21_o_n_segmented = wf.robust_binarize(p21_o_n,
                                                            in_channel='p21_o_n',
                                                            out_channel='p21_o_n_seg',
                                                            heterogeity_size=10, feature_size=100)

p21_seg_contacted = wf.in_contact(p21_o_n_segmented,
                                                       in_channel=['p21_o_n_seg', 'nuclei'],
                                                       out_channel=['p21_o_n_seg', '_'])

p21_o_n_filtered = wf.filter_labels(p21_seg_contacted,
                                                         in_channel=['vor_segment', 'p21_o_n_seg'],
                                                         out_channel=['extra_nuclear_p21'],
                                                         min_feature_size=30)

p21_en_eq = wf.label_based_aq(p21_o_n_filtered,
                                                   in_channel=['extra_nuclear_p21', 'p21_o_n'],
                                                   out_channel=['av_en_p21', 'av_en_p21_pad'])

running_render = examples.akshay_support.akshay_render(p21_en_eq,
                                                       in_channel=['name pattern', 'DAPI', 'p53', 'p21',
                                               'nuclei', 'vor_segment',
                                               'extra_nuclear_p53', 'av_p53_pad', 'av_en_p53_pad',
                                               'extra_nuclear_p21', 'nuc_p21_pad', 'av_en_p21_pad'],
                                                       out_channel='_',
                                                       save=False)

# after checking this pipeline, replace False with True

summary = examples.akshay_support.akshay_summarize(running_render,
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

