from csv import writer as csv_writer

import examples.akshay_support
import imagepipe.traversals as uf
import imagepipe.wrapped_functions as wf
from imagepipe import core_functions as cf

# source = uf.Akshay_traverse("L:\\Common\\AKN\\IF test images nuc vs cyto\\tiff")
# source = uf.Akshay_traverse("L:\\Common\\AKN\\IF test images nuc vs cyto\\11-21-16 RPE IF untreat mps nut images")
# source_directory = "L:\\Common\\AKN\\IF test images nuc vs cyto\\tiff"
# source_directory = "L:\\Users\\akshay\\Raw files HCT116p53YFP for Andrei quantification\\TIFF"
source_directory = "C:\\Users\\Andrei\\Desktop\\tmp\\TIFF"
# source_directory = '/run/user/1000/gvfs/smb-share:server=10.17.0.219,share=common/Common/AKN/IF test images nuc vs cyto/11-21-16 RPE IF untreat mps nut - Analysis run'  # because Linux
# source_directory = "L:\\Common\\AKN\\IF test images nuc vs cyto\\11-21-16 RPE IF untreat mps nut images"  # because Windows
# source_directory = "C:\\Users\\Andrei\Desktop\\tmp\\11-21-16 RPE IF untreat mps nut images"  # because windows and smb are not collaborating today


pipeline = uf.color_based_traversal(source_directory, coding_separator=' ', group_anchor=1)
pipeline = uf.name_channels(pipeline, ['transmission', 'DAPI', 'p53'])

pipeline = wf.gamma_stabilize(pipeline, in_channel='DAPI', floor_method='min', alpha_clean=.5)
pipeline = wf.smooth_2d(pipeline, in_channel='DAPI', smoothing_px=.5)  # corresponds to a pixel of 1500 nm

pipeline = wf.gamma_stabilize(pipeline, in_channel='p53', floor_method='min', alpha_clean=.0)
pipeline = wf.smooth_2d(pipeline, in_channel='p53', smoothing_px=.5)


pipeline = wf.robust_binarize(pipeline,
                                                           in_channel='DAPI',
                                                           out_channel=['nuclei'],
                                                           _dilation=0,
                                                           heterogeity_size=5, feature_size=50)

pipeline = wf.label_and_correct(pipeline,
                                                             in_channel=['nuclei', 'DAPI'],
                                                             out_channel='nuclei',
                                                             min_px_radius=15, min_intensity=20)

# p53 segmentation

pipeline = wf.label_based_aq(pipeline,
                                                in_channel=['nuclei', 'p53'],
                                                out_channel=['av_p53', 'av_p53_pad'])

pipeline = wf.exclude_region(pipeline,
                                                 in_channel=['nuclei', 'p53'],
                                                 out_channel='p53_o_n')

pipeline = wf.voronoi_segment_labels(pipeline,
                                                         in_channel='nuclei',
                                                         out_channel='vor_segment')

pipeline = wf.robust_binarize(pipeline,
                                                            in_channel='p53_o_n',
                                                            out_channel='p53_o_n_seg',
                                                            heterogeity_size=10, feature_size=100)

pipeline = wf.in_contact(pipeline,
                                                       in_channel=['p53_o_n_seg', 'nuclei'],
                                                       out_channel=['p53_o_n_seg', '_'])

# TODO: test the elimination

pipeline = wf.filter_labels(pipeline,
                                                         in_channel=['vor_segment', 'p53_o_n_seg'],
                                                         out_channel=['extra_nuclear_p53'],
                                                         min_feature_size=30)

pipeline = wf.label_based_aq(pipeline,
                                                   in_channel=['extra_nuclear_p53', 'p53_o_n'],
                                                   out_channel=['av_en_p53', 'av_en_p53_pad'])


pipeline = examples.akshay_support.akshay_partial_render(pipeline,
                                                in_channel=['name pattern', 'DAPI', 'p53',
                                                'nuclei', 'vor_segment',
                                                'extra_nuclear_p53', 'av_p53_pad', 'av_en_p53_pad'],
                                                out_channel='_',
                                                save=False)

# after checking this pipeline, replace False with True

pipeline = examples.akshay_support.akshay_partial_summarize(pipeline,
                                            in_channel=['name pattern', 'group id', 'av_p53', 'av_en_p53'],
                                            out_channel='_',
                                            output='akshay_analysis_results.csv')

with open('akshay_analysis_results.csv', 'wb') as output_file:
        writer = csv_writer(output_file)
        writer.writerow(['file', 'group id', 'cell no', 'nuclear p53',
                         'cellular p53', 'nuclear p21', 'cellular p21'])

for i, elt in enumerate(pipeline):
    print 'operation %s analyzed group %s - image %s' % (i, elt['group id'], elt['name pattern'])

