from csv import writer as csv_writer
from time import time, strftime

import examples.xi_support
import imagepipe.traversals as uf
import imagepipe.wrapped_functions as wf
from imagepipe import core_functions as cf

translator = {'C2': 0,
              'C1': 1
              }

source = examples.xi_support.xi_traverse('L:\\Users\\xi\\Quantification of colocalization',
                                         matching_map=translator)

# source = uf.xi_traverse('L:\\Users\\andrei\\2017\\Xi Data',
#                         matching_map=translator)

named_source = uf.name_channels(source, ['GFP', 'mCherry'])

stabilized_GFP = wf.gamma_stabilize(named_source, in_channel='GFP', floor_method='min', alpha_clean=.5)

smoothed_GFP = wf.smooth(stabilized_GFP, in_channel='GFP', smoothing_px=.5)

stabilized_mCh = wf.gamma_stabilize(smoothed_GFP, in_channel='mCherry', floor_method='min', alpha_clean=.5)

smoothed_mCh = wf.smooth(stabilized_mCh, in_channel='mCherry', smoothing_px=.5)

med_norm_GFP = wf.locally_normalize(smoothed_mCh,
                                                         in_channel='GFP',
                                                         )

med_norm_mCh = wf.locally_normalize(med_norm_GFP,
                                                         in_channel='mCherry',
                                                         )

projected_GFP = wf.max_projection(med_norm_mCh,
                                                       in_channel='GFP',
                                                       out_channel='projected_GFP')

projected_mCh = wf.max_projection(projected_GFP,
                                                       in_channel='mCherry',
                                                       out_channel='projected_mCh')


binarized_GFP = wf.robust_binarize(projected_mCh,
                                                        in_channel='projected_GFP',
                                                        out_channel='cell_tags',
                                                        heterogeity_size=5,
                                                        feature_size=70,
                                                        )

segmented_GFP = wf.improved_watershed(binarized_GFP,
                                                           in_channel=['cell_tags', 'projected_mCh'],
                                                           out_channel='pre_cell_labels',
                                                           expected_separation=100)

qualifying_GFP = wf.qualifying_gfp(segmented_GFP,
                                                        in_channel='projected_GFP',
                                                        out_channel='qualifying_GFP')

average_GFP = wf.average_qualifying_value_per_region(qualifying_GFP,
                                                                          in_channel=['pre_cell_labels', 'projected_GFP', 'qualifying_GFP'],
                                                                          out_channel=['average_GFP', 'average_GFP_pad'])

pre_render = examples.xi_support.xi_pre_render(average_GFP,
                                               in_channel=['name pattern',
                                         'projected_GFP', 'qualifying_GFP',
                                         'pre_cell_labels',
                                         'average_GFP_pad', 'projected_mCh',
                                         'mCherry', 'GFP', 'group id'],
                                               out_channel='_',
                                               save=True)

# since the resolution for the mitochondria is so much lower when compared to yeast cells, we will
# have to perform a bit more straightforward of cutting and

with open('xi_analys_results.csv', 'wb') as output_file:
        writer = csv_writer(output_file)
        writer.writerow(['file', 'time', 'correlation coeff',
                         'median GFP', 'average GFP', 'linreg slope', 'linreg rvalue',
                         'linreg pvalue'])

prev_time = time()

for primary_namespace in pre_render:
    print '%s - analyzed %s - %s in %s' % (strftime('%X %x'),
                                           primary_namespace['name pattern'],
                                           primary_namespace['group id'],
                                           time() - prev_time)
    prev_time = time()
