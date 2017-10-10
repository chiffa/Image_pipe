"""
This module essentially wraps the functions to be compatible with the usage inside the pipelines
"""

from imagepipe.core_functions import generator_wrapper
import imagepipe.raw_functions as rf


gamma_stabilize = generator_wrapper(in_dims=(None,))(rf.gamma_stabilize)
smooth = generator_wrapper(in_dims=(3,))(rf.smooth)
smooth_2d = generator_wrapper(in_dims=(2,))(rf.smooth_2d)
sum_projection = generator_wrapper(in_dims=(3,), out_dims=(2,))(rf.sum_projection)
max_projection = generator_wrapper(in_dims=(3,), out_dims=(2,))(rf.max_projection)
random_walker_binarize = generator_wrapper(in_dims=(2,))(rf.random_walker_binarize)
robust_binarize = generator_wrapper(in_dims=(2,), out_dims=(2,))(rf.robust_binarize)
filter_labels = generator_wrapper(in_dims=(2, 2), out_dims=(2,))(rf.filter_labels)
voronoi_segment_labels = generator_wrapper(in_dims=(2,))(rf.voronoi_segment_labels)
exclude_region = generator_wrapper(in_dims=(2, 2), out_dims=(2,))(rf.exclude_region)
in_contact = generator_wrapper(in_dims=(2, 2))(rf.in_contact)
improved_watershed = generator_wrapper(in_dims=(2, 2), out_dims=(2,))(rf.improved_watershed)
label_and_correct = generator_wrapper(in_dims=(2, 2,), out_dims=(2,))(rf.label_and_correct)
qualifying_gfp = generator_wrapper(in_dims=(2,))(rf.qualifying_gfp)
label_based_aq = generator_wrapper(in_dims=(2, 2), out_dims=(1, 2))(rf.label_based_aq)
average_qualifying_value_per_region = generator_wrapper(in_dims=(2, 2, 2), out_dims=(1, 2))(rf.average_qualifying_value_per_region)
detect_upper_outliers = generator_wrapper(in_dims=(1,), out_dims=(1, 1, None))(
    rf.detect_upper_outliers)
paint_mask = generator_wrapper(in_dims=(2, 1), out_dims=(2,))(rf.paint_mask)
mask_filter_2d = generator_wrapper(in_dims=(2, 2), out_dims=(2,))(rf.mask_filter_2d)
clear_based_on_2d_mask = generator_wrapper(in_dims=(3, 2), out_dims=(3,))(rf.clear_based_on_2d_mask)
binarize_3d = generator_wrapper()(rf.binarize_3d)
volume_mqvi = generator_wrapper(in_dims=(3, 3), out_dims=(None,))(rf.volume_mqvi)
volume_aqvi = generator_wrapper(in_dims=(3, 3), out_dims=(None,))(rf.volume_aqvi)
_3d_mask_from_2d_mask = generator_wrapper(in_dims=(3, 2), out_dims=(3,))(rf.otsu_tresholding)
binarize_2d = generator_wrapper(in_dims=(2,), out_dims=(2,))(rf.binarize_2d)
agreeing_skeletons = generator_wrapper(in_dims=(2, 2), out_dims=(2,))(rf.agreeing_skeletons)
classify_fragmentation_for_mitochondria = generator_wrapper(in_dims=(2, 2), out_dims=(None, 2, 2, 2))(
    rf.classify_fragmentation_for_mitochondria)
locally_normalize = generator_wrapper(in_dims=(3,), out_dims=(3,))(rf.locally_normalize)
