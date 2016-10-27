import debugger_skeleton
import filters as uf
import core_functions as cf
from matplotlib import pyplot as plt

# different debugger injection can be performed by doing a
# import module_of_interest
# module_of_interest.debugger = new_debugger

translator = {'w1488': 0,
              'w2561': 1
              }

source = uf.traverse_and_match("L:\\Users\\linghao\\Data for quantification\\Yeast\\NEW data for analysis",
                               matching_map=translator)

named_source = uf.name_channels(source, ['GFP', 'mCherry'])
GFP = cf.extract(named_source, 'GFP')
stabilized_GFP = cf.gamma_stabilize(GFP)
smoothed_GFP = cf.smooth(stabilized_GFP)
# should I even be considering wrapping the functions to make iterators out of them?
#  idea - since it's impossible to split the generators but we are interested only in stacks,
#  either build the pipeline within a round (filters always on all channels)
#  or iteratively modify the named_source


for payload in smoothed_GFP:
    plt.imshow(payload[6, :, :])
    plt.show()
