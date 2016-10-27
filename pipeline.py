import debugger_skeleton
import filters as uf

# different debugger injection can be performed by doing a
# import module_of_interest
# module_of_interest.debugger = new_debugger

translator = {'w1488': 0,
              'w2561': 1
              }

source = uf.traverse_and_match("L:\\Users\\linghao\\Data for quantification\\Yeast\\NEW data for analysis",
                   matching_map=translator)

named_source = uf.name_channels(source, ['GFP', 'mCherry'])

for payload in named_source:
    print payload['name pattern']
    print payload['GFP'].shape