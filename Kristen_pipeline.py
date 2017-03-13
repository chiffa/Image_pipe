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




