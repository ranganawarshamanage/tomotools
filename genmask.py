# making mask for G3 ribosome

import emda2.emda_methods2 as em
from emda2.core import iotools
import argparse
import numpy as np

parser = argparse.ArgumentParser(description='')
parser.add_argument('--firstmap', type=str, required=True, help='map1 file mrc/map')
parser.add_argument('--secondmap', type=str, required=True, help='map2 file mrc/map')
#parser.add_argument('--lowpass_resolution', type=float, required=True, help='Specify lowpass resolution')
#parser.add_argument('--outputmapname', type=str, default=None, required=False, help='outputmap name mrc/map')
#parser.add_argument('--binfactor', type=float, required=True, help='Specify the bin factor')
args = parser.parse_args()


m1 = iotools.Map(name=args.firstmap)
m1.read()

m2 = iotools.Map(name=args.secondmap)
m2.read()


# assuming that mask1 is the bigger mask

outputmapname = 'emda_genmask.mrc'
m3 = iotools.Map(name=outputmapname)
m3.cell = m1.cell
m3.arr = (m1.workarr - m2.workarr) > 0.9
m3.axorder = m1.axorder
m3.write()