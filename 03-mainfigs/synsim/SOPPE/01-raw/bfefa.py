# testscript for brute-force EFA

import os
from matplotlib import pyplot as plt
import itertools as it
import numpy as np

# get BIOXTAS-RAW module
from sys import path as syspath
syspath.insert(0, "/Users/jwinnikoff/bin/raw/bioxtasraw/")
import RAWAPI as raw

# load .dat files in the neighborhood
files_dat = [fname for fname in os.listdir('.') if (".dat" in fname) and ("up" in fname)]

# ensure they are properly ordered by pressure
def press_get(fname):
    "Get pressure from SAXS shot fname"
    return int([seg for seg in fname.split('_') if "MPa" in seg][0][:-3])
    
#print(files_dat)
files_dat.sort(key = press_get)
#print(files_dat)

#Load the profiles (in order!)
profiles = raw.load_profiles(files_dat)

# explicit number of components
n_comp = 3

# build metaparameter matrix
# try everything
comp_ranges = [range for range in it.product(range(len(files_dat)), range(len(files_dat)))]
# filter for reasonable ranges
# beg before end! Min length!
comp_ranges = [range for range in comp_ranges if range[1] >= range[0] + 2]

# these are all possible ranges per above criteria, EFA-valid or not!
all_ranges = it.product(*[comp_ranges]*n_comp)
# all ranges must be unique
all_ranges = [ranges for ranges in all_ranges if len(set(ranges)) == n_comp]
print("Trying {n} range specs...".format(n=len(list(all_ranges))))

# run some EFAs
for these_ranges in all_ranges:
    print(np.array(these_ranges))
    efa_profiles, converged, conv_data, rotation_data = raw.efa(profiles, np.array(these_ranges), niter = 500)
    if converged: print(conv_data)

#plt.plot(profiles[0].getQ(), profiles[0].getI())
#plt.show()