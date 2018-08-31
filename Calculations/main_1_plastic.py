import numpy as np
from mymain import my_main_1, preparation

eventfiles, background, scale, x_range = preparation()
eventfile_czt = eventfiles[0]
scale_czt = scale[0]
x_range_czt = x_range[0]

eventfile_coating = eventfiles[1]
scale_coating = scale[1]
x_range_coating = x_range[1]

eventfile_plastic = eventfiles[2]
scale_plastic = scale[2]
x_range_plastic = x_range[2]

my_main_1(eventfile_plastic, scale_plastic, x_range_plastic)
