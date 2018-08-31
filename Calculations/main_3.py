import numpy as np
from mymain import my_main_3, preparation

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

all_hists = []
all_hists.append(my_main_3(eventfile_plastic, scale_plastic, x_range_coating, background[0], returns=True))
all_hists.append(my_main_3(eventfile_coating, scale_coating, x_range_coating, background[1], returns=True))
all_hists.append(my_main_3(eventfile_czt, scale_czt, x_range_czt, background[2], returns=True))

hists = create_sumHist(all_hists, background)
all_sumhists = create_allsumHist(all_hists, background)
save_sumHist(hists, all_sumhists, background)
