import numpy as np
from mymain import my_main_3, preparation
from sim_events_3 import create_sumHist, create_allsumHist, save_sumHist

ev_dir  = ['rootfiles_noguard_06_09_2018', 
		'rootfiles_zcut_only_06_09_2018', 
		'rootfiles_guard_only_06_09_2018', 
		'rootfiles_guard+zcut_06_09_2018'
		]

savedir = ['noguard',
			'zcut_only',
			'guard_only',
			'guard+zcut'
			]

for i_dir in range(len(ev_dir)):
	eventfiles, background, scale, x_range = preparation(ev_dir[i_dir])
	eventfile_czt = eventfiles[2]
	scale_czt = scale[2]
	x_range_czt = x_range[2]
	
	eventfile_coating = eventfiles[1]
	scale_coating = scale[1]
	x_range_coating = x_range[1]
	
	eventfile_plastic = eventfiles[0]
	scale_plastic = scale[0]
	x_range_plastic = x_range[0]
	
	all_hists = []
	all_hists.append(my_main_3(eventfile_plastic, scale_plastic, x_range_coating, background[0], returns=True, savedir=savedir[i_dir]))
	all_hists.append(my_main_3(eventfile_coating, scale_coating, x_range_coating, background[1], returns=True, savedir=savedir[i_dir]))
	all_hists.append(my_main_3(eventfile_czt, scale_czt, x_range_czt, background[2], returns=True, savedir=savedir[i_dir]))
	
	hists = create_sumHist(all_hists, i_dir)
	all_sumhists = create_allsumHist(all_hists, i_dir, these_bins=250)
	save_sumHist(hists, all_sumhists, background, save_dir=savedir[i_dir])
