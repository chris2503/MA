import numpy as np
from mymain import my_main_1, preparation

ev_dir  = ['rootfiles_noguard_06_09_2018', 
		'rootfiles_zcut_only_06_09_2018', 
		'rootfiles_guard_only_06_09_2018', 
		'rootfiles_guard+zcut_06_09_2018'
		]

# hier in Bash Ersetzungen machen!!!
savedir = ['noguard',
			'zcut_only',
			'guard_only',
			'guard+zcut'
			]

eventfiles, background, scale, x_range = preparation(ev_dir[3])
background_czt = background[2]
eventfile_czt = eventfiles[2]
scale_czt = scale[2]
x_range_czt = x_range[2]

background_coating = background[1]
eventfile_coating = eventfiles[1]
scale_coating = scale[1]
x_range_coating = x_range[1]

background_plastic = background[0]
eventfile_plastic = eventfiles[0]
scale_plastic = scale[0]
x_range_plastic = x_range[0]

my_main_1(eventfile_template, scale_template, x_range_template, savedir=savedir[3])
