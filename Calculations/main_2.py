import numpy as np
from mymain import my_main_2, preparation
from Methods.my_methods import write_detailed_txtfile, transform2latex_tab_2

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
	
	#######################
	# CdZnTe 2v2b^-
	#######################
	
	
	N_at116Cd_czt, N_at130Te_czt = my_main_2(eventfile_czt, scale_czt, x_range_czt, background[2], Q_val_ret=True)
	sum_czt_at116Cd = []
	sum_czt_at130Te = []
	iso = []
	for i in range(len(eventfile_czt)):
		sim = eventfile_czt[i].split('/')
		sim = sim[len(sim)-1]
		sim = sim.split('entries_')
		sim = sim[len(sim)-1]
		sim = sim.split('.')
		sim = sim[0]
		dir = './calc_solutions/czt/%s/' %(sim)
		transform2latex_tab_2(N_at116Cd_czt[i], dir, 'at116Cd_detsec_latextab.tex')
		transform2latex_tab_2(N_at130Te_czt[i], dir, 'at130Te_detsec_latextab.tex')
	
		descriptions = ['Contributions at Qvalues', 'N in 1/kg/yr ', 'row: detector (1 to 9); columns: sector (1 to 4)']
		write_detailed_txtfile(N_at116Cd_czt[i], dir, 'at116Cd_detsec.txt',  var_names=None, description=descriptions)
		write_detailed_txtfile(N_at130Te_czt[i], dir, 'at130Te_detsec.txt',  var_names=None, description=descriptions)
	
		sum_det_at116Cd = np.sum(N_at116Cd_czt[i], axis=1)		# creates sum of all 4 sectors
		sum_det_at130Te = np.sum(N_at130Te_czt[i], axis=1)		# creates sum of all 4 sectors
		new_data = [sum_det_at116Cd, sum_det_at130Te]
		transform2latex_tab_2(np.transpose(new_data), dir, 'contrib_det_latextab.tex')
		descriptions = ['Contributions at Qvalues', 'N in 1/kg/yr ', 'row: sum of all sectors in detector (1 to 9)']
		var_names = ['N_116Cd', 'N_130Te']
		write_detailed_txtfile(np.transpose(new_data), dir, 'contrib_det.txt',  var_names=var_names, description=descriptions)
	
		sum_czt_at116Cd.append(np.sum(sum_det_at116Cd))
		sum_czt_at130Te.append(np.sum(sum_det_at130Te))
		iso.append(sim)
	
	thissum_at116Cd = []
	thissum_at130Te = []
	
	thissum_at116Cd.append(np.sum(sum_czt_at116Cd))
	thissum_at130Te.append(np.sum(sum_czt_at130Te))
	sum_czt_at116Cd.append(thissum_at116Cd[0])
	sum_czt_at130Te.append(thissum_at130Te[0])
	
	iso.append('sum')
	dir = './calc_solutions/%s/czt/' %(savedir[i_dir])
	new_data = [iso, sum_czt_at116Cd, sum_czt_at130Te]
	transform2latex_tab_2(np.transpose(new_data), dir, 'contrib_czt_latextab.tex')
	descriptions = ['Contributions at Qvalues', 'N in 1/kg/yr ']
	var_names = ['Isotope', 'N_116Cd', 'N_130Te']
	write_detailed_txtfile(np.transpose(new_data), dir, 'contrib_czt.txt',  var_names=var_names, description=descriptions)
	
	#####################
	# detector coating
	#####################
	
	N_coating_at116Cd, N_coating_at130Te = my_main_2(eventfile_coating, scale_coating, x_range_coating, background[1], Q_val_ret=True)
	sum_coating_at116Cd = []
	sum_coating_at130Te = []
	iso = []
	for i in range(len(eventfile_coating)):
		sim = eventfile_coating[i].split('/')
		sim = sim[len(sim)-1]
		sim = sim.split('entries_')
		sim = sim[len(sim)-1]
		sim = sim.split('.')
		sim = sim[0]
		dir = './calc_solutions/%s/coating/%s/' %(savedir[i_dir], sim)
		transform2latex_tab_2(N_coating_at116Cd[i], dir, 'at116Cd_latextab.tex')
		transform2latex_tab_2(N_coating_at130Te[i], dir, 'at130Te_latextab.tex')
	
		descriptions = ['Contributions at Qvalues', 'N in 1/kg/yr ', 'row: detector (1 to 9); columns: sector (1 to 4)']
		write_detailed_txtfile(N_coating_at116Cd[i], dir, 'at116Cd_detsec.txt',  var_names=None, description=descriptions)
		write_detailed_txtfile(N_coating_at130Te[i], dir, 'at130Te_detsec.txt',  var_names=None, description=descriptions)
	
		sum_det_at116Cd = np.sum(N_coating_at116Cd[i], axis=1)		# creates sum of all 4 sectors
		sum_det_at130Te = np.sum(N_coating_at130Te[i], axis=1)		# creates sum of all 4 sectors
		new_data = [sum_det_at116Cd, sum_det_at130Te]
		transform2latex_tab_2(np.transpose(new_data), dir, 'contrib_det_latextab.tex')
		descriptions = ['Contributions at Qvalues', 'N in 1/kg/yr ', 'row: sum of all sectors in detector (1 to 9)']
		var_names = ['N_116Cd', 'N_130Te']
		write_detailed_txtfile(np.transpose(new_data), dir, 'contrib_det.txt',  var_names=var_names, description=descriptions)
	
		sum_coating_at116Cd.append(np.sum(sum_det_at116Cd))
		sum_coating_at130Te.append(np.sum(sum_det_at130Te))
		iso.append(sim)
	
	# differentation between glyptal and epoxy
	sum_glyptal_at116Cd = np.sum([sum_coating_at116Cd[0], sum_coating_at116Cd[2], sum_coating_at116Cd[4]])
	sum_glyptal_at130Te = np.sum([sum_coating_at130Te[0], sum_coating_at130Te[2], sum_coating_at130Te[4]])
	sum_epoxy_at116Cd = np.sum([sum_coating_at116Cd[1], sum_coating_at116Cd[3], sum_coating_at116Cd[5]])
	sum_epoxy_at130Te = np.sum([sum_coating_at130Te[1], sum_coating_at130Te[3], sum_coating_at130Te[5]])
	dir = './calc_solutions/%s/coating/' %(savedir[i_dir])
	new_data = [['glyptal', 'epoxy'], [sum_glyptal_at116Cd, sum_glyptal_at130Te], [sum_epoxy_at116Cd, sum_epoxy_at130Te]]
	transform2latex_tab_2(np.transpose(new_data), dir, 'contrib_coantingtype_latextab.tex')
	descriptions = ['Contributions at Qvalues', 'N in 1/kg/yr ']
	var_names = ['Coating', 'N_116Cd', 'N_130Te']
	write_detailed_txtfile(np.transpose(new_data), dir, 'contrib_coatingtype.txt',  var_names=var_names, description=descriptions)
	
	
	thissum_at116Cd.append(np.sum(sum_coating_at116Cd))
	thissum_at130Te.append(np.sum(sum_coating_at130Te))
	sum_coating_at116Cd.append(thissum_at116Cd[1])
	sum_coating_at130Te.append(thissum_at130Te[1])
	
	iso.append('sum')
	dir = './calc_solutions/%s/coating/' %(savedir[i_dir])
	new_data = [iso, sum_coating_at116Cd, sum_coating_at130Te]
	transform2latex_tab_2(np.transpose(new_data), dir, 'contrib_coating_latextab.tex')
	descriptions = ['Contributions at Qvalues', 'N in 1/kg/yr ']
	var_names = ['Isotope', 'N_116Cd', 'N_130Te']
	write_detailed_txtfile(np.transpose(new_data), dir, 'contrib_coating.txt',  var_names=var_names, description=descriptions)
	
	#################
	# plastic holder
	#################
	
	N_plastic_at116Cd, N_plastic_at130Te = my_main_2(eventfile_plastic, scale_plastic, x_range_plastic, background[0], Q_val_ret=True)
	sum_plastic_at116Cd = []
	sum_plastic_at130Te = []
	iso = []
	for i in range(len(eventfile_plastic)):
		sim = eventfile_plastic[i].split('/')
		sim = sim[len(sim)-1]
		sim = sim.split('entries_')
		sim = sim[len(sim)-1]
		sim = sim.split('.')
		sim = sim[0]
		dir = './calc_solutions/%s/plastic/%s/' %(savedir[i_dir], sim)
		transform2latex_tab_2(N_plastic_at116Cd[i], dir, 'at116Cd_latextab.tex')
		transform2latex_tab_2(N_plastic_at130Te[i], dir, 'at130Te_latextab.tex')
	
		descriptions = ['Contributions at Qvalues', 'N in 1/kg/yr ', 'row: detector (1 to 9); columns: sector (1 to 4)']
		write_detailed_txtfile(N_plastic_at116Cd[i], dir, 'at116Cd_detsec.txt',  var_names=None, description=descriptions)
		write_detailed_txtfile(N_plastic_at130Te[i], dir, 'at130Te_detsec.txt',  var_names=None, description=descriptions)
	
		sum_det_at116Cd = np.sum(N_plastic_at116Cd[i], axis=1)		# creates sum of all 4 sectors
		sum_det_at130Te = np.sum(N_plastic_at130Te[i], axis=1)		# creates sum of all 4 sectors
		new_data = [sum_det_at116Cd, sum_det_at130Te]
		transform2latex_tab_2(np.transpose(new_data), dir, 'contrib_det_latextab.tex')
		descriptions = ['Contributions at Qvalues', 'N in 1/kg/yr ', 'row: sum of all sectors in detector (1 to 9)']
		var_names = ['N_116Cd', 'N_130Te']
		write_detailed_txtfile(np.transpose(new_data), dir, 'contrib_det.txt',  var_names=var_names, description=descriptions)
	
		sum_plastic_at116Cd.append(np.sum(sum_det_at116Cd))
		sum_plastic_at130Te.append(np.sum(sum_det_at130Te))
		iso.append(sim)
	
	# differentation between plate and screw
	sum_plasticplate_at116Cd = np.sum([sum_coating_at116Cd[0], sum_coating_at116Cd[1]])
	sum_plasticplate_at130Te = np.sum([sum_coating_at130Te[0], sum_coating_at130Te[1]])
	sum_plasticscrew_at116Cd = np.sum([sum_coating_at116Cd[2], sum_coating_at116Cd[3]])
	sum_plasticscrew_at130Te = np.sum([sum_coating_at130Te[2], sum_coating_at130Te[3]])
	dir = './calc_solutions/%s/plastic/' %(savedir[i_dir]) 
	new_data = [['plate', 'screw'], [sum_glyptal_at116Cd, sum_glyptal_at130Te], [sum_epoxy_at116Cd, sum_epoxy_at130Te]]
	transform2latex_tab_2(np.transpose(new_data), dir, 'contrib_plastictype_latextab.tex')
	descriptions = ['Contributions at Qvalues', 'N in 1/kg/yr ']
	var_names = ['plastic_el', 'N_116Cd', 'N_130Te']
	write_detailed_txtfile(np.transpose(new_data), dir, 'contrib_plastictype.txt',  var_names=var_names, description=descriptions)
	
	
	thissum_at116Cd.append(np.sum(sum_plastic_at116Cd))
	thissum_at130Te.append(np.sum(sum_plastic_at130Te))
	sum_plastic_at116Cd.append(thissum_at116Cd[2])
	sum_plastic_at130Te.append(thissum_at130Te[2])
	
	iso.append('sum')
	dir = './calc_solutions/%s/plastic/' %(savedir[i_dir])
	new_data = [iso, sum_plastic_at116Cd, sum_plastic_at130Te]
	transform2latex_tab_2(np.transpose(new_data), dir, 'contrib_plastic_latextab.tex')
	descriptions = ['Contributions at Qvalues', 'N in 1/kg/yr ']
	var_names = ['Isotope', 'N_116Cd', 'N_130Te']
	write_detailed_txtfile(np.transpose(new_data), dir, 'contrib_plastic.txt',  var_names=var_names, description=descriptions)
	
	
	
	####################
	# all contributions
	####################
	dir = './calc_solutions/%s/' %(savedir[i_dir])
	thissum_at116Cd.append(np.sum(thissum_at116Cd))
	thissum_at130Te.append(np.sum(thissum_at130Te))
	new_data = [['CdZnTe', 'coating', 'plastic', 'sum'], thissum_at116Cd, thissum_at130Te]
	var_names = ['Detector_el', 'N_116Cd', 'N_130Te']
	transform2latex_tab_2(np.transpose(new_data), dir, 'contrib_all_latextab.tex')
	descriptions = ['Contributions at Qvalues', 'N in 1/kg/yr ']
	var_names = ['Isotope', 'N_116Cd', 'N_130Te']
	write_detailed_txtfile(np.transpose(new_data), dir, 'contrib_all.txt',  var_names=var_names, description=descriptions)
