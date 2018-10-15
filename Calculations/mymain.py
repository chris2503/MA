import numpy as np
from sim_events_3 import *
from Methods.my_methods import read_File, convert_str2num, write_txtfile, write_detailed_txtfile
import uncertainties
import uncertainties.unumpy as unp
from uncertainties.unumpy import(nominal_values as noms, std_devs as stds)

##############################################
# main
##############################################
def preparation(thisdir=None):
	kprint=None
	if thisdir:
		eventfile_coating = ['../%s/edep_entries/entries_40K_glyptal.txt' %(thisdir),
							 '../%s/edep_entries/entries_40K_epoxy.txt' %(thisdir),
							 '../%s/edep_entries/entries_232Th_glyptal.txt' %(thisdir),
							 '../%s/edep_entries/entries_232Th_epoxy.txt' %(thisdir),
							 '../%s/edep_entries/entries_238U_glyptal.txt' %(thisdir),
							 '../%s/edep_entries/entries_238U_epoxy.txt' %(thisdir)]
		
		eventfile_czt = ['../%s/edep_entries/entries_114Cd.txt' %(thisdir),
						 '../%s/edep_entries/entries_116Cd.txt' %(thisdir),
						 '../%s/edep_entries/entries_70Zn.txt' %(thisdir),
						 '../%s/edep_entries/entries_128Te.txt' %(thisdir),
						 '../%s/edep_entries/entries_130Te.txt' %(thisdir)]
		
		eventfile_plastic = ['../%s/edep_entries/entries_232Th_dplate.txt' %(thisdir),
							 '../%s/edep_entries/entries_238U_dplate.txt' %(thisdir),
							 '../%s/edep_entries/entries_232Th_dscrew.txt' %(thisdir),
							 '../%s/edep_entries/entries_238U_dscrew.txt' %(thisdir)]

	else:
		eventfile_coating = ['../rootfiles/edep_entries/entries_40K_glyptal.txt',
							 '../rootfiles/edep_entries/entries_40K_epoxy.txt',
							 '../rootfiles/edep_entries/entries_232Th_glyptal.txt',
							 '../rootfiles/edep_entries/entries_232Th_epoxy.txt',
							 '../rootfiles/edep_entries/entries_238U_glyptal.txt',
							 '../rootfiles/edep_entries/entries_238U_epoxy.txt']
		
		eventfile_czt = ['../rootfiles/edep_entries/entries_114Cd.txt',
						 '../rootfiles/edep_entries/entries_116Cd.txt',
						 '../rootfiles/edep_entries/entries_70Zn.txt',
						 '../rootfiles/edep_entries/entries_128Te.txt',
						 '../rootfiles/edep_entries/entries_130Te.txt']
		
		eventfile_plastic = ['../rootfiles/edep_entries/entries_232Th_dplate.txt',
							 '../rootfiles/edep_entries/entries_238U_dplate.txt',
							 '../rootfiles/edep_entries/entries_232Th_dscrew.txt',
							 '../rootfiles/edep_entries/entries_238U_dscrew.txt']
	
	
	x_range_coating = 1e4
	x_range_czt = 3e3
	
	years = 60*60*24*365 							# 1 year
	
	m_det_or = np.array([34.21, 35.15, 34.63, 33.91, 35.50, 35.50, 35.50, 35.50, 35.50]) * 1e-3			#kg
	m_det = np.mean(m_det_or)
	
	
	######################
	# Delrin layer and screws
	######################
	# 1) Layer
	V_T = 2*(47.8-2)* 2*46 * 1 - 9* (5+5.2)**2 * 1 	# mm^3		# top plate, upper bottom plate
	V_B = 2*47.8 * 2*46 * 2 - 9*5**2				# mm^3		# bottom plate, implemented as one plate instead of two single plates
	V_Layer = (2*V_T + V_B) *1e-9 # m^3
	
	rho_Layer = 1.43 * 1e3 # kg/m^3, Polyoxymethylene
	m_layer = V_Layer * rho_Layer
	
	# 2) screws
	V_screws = 4*(1/2)**2*np.pi*17 *1e-9 # m^3			# plastic holder screws
	
	rho_screws = 1.14 * 1e3 # kg/m^3
	m_screws = V_screws * rho_screws
	
	
	datafile = './Lists/activities_plastic.txt'
	this_data = read_File(datafile)
	iso_list = this_data[:,0]									# get isotope
	sp_ac_dplate = convert_str2num(this_data[:,1])				# get specific activity
	sp_ac_dscrews = convert_str2num(this_data[:,2])				# get specific activity
	
	N_norm_plastic = []
	scale_plastic = []
	n_chain = [10, 14]										# 232Th-, 238U-chain
	for i in range(len(iso_list)):
		N_norm_plastic.append(sp_ac_dplate[i] * m_layer * years * n_chain[i] * 1/(9*m_det))			# m_det = [kg]
	for i in range(len(iso_list)):
		N_norm_plastic.append(sp_ac_dscrews[i] * m_screws * years * n_chain[i] * 1/(9*m_det))
	
	N_simEv = 1e6										# 1 Mio simulated Events
	for i in range(len(N_norm_plastic)):
		scale_plastic.append(1/4 * N_norm_plastic[i]/N_simEv)
	
	################
	#### Coating
	################
	
	m_detpaint_glyptal_or = np.array([0.06, 0.07, 0.06, 0.07]) * 1e-3									#kg
	m_detpaint_epoxy_or = np.array([0.14, 0.13, 0.14, 0.12, 0.16]) * 1e-3								#kg
	
	
	m_detpaint_epoxy = np.mean(m_detpaint_epoxy_or)
	m_detpaint_glyptal = np.mean(m_detpaint_glyptal_or)
	if kprint:
		print("\n"+10*"#")
		print("\nDetector mean mass: %f" %(m_det))
		print("\nepoxy mean mass: %f" %(m_detpaint_epoxy))
		print("glyptal mean mass: %f" %(m_detpaint_glyptal))
	
	rho_glyptal= 1.441*1e3		# kg/m^3
	# source: http://www.physics.purdue.edu/primelab/safety/MSDS/SDS/epoxy%20Red%20Enamel%20%20%201201B.pdf
	
	rho_epoxy = 1.2*1e3 		# kg/m^3
	# source: https://www.google.de/url?sa=t&rct=j&q=&esrc=s&source=web&cd=4&ved=0ahUKEwixx9XCpKrZAhXL2KQKHUu_COYQFghCMAM&url=https%3A%2F%2Fwww.epoxies.com%2F_resources%2Fcommon%2Fuserfiles%2Ffile%2F20-3001NC.pdf&usg=AOvVaw2UKPhVUDbO2-M9trv9BrR-
	args_0 = (m_detpaint_epoxy, rho_epoxy)
	
	x_D = 10.2*1e-3				# m
	z_D = 16.0*1e-3				# m
	
	
	O = x_D**2 + z_D*x_D*4
	V_epoxy = m_detpaint_epoxy/rho_epoxy
	d_epoxy = V_epoxy/O
	
	V_glyptal = m_detpaint_glyptal/rho_glyptal
	d_glyptal = V_glyptal/O
	
	n_det_glyptal = 4										# number of with glyptal coated detectors
	n_det_epoxy = 5											# number of with epoxy coated detectors
	n_chain = [1, 10, 14]									# 40K-, 232Th-, 238U-chain
	
	datafile = './Lists/activities_detpaint.txt'
	this_data = read_File(datafile)
	iso_list = this_data[:,0]									# get isotope
	sp_ac_epoxy = convert_str2num(this_data[:,1])				# get specific activity
	sp_ac_glyptal = convert_str2num(this_data[:,2])				# get specific activity
	
	N_norm_coating = []
	scale_coating = []
	for i in range(len(iso_list)):
		N_norm_coating.append(sp_ac_glyptal[i] * m_detpaint_glyptal * years * n_chain[i] * 1/(n_det_glyptal* m_det))
		N_norm_coating.append(sp_ac_epoxy[i] * m_detpaint_epoxy * years * n_chain[i] * 1/(n_det_epoxy* m_det))
	
	N_simEv = 1e6											# 1 Mio simulated Events
	for i in range(len(N_norm_coating)):
		scale_coating.append(1/4 * N_norm_coating[i]/N_simEv)
	
	###########################
	#### CZT
	###########################
	datafile = './calc_solutions/calculated_events_2.txt'
	this_data = read_File(datafile)
	
	iso_list = this_data[:,0]								# get isotope
	N_norm_czt = convert_str2num(this_data[:,1])			# get norming factors, number convertion necessary
	N_simEv = 1e6											# 1 Mio simulated Events
	scale_czt = []								
	
	for i in range(len(N_norm_czt)):
		scale_czt.append(1/4 * N_norm_czt[i] / N_simEv) # It was calculated for 9 detectors

	background = ['plastic', 'coating', 'czt']
	eventfile = [eventfile_plastic, eventfile_coating, eventfile_czt]
	scale = [scale_plastic, scale_coating, scale_czt]
	x_range = [x_range_coating, x_range_coating, x_range_czt]

	return eventfile, background, scale, x_range 

#########################################
#########################################
def my_main_1(eventfile, scale, x_range, background=None, savedir=None):
	# creating plots
	# 1) single plots
	bins = 250
	if background:
		for i_file in range(len(eventfile)):
			data = read_evData(eventfile[i_file])
			my_secHists = create_sector_hists(data, scale[i_file], k=background)

			my_detHists = create_det_hists(my_secHists, these_bins=bins, k=background)
			my_sumhists = create_sumsecHist(my_secHists, these_bins=bins, k=background)
			
			iso_hist = create_iso_hist(my_detHists, eventfile[i_file], these_bins=bins, k=background)
			iso_sumhist = create_sumdetHist(my_detHists, these_bins=bins, k=background)

			if savedir:
				save_single_sector_hists(my_secHists, eventfile[i_file], x_range, k=background, save_dir=savedir)
				save_single_det_hists(my_detHists, my_sumhists, eventfile[i_file], x_range, k=background, save_dir=savedir)
				save_single_iso_hists(iso_hist, iso_sumhist, eventfile[i_file], x_range, k=background, save_dir=savedir)

			else:
				save_single_sector_hists(my_secHists, eventfile[i_file], x_range, k=background)
				save_single_det_hists(my_detHists, my_sumhists, eventfile[i_file], x_range, k=background)
				save_single_iso_hists(iso_hist, iso_sumhist, eventfile[i_file], x_range, k=background)

			delete_iso_sumHist(iso_sumhist, k=background)
			delete_all_sumdetHists(my_sumhists, k=background)
			delete_iso_Hist(iso_hist, eventfile[i_file], k=background)
			delete_all_detHists(my_detHists, k=background)
			delete_all_sectorHists(my_secHists, k=background)
		 	
		for i_file in range(len(eventfile)):
			data = read_evData(eventfile[i_file])
			my_depHists = create_dep_secHist(data, scale[i_file], these_bins=bins, k=background)
			if savedir:
				save_dep_sec_hists(my_depHists, eventfile[i_file], k=background, save_dir=savedir)
			else:
				save_dep_sec_hists(my_depHists, eventfile[i_file], k=background)
			delete_all_detHists(my_depHists, k=background)

			my_depHist = create_dep_detHist(data, eventfile[i_file], scale[i_file], k=background)
			if savedir:
				save_dep_det_hists(my_depHist, eventfile[i_file], k=background, save_dir=savedir)
			else:
				save_dep_det_hists(my_depHist, eventfile[i_file], k=background)
			delete_iso_Hist(my_depHist, eventfile[i_file], k=background)

			my_depHists_2 = create_dep_secHist(data, scale[i_file], k=background)
			if savedir:
				save_dep_sec_hists_2(my_depHists_2, eventfile[i_file], k=background, save_dir=savedir)
			else:
				save_dep_sec_hists_2(my_depHists_2, eventfile[i_file], k=background)
			delete_all_detHists(my_depHists_2, k=background)

	else:
		for i_file in range(len(eventfile)):
			data = read_evData(eventfile[i_file])
			my_secHists = create_sector_hists(data, scale[i_file])
			if savedir:
				save_single_sector_hists(my_secHists, eventfile[i_file], x_range, save_dir=savedir)
			else:
				save_single_sector_hists(my_secHists, eventfile[i_file], x_range)
			delete_all_sectorHists(my_secHists)

			my_secHists = create_sector_hists(data, scale[i_file])
			my_detHists = create_det_hists(my_secHists)
			my_sumhists = create_sumsecHist(my_secHists)
			if savedir:
				save_single_det_hists(my_detHists, my_sumhists, eventfile[i_file], x_range, save_dir=savedir)
			else:
				save_single_det_hists(my_detHists, my_sumhists, eventfile[i_file], x_range)
			delete_all_sumdetHists(my_sumhists)
			delete_all_detHists(my_detHists)
			delete_all_sectorHists(my_secHists)

			my_secHists = create_sector_hists(data, scale[i_file])
			my_detHists = create_sumsecHist(my_secHists, hcolor=True)
			iso_hist = create_iso_hist(my_detHists, eventfile[i_file])
			iso_sumhist = create_sumdetHist(my_detHists)
			if savedir:
				save_single_iso_hists(iso_hist, iso_sumhist, eventfile[i_file], x_range, save_dir=savedir)
			else:
				save_single_iso_hists(iso_hist, iso_sumhist, eventfile[i_file], x_range, save_dir=savedir)
			delete_iso_Hist(iso_hist, eventfile[i_file])
			delete_all_detHists(my_detHists)
			delete_all_sectorHists(my_secHists)
			delete_iso_sumHist(iso_sumhist)
		 	
		for i_file in range(len(eventfile)):
			data = read_evData(eventfile[i_file])
			my_depHists = create_dep_secHist(data, scale[i_file], k=background)
			if savedir:
				save_dep_sec_hists(my_depHists, eventfile[i_file], k=background, save_dir=savedir)
			else:
				save_dep_sec_hists(my_depHists, eventfile[i_file], k=background)
			delete_all_detHists(my_depHists)

			my_depHist = create_dep_detHist(data, eventfile[i_file], scale[i_file], k=background)
			if savedir:
				save_dep_det_hists(my_depHist, eventfile[i_file], k=background, save_dir=savedir)
			else:
				save_dep_det_hists(my_depHist, eventfile[i_file], k=background, save_dir=savedir)
			delete_iso_Hist(my_depHist, eventfile[i_file])

			my_depHists_2 = create_dep_secHist(data, scale[i_file], k=background)
			if savedir:
				save_dep_sec_hists_2(my_depHists_2, eventfile=eventfile[i_file], k=background, save_dir=savedir)
			else:
				save_dep_sec_hists_2(my_depHists_2, eventfile=eventfile[i_file], k=background)
			delete_all_detHists(my_depHists_2)

	# Heatmaps
	for i_file in range(len(eventfile)):
		data = read_evData(eventfile[i_file])
		my_depHists_heat = create_dep_secHist(data, scale[i_file])
		if savedir:
			save_dep_heatmap(my_depHists_heat, eventfile[i_file], save_dir=savedir)
		else:
			save_dep_heatmap(my_depHists_heat, eventfile[i_file])
		delete_all_detHists(my_depHists_heat)


def my_main_2(eventfile, scale, x_range, background, Q_val_ret=None, Q_val_ret_2=None, savedir=None):
	# 2) Values
	all_contrib_at116Cd = []
	all_contrib_at130Te = []

	if not Q_val_ret_2 and not Q_val_ret:
		for i_file in range(len(eventfile)):
			thiscase = background+str(i_file)
			data = read_evData(eventfile[i_file])

			my_secHists , contrib_at116Cd, contrib_at116Cd_err, contrib_at130Te, contrib_at130Te_err = create_sector_hists(data, scale[i_file], k=thiscase, Q_val_returns=True)
			sum = 0
			contrib_at116Cd = unp.uarray(contrib_at116Cd, contrib_at116Cd_err)
			#print(type(contrib_at116Cd))
			#print('116Cd Beiträge \n', contrib_at116Cd)
			for i in range(len(contrib_at116Cd)):
				for j in range(len(contrib_at116Cd[i])):
					sum = sum + contrib_at116Cd[i][j]
			all_contrib_at116Cd.append(sum)
			sum = 0
			contrib_at130Te = unp.uarray(contrib_at130Te, contrib_at130Te_err)
			#print('130Te Beiträge \n', contrib_at130Te)
			for i in range(len(contrib_at130Te)):
				for j in range(len(contrib_at130Te[i])):
					sum = sum + contrib_at130Te[i][j]
			all_contrib_at130Te.append(sum)
			delete_all_sectorHists(my_secHists)

		new_data = np.array([eventfile, all_contrib_at116Cd, all_contrib_at130Te])
		descriptions = ['Contributions at Qvalues', 'N in 1/kg/keV/yr ']
		var_names = ['File', 'N_at116Cd', 'N_at130Te']
		if savedir:
			write_detailed_txtfile(np.transpose(new_data), var_names, descriptions, './calc_solutions/%s/' %(savedir), 'events_at_Qvalues_%s.txt' %(background))
		else:
			write_detailed_txtfile(np.transpose(new_data), var_names, descriptions, './calc_solutions/', 'events_at_Qvalues_%s.txt' %(background))

	elif Q_val_ret:
		for i_file in range(len(eventfile)):
			thiscase = background+str(i_file)
			data = read_evData(eventfile[i_file])

			my_secHists , contrib_at116Cd, contrib_at116Cd_err, contrib_at130Te, contrib_at130Te_err = create_sector_hists(data, scale[i_file], k=thiscase, Q_val_returns=True)
			sum = 0
			contrib_at116Cd = unp.uarray(contrib_at116Cd, contrib_at116Cd_err)
			all_contrib_at116Cd.append(contrib_at116Cd)

			contrib_at130Te = unp.uarray(contrib_at130Te, contrib_at130Te_err)
			all_contrib_at130Te.append(contrib_at130Te)
			delete_all_sectorHists(my_secHists)
		return all_contrib_at116Cd, all_contrib_at130Te

	elif Q_val_ret_2:
		for i_file in range(len(eventfile)):
			thiscase = background+str(i_file)
			data = read_evData(eventfile[i_file])

			my_secHists , contrib_at116Cd, contrib_at116Cd_err, contrib_at130Te, contrib_at130Te_err = create_sector_hists(data, scale[i_file], k=thiscase, Q_val_returns=True)
			sum = 0
			contrib_at116Cd = unp.uarray(contrib_at116Cd, contrib_at116Cd_err)
			#print(type(contrib_at116Cd))
			#print('116Cd Beiträge \n', contrib_at116Cd)
			for i in range(len(contrib_at116Cd)):
				for j in range(len(contrib_at116Cd[i])):
					sum = sum + contrib_at116Cd[i][j]
			all_contrib_at116Cd.append(sum)
			sum = 0
			contrib_at130Te = unp.uarray(contrib_at130Te, contrib_at130Te_err)
			#print('130Te Beiträge \n', contrib_at130Te)
			for i in range(len(contrib_at130Te)):
				for j in range(len(contrib_at130Te[i])):
					sum = sum + contrib_at130Te[i][j]
			all_contrib_at130Te.append(sum)
			delete_all_sectorHists(my_secHists)
		return all_contrib_at116Cd, all_contrib_at130Te
		


def my_main_3(eventfile, scale, x_range, background, returns=None, savedir=None):
	# 3) combined plots
	bins = 250
	all_isohist =[]
	for i_file in range(len(eventfile)):
		thiscase = background+str(i_file)
		data = read_evData(eventfile[i_file])

		my_secHists = create_sector_hists(data, scale[i_file], these_bins=bins, k=thiscase)
		my_detHists = create_sumsecHist(my_secHists, hcolor=True, these_bins=bins, k=thiscase)
		iso_hist = create_iso_hist(my_detHists, eventfile[i_file], these_bins=bins, k=thiscase)
		iso_sumhist = create_sumdetHist(my_detHists, these_bins=bins, k=thiscase)
		all_isohist.append(iso_sumhist)

	all_mat_hist = create_material_hist(all_isohist, background)
	sum_all_mat_hist = create_summaterial_hist(all_isohist, background, these_bins=bins)
	if savedir:
		save_material_hists(all_mat_hist, sum_all_mat_hist, background, save_dir=savedir)
	else:
		save_material_hists(all_mat_hist, sum_all_mat_hist, background)

	if returns:
		return sum_all_mat_hist


def compare_versions(ev_old, ev_new, scale, x_range, background, background_2):
	for i_file in range(len(ev_old)):
		data = read_evData(ev_old[i_file])
		sim = ev_old[i_file].split('/')
		sim = sim[len(sim)-1]
		sim = sim.split('entries_')
		sim = sim[len(sim)-1]
		sim = sim.split('.')
		sim = sim[0]


		thiscase = background+'_'+sim
		thiscase_2 = background_2+'_'+sim

		my_secHists = create_sector_hists(data, scale[i_file], k=thiscase)
		my_detHists = create_det_hists(my_secHists, k=thiscase)
				
		data_2 = read_evData(ev_new[i_file])
		my_secHists_2 = create_sector_hists(data_2, scale[i_file], k=thiscase_2)
		my_detHists_2 = create_det_hists(my_secHists_2, k=thiscase_2)
		make_the_difference(my_detHists, my_detHists_2, thiscase, thiscase_2, sim)
