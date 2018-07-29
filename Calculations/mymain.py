import numpy as np
from sim_events_3 import *
##############################################
# main
##############################################
kprint=None

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



x_range_coating = 1e4
x_range_czt = 3e3

######################
# Delrin layer and screws
######################



################
#### Coating
################
m_det_or = np.array([34.21, 35.15, 34.63, 33.91, 35.50, 35.50, 35.50, 35.50, 35.50]) * 1e-3			#kg
m_detpaint_glyptal_or = np.array([0.06, 0.07, 0.06, 0.07]) * 1e-3									#kg
m_detpaint_epoxy_or = np.array([0.14, 0.13, 0.14, 0.12, 0.16]) * 1e-3								#kg

m_det = np.mean(m_det_or)
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

years = 60*60*24*365
n_det_glyptal = 4
n_det_epoxy = 5
n_chain = [1, 10, 14]

datafile = './Lists/activities_detpaint.txt'
data = read_File(datafile)
iso_list = data[:,0]								# get isotope
sp_ac_epoxy = convert_str2num(data[:,1])				# get norming factors, number convertion necessary
sp_ac_glyptal = convert_str2num(data[:,2])				# get norming factors, number convertion necessary

N_norm_coating = []
scale_coating = []
for i in range(len(iso_list)):
	N_norm_coating.append(sp_ac_glyptal[i] * m_detpaint_glyptal * years * n_chain[i] * 1/(n_det_glyptal* m_det))
	N_norm_coating.append(sp_ac_epoxy[i] * m_detpaint_epoxy * years * n_chain[i] * 1/(n_det_epoxy* m_det))

N_simEv = 1e6
scale_czt = []
for i in range(len(N_norm_coating)):
	scale_coating.append(1/4 * N_norm_coating[i]/N_simEv)

###########################
#### CZT
###########################
datafile = './calc_solutions/calculated_events.txt'
data = read_File(datafile)

iso_list = data[:,0]							# get isotope
N_norm_czt = convert_str2num(data[:,4])			# get norming factors, number convertion necessary
N_simEv = 1e6									# 1 Mio simulated Events

for i in range(len(N_norm_czt)):
	scale_czt.append(1/36 * N_norm_czt[i] / N_simEv) # It was calculated for 9 detectors


def my_main(eventfile, scale, x_range, background, returns = None):
	# creating plots
	# 1) single plots
#	for i_file in range(len(eventfile)):
#		data = read_evData(eventfile[i_file])
#		my_secHists = create_sector_hists(data, scale[i_file])
#		save_single_sector_hists(my_secHists, eventfile[i_file], x_range)
#		delete_all_sectorHists(my_secHists)
#
#		my_secHists = create_sector_hists(data, scale[i_file])
#		my_detHists = create_det_hists(my_secHists)
#		my_sumhists = create_sumsecHist(my_secHists)
#		save_single_det_hists(my_detHists, my_sumhists, eventfile[i_file], x_range)
#		delete_all_sumdetHists(my_sumhists)
#		delete_all_detHists(my_detHists)
#		delete_all_sectorHists(my_secHists)
#
#		my_secHists = create_sector_hists(data, scale[i_file])
#		my_detHists = create_sumsecHist(my_secHists, hcolor=True)
#		iso_hist = create_iso_hist(my_detHists, eventfile[i_file])
#		iso_sumhist = create_sumdetHist(my_detHists)
#		save_single_iso_hists(iso_hist, iso_sumhist, eventfile[i_file], x_range)
#		delete_iso_Hist(iso_hist, eventfile[i])
#		delete_all_detHists(my_detHists)
#		delete_all_sectorHists(my_secHists)
#		delete_iso_sumHist(iso_sumhist)


	#for i_file in range(len(eventfile)):
	#	delete_iso_Hist(iso_hist, eventfile[i], k=i_file)
	#	delete_all_detHists(my_detHists, k=i_file)
	#	delete_all_sectorHists(my_secHists, k=i_file)
	#	delete_iso_sumHist(iso_sumhist, k=i_file)


#	for i_file in range(len(eventfile)):
#		data = read_evData(eventfile[i_file])
#		my_depHists = create_dep_secHist(data, scale[i_file])
#		save_dep_sec_hists(my_depHists, eventfile[i_file])
#		delete_all_detHists(my_depHists)
#
#		my_depHist = create_dep_detHist(data, eventfile[i_file], scale[i_file])
#		save_dep_det_hists(my_depHist, eventfile[i_file])
#		delete_iso_Hist(my_depHist, eventfile[i])
#
#		my_depHists = create_dep_secHist(data, scale[i_file])
#		save_dep_sec_hists_2(my_depHists, eventfile[i_file], None)
#		delete_all_detHists(my_depHists)
#
#	for i_file in range(len(eventfile)):
#		data = read_evData(eventfile[i_file])
#		my_depHists_heat = create_dep_secHist(data, scale[i_file])
#		save_dep_heatmap(my_depHists_heat, eventfile[i_file])
#		delete_all_detHists(my_depHists_heat)
#

	# 2) Values
	for i_file in range(len(eventfile)):
		thiscase = background+str(i_file)
		data = read_evData(eventfile[i_file])

		my_secHists , contrib_at116Cd, contrib_at116Cd_err, contrib_at130Te, contrib_at130Te_err = create_sector_hists(data, scale[i_file], k=thiscase, Q_val_returns=True)
		print(contrib_at116Cd.size)
		save_single_sector_hists(my_secHists, eventfile[i_file], x_range)
		delete_all_sectorHists(my_secHists)


	# 3) combined plots
	all_isohist =[]
	for i_file in range(len(eventfile)):
		thiscase = background+str(i_file)
		data = read_evData(eventfile[i_file])

		my_secHists = create_sector_hists(data, scale[i_file], k=thiscase)
		my_detHists = create_sumsecHist(my_secHists, hcolor=True, k=thiscase)
		iso_hist = create_iso_hist(my_detHists, eventfile[i_file], k=thiscase)
		iso_sumhist = create_sumdetHist(my_detHists, k=thiscase)
		all_isohist.append(iso_sumhist)

	all_mat_hist = create_material_hist(all_isohist, background)
	sum_all_mat_hist = create_summaterial_hist(all_isohist, background)
	#save_material_hists(all_mat_hist, sum_all_mat_hist, background)




	if returns:
		return sum_all_mat_hist





background = ['coating', 'czt']
#my_main(eventfile_coating, scale_coating, x_range_coating, background[0])
#my_main(eventfile_czt, scale_czt, x_range_czt, background[1])
all_hists = []
all_hists.append(my_main(eventfile_coating, scale_coating, x_range_coating, background[0], returns=True))
all_hists.append(my_main(eventfile_czt, scale_czt, x_range_czt, background[1], returns=True))

hists = create_sumHist(all_hists, background)
all_sumhists = create_allsumHist(all_hists, background)
save_sumHist(hists, all_sumhists, background)
