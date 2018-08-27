import numpy as np
from sim_events_3 import *
from Methods.my_methods import read_File, convert_str2num, write_txtfile, write_detailed_txtfile
import uncertainties
import uncertainties.unumpy as unp
from uncertainties.unumpy import(nominal_values as noms, std_devs as stds)
from mymain import my_main_1

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

eventfile_plastic = ['../rootfiles/edep_entries/entries_232Th_dplate.txt',
								'../rootfiles/edep_entries/entries_238U_dplate.txt',
								'../rootfiles/edep_entries/entries_232Th_dscrew.txt',
								'../rootfiles/edep_entries/entries_238U_dscrew.txt']


x_range_plastic = 1e4
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

iso_list = this_data[:,0]							# get isotope
N_norm_czt = convert_str2num(this_data[:,1])			# get norming factors, number convertion necessary
N_simEv = 1e6									# 1 Mio simulated Events
scale_czt = []								

for i in range(len(N_norm_czt)):
	scale_czt.append(1/4 * N_norm_czt[i] / N_simEv) # It was calculated for 9 detectors

#background = ['plastic', 'coating', 'czt']
#eventfile = [eventfile_plastic, eventfile_coating, eventfile_czt]
#scale = [scale_plastic, scale_coating, scale_czt]
#x_range = [x_range_coating, x_range_coating, x_range_czt]

my_main_1(eventfile_czt, scale_czt, x_range_czt)
