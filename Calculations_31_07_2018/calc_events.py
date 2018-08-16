import numpy as np
import scipy.constants as const
import re
import os
import sys

#sys.path.append('../')
#sys.path.append(os.path.abspath('../'))
#print (sys.path)
from Methods.my_methods import read_File, convert_str2num, write_txtfile, write_detailed_txtfile


# make list of elements from given isotopes
def list_elements(iso) :
	counter = 0
	for isotopes in iso_list:
		temp = re.search("([a-zA-Z]+)", isotopes)
		if counter == 0:									# list element of 1. isotope
			element_list = [temp.group(0)]
			counter = counter + 1
		else:
			if any(temp.group(0) in counter_2 for counter_2 in element_list):	# if isotope is already listed
				continue														# reject equal elements 
			else:
				element_list = np.append([element_list],[temp.group(0)])
	return element_list

# evaluates total percentage of the istopes' percentage of the element in e.g. CZT-detector
# for other elements, please add info
def calc_tot_perc(iso_list, el_list, el_fractions, perc):
	counter = 0
	tot_perc_list = []
	for isotope in iso_list:
		temp = re.search("([a-zA-Z]+)", isotope)
		if temp.group(0) == 'Cd':
			this_perc = el_fractions[0] * perc[counter]
		elif temp.group(0) == 'Zn':
			this_perc = el_fractions[1] * perc[counter]
		elif temp.group(0) == 'Te':
			this_perc = el_fractions[2] * perc[counter]
		else:
			if any(temp.group(0) in counter_2 for counter_2 in el_list):
				print('Warning: No volume percentage of element "%s" listed!' %(temp.group(0)))
			else:
				print('Warning: No information about "%s"!' %(isotope))
			this_perc = 0

		tot_perc_list = np.append([tot_perc_list], this_perc)
		counter = counter + 1
	return np.float32(tot_perc_list)

# calculate masses from the detector components
def calc_det_compMass(rho, V):			#	V: volume, rho: density
	mass = []
	for counter in range(len(V)):
		mass= np.append(mass, rho[counter] * V[counter])	# 9 detectors
		counter = counter + 1  
	return np.float32(mass)

# calculate the total mass of the isotope in the detector components
def calc_tot_mass(mass_array, tot_perc, iso_list, el_fractions, el_molmass, el_list):		# m: masses of components of the detector
	mass = []
	counter = 0
	M_Det = np.sum(el_fractions*el_molmass)
	for isotope in iso_list:
		temp = re.search("([a-zA-Z]+)", isotope)
		if temp.group(0) == 'Cd':
			mass = np.append(mass, mass_array[0] * tot_perc[counter] *el_molmass[0]/M_Det)
		elif temp.group(0) == 'Zn':
			mass = np.append(mass, mass_array[0] * tot_perc[counter] *el_molmass[1]/M_Det)
		elif temp.group(0) == 'Te':
			mass = np.append(mass, mass_array[0] * tot_perc[counter] *el_molmass[2]/M_Det)
		else: 
			mass = np.append(mass, mass_array[1] * tot_perc[counter])
		counter = counter + 1
	return np.float32(mass)

# calculate the activity of the selected istotope
def calc_spec_activity(t_12, M):
	N_A = const.N_A
	return np.float32(np.log(2)/t_12 * N_A/M) 	# t_12: half life, N_A: Avogadro constant, M: molar mass

# calculate the real number of events of the selected istotope for a period of time (e.g 1year for scaleing ;) )
def calc_events(t, m, a):	
	return np.float32(t* m * a) 			# t: time, A: activity, m: mass




#######################
# run program
#######################
# defining constants
time = 60*60*24*365				#	[s] -> 1 year
N_A = const.N_A					#	Avogadro constant

# volume of the detector:
V_det = 9*(20.4 * 10**(-3))**2 * 15*10**(-3)		# [m^3]		# 9 detectors !
# volume of the paint around the detector (whole detector surface exept from bottom)
V_comp = [V_det]

# densities
rho_CZT = 5.78 * 1e-3/1e-6						# [g/cm^3] -> [kg/m^3]
rho_varnish = 1.18 * 1e-3/1e-6					# [g/cm^3] -> [kg/m^3],	Polymethyl_methacrylate
rho_comp = [rho_CZT, rho_varnish]


print('V_CZT = ', (V_det)) 
print('m_CZT = ', (rho_CZT * V_det))
print()



# fractions of Cd, Zn and Te in CZT in this order:
CdZnTe_fractions = np.array([0.45, 0.05, 0.50])
CdZnTe_molmass = np.array([112.41, 65.38, 127.6]) *1e-3 # http://www.periodensystem.info/periodensystem/


# load isotope properties from list
data = read_File("./Lists/list_isotope_properties_czt.txt")
iso_list = data[:,0]							# these arrays are socalled 'ndarrays'
t_12 = convert_str2num(data[:,1])	* time		# [s]	t_12: half life
iso_perc = convert_str2num(data[:,2])			# [unitless]	iso_perc: percentage of the isotope's abundance in the element 
mol_mass = convert_str2num(data[:,3]) * 1e-3	# [g/mole] -> [kg/mole]	mol_mass: molar mass

print('Isotopes:\n %s' %(iso_list))
print('')


element_list = list_elements(iso_list)
print('List of Elements: \n %s' %(element_list))
print('')

tot_iso_perc = calc_tot_perc(iso_list, element_list, CdZnTe_fractions, iso_perc)
print('Total fraction of the isotopes in the chemical compound: \n %s' %(tot_iso_perc))
print('')

mass_det_comp = calc_det_compMass(rho_comp, V_comp)
print('Masses of the components: %s' %(mass_det_comp))
print('')

tot_iso_mass = calc_tot_mass(mass_det_comp, tot_iso_perc, iso_list, CdZnTe_fractions, CdZnTe_molmass, element_list)
print('Masses of the isotopes in the components: \n %s' %(tot_iso_mass))
print('')

spec_ac = calc_spec_activity(t_12, mol_mass)
print('Special activities of the istotopes: \n %s' %(spec_ac))
print('')

# scale to 1year
N_iso_ev = calc_events(time, tot_iso_mass, spec_ac)
print('Total number of events after 1 year: \n %s' %(N_iso_ev))
print('')

# TOTAL NUMBER OF EVENTS Scaled TO detector mass [kg] for 1 year
N_iso_scaled = N_iso_ev * 1/(mass_det_comp[0])					#	9 detectors
print('Total number of events after 1 year scaled to detector mass: \n %s' %(N_iso_scaled))
print('')	


### export data to a txt-file
# for saving data to a table, data must be transposed
new_data = np.array([iso_list, tot_iso_perc, tot_iso_mass, spec_ac, N_iso_scaled])
#write_txtfile(np.transpose(new_data), '../calc_solutions/', 'calculated_events.txt')
descriptions = ['List of calculated variables',
				'p: total fraction of the isotope in the detector volume',
				'm: total mass fraction of the isotope in the detector volume',
				'a: specific activity of the isotope',
				'N_scaled: REAL NUMBER of events scaled per 1kg detector mass AND 1year']
var_names = ['isotope', 'p', 'm / [kg]', 'a / [Bq/kg]', 'N_scaled [1/(kg*yr)]']
write_detailed_txtfile(np.transpose(new_data), var_names, descriptions, './calc_solutions/', 'calculated_events.txt')
