import numpy as np
from mymain import compare_versions, preparation

eventfiles, background, scale, x_range = preparation()

eventfile_coating_1 = ['../rootfiles_old_geant_29_08_2018/edep_entries/entries_40K_glyptal.txt',
							'../rootfiles_old_geant_29_08_2018/edep_entries/entries_40K_epoxy.txt',
							'../rootfiles_old_geant_29_08_2018/edep_entries/entries_232Th_glyptal.txt',
							'../rootfiles_old_geant_29_08_2018/edep_entries/entries_232Th_epoxy.txt',
							'../rootfiles_old_geant_29_08_2018/edep_entries/entries_238U_glyptal.txt',
							'../rootfiles_old_geant_29_08_2018/edep_entries/entries_238U_epoxy.txt']
	
eventfile_czt_1 = ['../rootfiles_old_geant_29_08_2018/edep_entries/entries_114Cd.txt',
						'../rootfiles_old_geant_29_08_2018/edep_entries/entries_116Cd.txt',
						'../rootfiles_old_geant_29_08_2018/edep_entries/entries_70Zn.txt',
						'../rootfiles_old_geant_29_08_2018/edep_entries/entries_128Te.txt',
						'../rootfiles_old_geant_29_08_2018/edep_entries/entries_130Te.txt']
	
eventfile_plastic_1 = ['../rootfiles_old_geant_29_08_2018/edep_entries/entries_232Th_dplate.txt',
									'../rootfiles_old_geant_29_08_2018/edep_entries/entries_238U_dplate.txt',
									'../rootfiles_old_geant_29_08_2018/edep_entries/entries_232Th_dscrew.txt',
									'../rootfiles_old_geant_29_08_2018/edep_entries/entries_238U_dscrew.txt']

eventfile_coating_2 = ['../rootfiles_geant4.10.4/edep_entries/entries_40K_glyptal.txt',
							'../rootfiles_geant4.10.4/edep_entries/entries_40K_epoxy.txt',
							'../rootfiles_geant4.10.4/edep_entries/entries_232Th_glyptal.txt',
							'../rootfiles_geant4.10.4/edep_entries/entries_232Th_epoxy.txt',
							'../rootfiles_geant4.10.4/edep_entries/entries_238U_glyptal.txt',
							'../rootfiles_geant4.10.4/edep_entries/entries_238U_epoxy.txt']
	
eventfile_czt_2 = ['../rootfiles_geant4.10.4/edep_entries/entries_114Cd.txt',
						'../rootfiles_geant4.10.4/edep_entries/entries_116Cd.txt',
						'../rootfiles_geant4.10.4/edep_entries/entries_70Zn.txt',
						'../rootfiles_geant4.10.4/edep_entries/entries_128Te.txt',
						'../rootfiles_geant4.10.4/edep_entries/entries_130Te.txt']
	
eventfile_plastic_2 = ['../rootfiles_geant4.10.4/edep_entries/entries_232Th_dplate.txt',
									'../rootfiles_geant4.10.4/edep_entries/entries_238U_dplate.txt',
									'../rootfiles_geant4.10.4/edep_entries/entries_232Th_dscrew.txt',
									'../rootfiles_geant4.10.4/edep_entries/entries_238U_dscrew.txt']

background_2=background
for i in range(len(background_2)):
	background_2[i]=background_2[i]+"_2"

scale_czt = scale[0]
x_range_czt = x_range[0]

scale_coating = scale[1]
x_range_coating = x_range[1]

scale_plastic = scale[2]
x_range_plastic = x_range[2]

background_plastic = background[0]
background_coating = background[1]
background_czt = background[2]

background_2_plastic = background_2[0]
background_2_coating = background_2[1]
background_2_czt = background_2[2]

compare_versions(eventfile_czt_1, eventfile_czt_2, scale_czt, x_range_czt, background_czt, background_2_czt)
