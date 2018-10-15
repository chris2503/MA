import numpy as np
from sim_events_3 import *
from mymain import preparation

eventfiles, background, scale, x_range = preparation()

eventfile_coating_noG_noZ = [
						'../rootfiles_noguard_06_09_2018/edep_entries/entries_40K_glyptal.txt',
						'../rootfiles_noguard_06_09_2018/edep_entries/entries_40K_epoxy.txt',
						'../rootfiles_noguard_06_09_2018/edep_entries/entries_232Th_glyptal.txt',
						'../rootfiles_noguard_06_09_2018/edep_entries/entries_232Th_epoxy.txt',
						'../rootfiles_noguard_06_09_2018/edep_entries/entries_238U_glyptal.txt',
						'../rootfiles_noguard_06_09_2018/edep_entries/entries_238U_epoxy.txt']
	
#eventfile_czt_noG_noZ = [
#					'../rootfiles_noguard_06_09_2018/edep_entries/entries_114Cd.txt',
#					'../rootfiles_noguard_06_09_2018/edep_entries/entries_116Cd.txt',
#					'../rootfiles_noguard_06_09_2018/edep_entries/entries_70Zn.txt',
#					'../rootfiles_noguard_06_09_2018/edep_entries/entries_128Te.txt',
#					'../rootfiles_noguard_06_09_2018/edep_entries/entries_130Te.txt']
	
eventfile_POM_noG_noZ = [
						'../rootfiles_noguard_06_09_2018/edep_entries/entries_232Th_dplate.txt',
						'../rootfiles_noguard_06_09_2018/edep_entries/entries_238U_dplate.txt',
						]

eventfile_PA_noG_noZ = [
						'../rootfiles_noguard_06_09_2018/edep_entries/entries_232Th_dscrew.txt',
						'../rootfiles_noguard_06_09_2018/edep_entries/entries_238U_dscrew.txt'
						]

###############################################################################################

eventfile_coating_noG_wZ = [
						'../rootfiles_zcut_only_06_09_2018/edep_entries/entries_40K_glyptal.txt',
						'../rootfiles_zcut_only_06_09_2018/edep_entries/entries_40K_epoxy.txt',
						'../rootfiles_zcut_only_06_09_2018/edep_entries/entries_232Th_glyptal.txt',
						'../rootfiles_zcut_only_06_09_2018/edep_entries/entries_232Th_epoxy.txt',
						'../rootfiles_zcut_only_06_09_2018/edep_entries/entries_238U_glyptal.txt',
						'../rootfiles_zcut_only_06_09_2018/edep_entries/entries_238U_epoxy.txt']
	
#eventfile_czt_noG_wZ = [
#					'../rootfiles_zcut_only_06_09_2018/edep_entries/entries_114Cd.txt',
#					'../rootfiles_zcut_only_06_09_2018/edep_entries/entries_116Cd.txt',
#					'../rootfiles_zcut_only_06_09_2018/edep_entries/entries_70Zn.txt',
#					'../rootfiles_zcut_only_06_09_2018/edep_entries/entries_128Te.txt',
#					'../rootfiles_zcut_only_06_09_2018/edep_entries/entries_130Te.txt']
	
eventfile_POM_noG_wZ = [
						'../rootfiles_zcut_only_06_09_2018/edep_entries/entries_232Th_dplate.txt',
						'../rootfiles_zcut_only_06_09_2018/edep_entries/entries_238U_dplate.txt'
						]
eventfile_PA_noG_wZ = [
						'../rootfiles_zcut_only_06_09_2018/edep_entries/entries_232Th_dscrew.txt',
						'../rootfiles_zcut_only_06_09_2018/edep_entries/entries_238U_dscrew.txt'
						]

###############################################################################################

eventfile_coating_wG_noZ = [
						'../rootfiles_guard_only_06_09_2018/edep_entries/entries_40K_glyptal.txt',
						'../rootfiles_guard_only_06_09_2018/edep_entries/entries_40K_epoxy.txt',
						'../rootfiles_guard_only_06_09_2018/edep_entries/entries_232Th_glyptal.txt',
						'../rootfiles_guard_only_06_09_2018/edep_entries/entries_232Th_epoxy.txt',
						'../rootfiles_guard_only_06_09_2018/edep_entries/entries_238U_glyptal.txt',
						'../rootfiles_guard_only_06_09_2018/edep_entries/entries_238U_epoxy.txt']
	
#eventfile_czt_wG_noZ = [
#					'../rootfiles_guard_only_06_09_2018/edep_entries/entries_114Cd.txt',
#					'../rootfiles_guard_only_06_09_2018/edep_entries/entries_116Cd.txt',
#					'../rootfiles_guard_only_06_09_2018/edep_entries/entries_70Zn.txt',
#					'../rootfiles_guard_only_06_09_2018/edep_entries/entries_128Te.txt',
#					'../rootfiles_guard_only_06_09_2018/edep_entries/entries_130Te.txt']
	
eventfile_POM_wG_noZ = [
						'../rootfiles_guard_only_06_09_2018/edep_entries/entries_232Th_dplate.txt',
						'../rootfiles_guard_only_06_09_2018/edep_entries/entries_238U_dplate.txt'
						]
eventfile_PA_wG_noZ = [
						'../rootfiles_guard_only_06_09_2018/edep_entries/entries_232Th_dscrew.txt',
						'../rootfiles_guard_only_06_09_2018/edep_entries/entries_238U_dscrew.txt'
						]

###############################################################################################

eventfile_coating_wG_wZ = [
						'../rootfiles_guard+zcut_06_09_2018/edep_entries/entries_40K_glyptal.txt',
						'../rootfiles_guard+zcut_06_09_2018/edep_entries/entries_40K_epoxy.txt',
						'../rootfiles_guard+zcut_06_09_2018/edep_entries/entries_232Th_glyptal.txt',
						'../rootfiles_guard+zcut_06_09_2018/edep_entries/entries_232Th_epoxy.txt',
						'../rootfiles_guard+zcut_06_09_2018/edep_entries/entries_238U_glyptal.txt',
						'../rootfiles_guard+zcut_06_09_2018/edep_entries/entries_238U_epoxy.txt']
	
#eventfile_czt_wG_wZ = [
#					'../rootfiles_guard+zcut_06_09_2018/edep_entries/entries_114Cd.txt',
#					'../rootfiles_guard+zcut_06_09_2018/edep_entries/entries_116Cd.txt',
#					'../rootfiles_guard+zcut_06_09_2018/edep_entries/entries_70Zn.txt',
#					'../rootfiles_guard+zcut_06_09_2018/edep_entries/entries_128Te.txt',
#					'../rootfiles_guard+zcut_06_09_2018/edep_entries/entries_130Te.txt']
	
eventfile_POM_wG_wZ = [
						'../rootfiles_guard+zcut_06_09_2018/edep_entries/entries_232Th_dplate.txt',
						'../rootfiles_guard+zcut_06_09_2018/edep_entries/entries_238U_dplate.txt'

						]
eventfile_PA_wG_wZ = [
						'../rootfiles_guard+zcut_06_09_2018/edep_entries/entries_232Th_dscrew.txt',
						'../rootfiles_guard+zcut_06_09_2018/edep_entries/entries_238U_dscrew.txt'
						]

################################
# no plot stuff
################################
bins = 250
case = ['nG_nZ', 'wG_nZ', 'nG_wZ', 'wG_WZ']
scale_POM = [scale[0][0], scale[0][1]]
scale_PA = [scale[0][2], scale[0][3]]
# no guard and no z-cut
all_isohist_10 = []
all_isohist_11 = []
all_isohist_12 = []

for i_file in range(len(eventfile_POM_noG_noZ)):
	data = read_evData(eventfile_POM_noG_noZ[i_file])
	thiscase = case[0] + str(i_file) + '_POM'
	my_secHists_10 = create_sector_hists(data, scale_POM[i_file], these_bins=bins, k=thiscase)
	my_detHists_10 = create_det_hists(my_secHists_10, these_bins=bins, k=thiscase)
	iso_hist_10 = create_iso_hist(my_detHists_10, eventfile_POM_noG_noZ[i_file], these_bins=bins, k=thiscase)
	all_isohist_10.append(iso_hist_10)
mathist_10 = create_summaterial_hist(all_isohist_10, thiscase, these_bins=bins)

for i_file in range(len(eventfile_PA_noG_noZ)):
	data = read_evData(eventfile_PA_noG_noZ[i_file])
	thiscase = case[0] + str(i_file) + '_PA'
	my_secHists_11 = create_sector_hists(data, scale_PA[i_file], these_bins=bins, k=thiscase)
	my_detHists_11 = create_det_hists(my_secHists_11, these_bins=bins, k=thiscase)
	iso_hist_11 = create_iso_hist(my_detHists_11, eventfile_PA_noG_noZ[i_file], these_bins=bins, k=thiscase)
	all_isohist_11.append(iso_hist_11)
mathist_11 = create_summaterial_hist(all_isohist_11, thiscase, these_bins=bins)

for i_file in range(len(eventfile_coating_noG_noZ)):
	data = read_evData(eventfile_coating_noG_noZ[i_file])
	thiscase = case[0] + str(i_file)  + '_coating'
	my_secHists_12 = create_sector_hists(data, scale[1][i_file], these_bins=bins, k=thiscase)
	my_detHists_12 = create_det_hists(my_secHists_12, these_bins=bins, k=thiscase)
	iso_hist_12 = create_iso_hist(my_detHists_12, eventfile_coating_noG_noZ[i_file], these_bins=bins, k=thiscase)
	all_isohist_12.append(iso_hist_12)
mathist_12 = create_summaterial_hist(all_isohist_12, thiscase, these_bins=bins)


# with guard and no z-cut
all_isohist_20 = []
all_isohist_21 = []
all_isohist_22 = []

for i_file in range(len(eventfile_POM_wG_noZ)):
	data = read_evData(eventfile_POM_wG_noZ[i_file])
	thiscase = case[1] + str(i_file) + '_POM'
	my_secHists_20 = create_sector_hists(data, scale_POM[i_file], these_bins=bins, k=thiscase)
	my_detHists_20 = create_det_hists(my_secHists_20, these_bins=bins, k=thiscase)
	iso_hist_20 = create_iso_hist(my_detHists_20, eventfile_POM_wG_noZ[i_file], these_bins=bins, k=thiscase)
	all_isohist_20.append(iso_hist_20)
mathist_20 = create_summaterial_hist(all_isohist_20, thiscase, these_bins=bins)

for i_file in range(len(eventfile_PA_wG_noZ)):
	data = read_evData(eventfile_PA_wG_noZ[i_file])
	thiscase = case[1] + str(i_file)  + '_PA'
	my_secHists_21 = create_sector_hists(data, scale_PA[i_file], these_bins=bins, k=thiscase)
	my_detHists_21 = create_det_hists(my_secHists_21, these_bins=bins, k=thiscase)
	iso_hist_21 = create_iso_hist(my_detHists_21, eventfile_PA_wG_noZ[i_file], these_bins=bins, k=thiscase)
	all_isohist_21.append(iso_hist_21)
mathist_21 = create_summaterial_hist(all_isohist_21, thiscase, these_bins=bins)

for i_file in range(len(eventfile_coating_wG_noZ)):
	data = read_evData(eventfile_coating_wG_noZ[i_file])
	thiscase = case[1] + str(i_file)  + '_coating'
	my_secHists_22 = create_sector_hists(data, scale[1][i_file], these_bins=bins, k=thiscase)
	my_detHists_22 = create_det_hists(my_secHists_21, these_bins=bins, k=thiscase)
	iso_hist_22 = create_iso_hist(my_detHists_22, eventfile_coating_wG_noZ[i_file], these_bins=bins, k=thiscase)
	all_isohist_22.append(iso_hist_22)
mathist_22 = create_summaterial_hist(all_isohist_22, thiscase, these_bins=bins)


# no guard and with z-cut
all_isohist_30 = []
all_isohist_31 = []
all_isohist_32 = []

for i_file in range(len(eventfile_POM_noG_wZ)):
	data = read_evData(eventfile_POM_noG_wZ[i_file])
	thiscase = case[2] + str(i_file) + '_POM'
	my_secHists_30 = create_sector_hists(data, scale_POM[i_file], these_bins=bins, k=thiscase)
	my_detHists_30 = create_det_hists(my_secHists_30, these_bins=bins, k=thiscase)
	iso_hist_30 = create_iso_hist(my_detHists_30, eventfile_POM_noG_wZ[i_file], these_bins=bins, k=thiscase)
	all_isohist_30.append(iso_hist_30)
mathist_30 = create_summaterial_hist(all_isohist_30, thiscase, these_bins=bins)

for i_file in range(len(eventfile_PA_noG_wZ)):
	data = read_evData(eventfile_PA_noG_wZ[i_file])
	thiscase = case[2] + str(i_file)  + '_PA'
	my_secHists_31 = create_sector_hists(data, scale_PA[i_file], these_bins=bins, k=thiscase)
	my_detHists_31 = create_det_hists(my_secHists_31, these_bins=bins, k=thiscase)
	iso_hist_31 = create_iso_hist(my_detHists_31, eventfile_PA_noG_wZ[i_file], these_bins=bins, k=thiscase)
	all_isohist_31.append(iso_hist_31)
mathist_31 = create_summaterial_hist(all_isohist_31, thiscase, these_bins=bins)

for i_file in range(len(eventfile_coating_noG_wZ)):
	data = read_evData(eventfile_coating_noG_wZ[i_file])
	thiscase = case[2] + str(i_file)  + '_coating'
	my_secHists_32 = create_sector_hists(data, scale[1][i_file], these_bins=bins, k=thiscase)
	my_detHists_32 = create_det_hists(my_secHists_32, these_bins=bins, k=thiscase)
	iso_hist_32 = create_iso_hist(my_detHists_32, eventfile_coating_noG_wZ[i_file], these_bins=bins, k=thiscase)
	all_isohist_32.append(iso_hist_32)
mathist_32 = create_summaterial_hist(all_isohist_32, thiscase, these_bins=bins)


# with guard and with z-cut
all_isohist_40 = []
all_isohist_41 = []
all_isohist_42 = []

for i_file in range(len(eventfile_POM_wG_wZ)):
	data = read_evData(eventfile_POM_wG_wZ[i_file])
	thiscase = case[3] + str(i_file) + '_POM'
	my_secHists_40 = create_sector_hists(data, scale_POM[i_file], these_bins=bins, k=thiscase)
	my_detHists_40 = create_det_hists(my_secHists_40, these_bins=bins, k=thiscase)
	iso_hist_40 = create_iso_hist(my_detHists_40, eventfile_POM_wG_wZ[i_file], these_bins=bins, k=thiscase)
	all_isohist_40.append(iso_hist_40)
mathist_40 = create_summaterial_hist(all_isohist_40, thiscase, these_bins=bins)

for i_file in range(len(eventfile_PA_wG_wZ)):
	data = read_evData(eventfile_PA_wG_wZ[i_file])
	thiscase = case[3] + str(i_file)  + '_PA'
	my_secHists_41 = create_sector_hists(data, scale_PA[i_file], these_bins=bins, k=thiscase)
	my_detHists_41 = create_det_hists(my_secHists_41, these_bins=bins, k=thiscase)
	iso_hist_41 = create_iso_hist(my_detHists_41, eventfile_PA_wG_wZ[i_file], these_bins=bins, k=thiscase)
	all_isohist_41.append(iso_hist_41)
mathist_41 = create_summaterial_hist(all_isohist_41, thiscase, these_bins=bins)

for i_file in range(len(eventfile_coating_wG_wZ)):
	data = read_evData(eventfile_coating_wG_wZ[i_file])
	thiscase = case[3] + str(i_file)  + '_coating'
	my_secHists_42 = create_sector_hists(data, scale[1][i_file], these_bins=bins, k=thiscase)
	my_detHists_42 = create_det_hists(my_secHists_42, these_bins=bins, k=thiscase)
	iso_hist_42 = create_iso_hist(my_detHists_42, eventfile_coating_wG_wZ[i_file], these_bins=bins, k=thiscase)
	all_isohist_42.append(iso_hist_42)
mathist_42 = create_summaterial_hist(all_isohist_42, thiscase, these_bins=bins)

mathists_1 = create_allsumHist([mathist_10, mathist_11, mathist_12], '1', these_bins=bins)
mathists_2 = create_allsumHist([mathist_20, mathist_21, mathist_22], '2', these_bins=bins)
mathists_3 = create_allsumHist([mathist_30, mathist_31, mathist_32], '3', these_bins=bins)
mathists_4 = create_allsumHist([mathist_40, mathist_41, mathist_42], '4', these_bins=bins)

all_mathists = [mathists_1, mathists_2, mathists_3, mathists_4]
create_comp_plots(all_mathists)
