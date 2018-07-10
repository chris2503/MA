import numpy as np
import uncertainties
import uncertainties.unumpy as unp
from uncertainties.unumpy import(nominal_values as noms, std_devs as stds)
import scipy.constants as const
import re
import os
import sys
from time import sleep

import ROOT

sys.path.append(os.path.abspath('./Methods/'))
#print (sys.path)
from Methods.my_methods import read_File, convert_str2num, transform2latex_tab_2

def f(x):
	x_D = 10.2*1e-3				# m													
	z_D = 16.0*1e-3				# m										
	return 4*x*x*x + 4*x*x*(x_D+z_D) + (x_D*x_D + 4*x_D * z_D)

m_det_or = np.array([34.21, 35.15, 34.63, 33.91, 35.50, 35.50, 35.50, 35.50, 35.50]) * 1e-3	#kg
m_detpaint_glyptal_or = np.array([0.06, 0.07, 0.06, 0.07]) * 1e-3							#kg
m_detpaint_epoxy_or= np.array([0.14, 0.13, 0.14, 0.12]) * 1e-3								#kg

m_det = np.mean(m_det_or)
m_detpaint_epoxy = np.mean(m_detpaint_epoxy_or)
m_detpaint_glyptal = np.mean(m_detpaint_glyptal_or)

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

print("\nepoxy-Lackdicke: %f" %(d_epoxy))
print("glyptal-Lackdicke: %f" %(d_glyptal))
print(10*"#")

def pause():
	print()
	wait = input("PRESS ENTER TO CONTINUE")
	print()


def sigma_E(E):
	a = 1
	b = 1
	c = 1
	E = np.sqrt(a**2 + (E/1000)*b**2 + (E/1000)**2 * c**2)	# E: energy in ! MeV !; a, b, c: parameters
	return E
	# source: https://wiki.cobra-experiment.org/Alternative_to_standard_linear_resolution_function



def background_coating(eventfile, datafile, savename):
	################################
	# Unstacked Histogram
	################################
	years = 60*60*24*365
	n_det_glyptal = 5
	n_det_epoxy = 4
	n_chain = [1, 10, 14]


	data = read_File(datafile)
	iso_list = data[:,0]								# get isotope
	sp_ac_epoxy = convert_str2num(data[:,1])				# get norming factors, number convertion necessary
	sp_ac_glyptal = convert_str2num(data[:,2])				# get norming factors, number convertion necessary

	N_norm = []
	for i in range(len(iso_list)):
		N_norm.append(sp_ac_glyptal[i] * m_detpaint_glyptal * years * 1/(n_det_glyptal* m_det))
		N_norm.append(sp_ac_epoxy[i] * m_detpaint_epoxy * years * 1/(n_det_epoxy* m_det))

	#print('Data:\n%s \n' %(data))
	#print('Isotopes:\n%s \n' %(iso_list))
	print('Number of events per kilogram and year: \n %s\n' %(N_norm))
	
	all_events = []
	particle_id = []
	print('Reading files:')
	for i in range(len(eventfile)):
		print('\t"%s"' %(eventfile[i]))
	# reading all eventdata and: MeV -> keV
		ev_data = np.genfromtxt(eventfile[i], unpack=True)
		events = 1000*ev_data[0] 				# keV
		particle_id_temp = ev_data[2]
		particle_id.append(particle_id_temp) 
		all_events.append(events)
	print()

	bins = 1000
	E_range = 10000

	# entries at Qvalues:
	Q_116Cd = 2813
	Q_130Te = 2527
	N_at116Cd = []
	N_at130Te = []
	N_at116Cd_particles = []
	N_at130Te_particles = []
	stat_err_at116Cd = []
	stat_err_at130Te = []

	atQ_116Cd = np.int(np.around(Q_116Cd*bins/E_range, decimals=0))
	#print(atQ_116Cd)
	atQ_130Te = np.int(np.around(Q_130Te*bins/E_range, decimals=0))

	canvas = ROOT.TCanvas('canv', 'Histogramm')
	hists = ROOT.THStack("hists", "")
	rand_gen = ROOT.TRandom3()

	# ROOT dislikes numpy-objects so it's necessary to cast to python-objects
	# all_events[i]-entries: np.ndarray is casted to list, 
	# all_events[i][j]-entries: np.float64 is casted to float (automatically)
	for i in range(len(all_events)):
		all_events[i] = all_events[i].tolist()
	for i in range(len(particle_id)):
		particle_id[i] = particle_id[i].tolist()


	# Consider energy resolution of the detector
	#resolution = all_events
	#for i in range(len(all_events)):
	#	for j in range(len(all_events[i])):
	#		resolution[i][j] = rand_gen.Gaus(all_events[i][j], sigma_E(all_events[i][j]))

	
	# creating every single histogram
	hist_1 = ROOT.TH1F('hist_1', '', bins, 0, E_range)
	bin_width = hist_1.GetBinWidth(1)
	for j in range(len(all_events[0])):
		hist_1.Fill(all_events[0][j])
	#hist_1.SetFillColor(2)
	#hist_1.SetFillStyle(1001)
	hist_1.SetLineColor(2)
	hist_1.SetLineWidth(2)
	# Calculating statistical errors
	stat_err_at116Cd_temp = np.sqrt(hist_1.GetBinContent(atQ_116Cd))
	stat_err_at130Te_temp = np.sqrt(hist_1.GetBinContent(atQ_130Te))
	#bin_width = hist_1.GetBinWidth(1)

	# Scaling
	scale = n_chain[0] * N_norm[0]/(1e6)		# for 1 million events # divide per bin width for per keV
	stat_err_at116Cd_temp = scale * stat_err_at116Cd_temp
	stat_err_at116Cd.append(stat_err_at116Cd_temp)
	stat_err_at130Te_temp = scale * stat_err_at130Te_temp
	stat_err_at130Te.append(stat_err_at130Te_temp)
	hist_1.Scale(scale)

	hist_2 = ROOT.TH1F('hist_2', '', bins, 0, E_range)
	for j in range(len(all_events[1])):
		hist_2.Fill(all_events[1][j])
	# Statistical errors
	stat_err_at116Cd_temp = np.sqrt(hist_2.GetBinContent(atQ_116Cd))
	stat_err_at130Te_temp = np.sqrt(hist_2.GetBinContent(atQ_130Te))
	hist_2.SetLineColor(3)
	hist_2.SetLineWidth(2)
	# Scaling
	scale = n_chain[0] * N_norm[1]/(1e6)		# for 1 million events # divide per bin width for per keV
	stat_err_at116Cd_temp = scale * stat_err_at116Cd_temp
	stat_err_at116Cd.append(stat_err_at116Cd_temp)
	stat_err_at130Te_temp = scale * stat_err_at130Te_temp
	stat_err_at130Te.append(stat_err_at130Te_temp)
	hist_2.Scale(scale)

	hist_3 = ROOT.TH1F('hist_3', '', bins, 0, E_range)
	for j in range(len(all_events[2])):
		hist_3.Fill(all_events[2][j])
	stat_err_at116Cd_temp = np.sqrt(hist_3.GetBinContent(atQ_116Cd))
	stat_err_at130Te_temp = np.sqrt(hist_3.GetBinContent(atQ_130Te))
	hist_3.SetLineColor(4)
	hist_3.SetLineWidth(2)
	# Scaling
	scale = n_chain[1] * N_norm[2]/(1e6)		# for 1 million events # divide per bin width for per keV
	stat_err_at116Cd_temp = scale * stat_err_at116Cd_temp
	stat_err_at116Cd.append(stat_err_at116Cd_temp)
	stat_err_at130Te_temp = scale * stat_err_at130Te_temp
	stat_err_at130Te.append(stat_err_at130Te_temp)
	hist_3.Scale(scale)

	hist_4 = ROOT.TH1F('hist_4', '', bins, 0, E_range)
	for j in range(len(all_events[3])):
		hist_4.Fill(all_events[3][j])
	stat_err_at116Cd_temp = np.sqrt(hist_4.GetBinContent(atQ_116Cd))
	stat_err_at130Te_temp = np.sqrt(hist_4.GetBinContent(atQ_130Te))
	hist_4.SetLineColor(5)
	hist_4.SetLineWidth(2)
	# Scaling
	scale = n_chain[1] * N_norm[3]/(1e6)		# for 1 million events # divide per bin width for per keV
	stat_err_at116Cd_temp = scale * stat_err_at116Cd_temp
	stat_err_at116Cd.append(stat_err_at116Cd_temp)
	stat_err_at130Te_temp = scale * stat_err_at130Te_temp
	stat_err_at130Te.append(stat_err_at130Te_temp)
	hist_4.Scale(scale)

	hist_5 = ROOT.TH1F('hist_5', '', bins, 0, E_range)
	for j in range(len(all_events[4])):
		hist_5.Fill(all_events[4][j])
	stat_err_at116Cd_temp = np.sqrt(hist_5.GetBinContent(atQ_116Cd))
	stat_err_at130Te_temp = np.sqrt(hist_5.GetBinContent(atQ_130Te))
	hist_5.SetLineColor(6)
	hist_5.SetLineWidth(2)
	# Scaling
	scale = n_chain[2] *  N_norm[4]/(1e6)		# for 1 million events # divide per bin width for per keV
	stat_err_at116Cd_temp = scale * stat_err_at116Cd_temp
	stat_err_at116Cd.append(stat_err_at116Cd_temp)
	stat_err_at130Te_temp = scale * stat_err_at130Te_temp
	stat_err_at130Te.append(stat_err_at130Te_temp)
	hist_5.Scale(scale)

	hist_6 = ROOT.TH1F('hist_6', '', bins, 0, E_range)
	for j in range(len(all_events[5])):
		hist_6.Fill(all_events[5][j])
	stat_err_at116Cd_temp = np.sqrt(hist_6.GetBinContent(atQ_116Cd))
	stat_err_at130Te_temp = np.sqrt(hist_6.GetBinContent(atQ_130Te))
	hist_6.SetLineColor(7)
	hist_6.SetLineWidth(2)
	# Scaling
	scale = n_chain[2] * N_norm[5]/(1e6)		# for 1 million events # divide per bin width for per keV
	stat_err_at116Cd_temp = scale * stat_err_at116Cd_temp
	stat_err_at116Cd.append(stat_err_at116Cd_temp)
	stat_err_at130Te_temp = scale * stat_err_at130Te_temp
	stat_err_at130Te.append(stat_err_at130Te_temp)
	hist_6.Scale(scale)

	# Getting number of events at Q_values:
	N_at116Cd.append(hist_1.GetBinContent(atQ_116Cd))
	N_at116Cd.append(hist_2.GetBinContent(atQ_116Cd))
	N_at116Cd.append(hist_3.GetBinContent(atQ_116Cd))
	N_at116Cd.append(hist_4.GetBinContent(atQ_116Cd))
	N_at116Cd.append(hist_5.GetBinContent(atQ_116Cd))
	N_at116Cd.append(hist_6.GetBinContent(atQ_116Cd))
	N_at116Cd = unp.uarray(N_at116Cd, stat_err_at116Cd)

	N_at130Te.append(hist_1.GetBinContent(atQ_130Te))
	N_at130Te.append(hist_2.GetBinContent(atQ_130Te))
	N_at130Te.append(hist_3.GetBinContent(atQ_130Te))
	N_at130Te.append(hist_4.GetBinContent(atQ_130Te))
	N_at130Te.append(hist_5.GetBinContent(atQ_130Te))
	N_at130Te.append(hist_6.GetBinContent(atQ_130Te))
	N_at130Te = unp.uarray(N_at130Te, stat_err_at130Te)

	stat_err_at116Cd = []
	stat_err_at130Te = []


	############################################
	### every known contributions of paint
	############################################
	print('\nCreating plot %s.pdf \n' %(savename[0]))
	# summing all histograms
	hist_comb  = ROOT.TH1F('hist_comb', '', bins, 0, E_range)
	hist_comb.Add(hist_comb, hist_1)
	hist_comb.Add(hist_comb, hist_2)
	hist_comb.Add(hist_comb, hist_3)
	hist_comb.Add(hist_comb, hist_4)
	hist_comb.Add(hist_comb, hist_5)
	hist_comb.Add(hist_comb, hist_6)
	hist_comb.SetLineColor(1)
	hist_comb.SetLineWidth(1)

	# adding all histograms
	hists.Add(hist_1)
	hists.Add(hist_2)
	hists.Add(hist_3)
	hists.Add(hist_4)
	hists.Add(hist_5)
	hists.Add(hist_6)
	hists.Add(hist_comb)

	hists.Draw("nostack")

	# nice labeling
	hists.GetXaxis().SetTitleSize(0.03)
	hists.GetXaxis().SetTitle("E / keV")
	hists.GetXaxis().SetTickLength(0.02)
	hists.GetXaxis().SetTicks("-")
	hists.GetXaxis().SetLabelOffset(-0.05)
	hists.GetXaxis().SetTitleOffset(-1.4)

	hists.GetYaxis().SetTitleSize(0.025)
	hists.GetYaxis().SetTitle('counts/kg keV yr')
	hists.GetYaxis().SetTickLength(0.02)
	hists.GetYaxis().SetTicks("+")
	hists.GetYaxis().SetLabelOffset(-0.025)
	hists.GetYaxis().SetTitleOffset(-1.8)

	hists.SetMinimum(1e-3)			# Set minimum of y-axis
	hists.SetMaximum(1e4)			# Set maximum of y-axis


	leg = ROOT.TLegend(0.375, 0.75, 0.9, 0.9)
	# Legend
	leg.SetHeader("Isotopes")
	leg.SetNColumns(4)
	leg.SetTextSize(0.025)

	leg.AddEntry('hist_1', "40K_glyptal", "f")
	leg.AddEntry('hist_2', "40K_epoxy", "f")
	leg.AddEntry('hist_3', "232Th_glyptal", "f")
	leg.AddEntry('hist_4', "232Th_epoxy", "f")
	leg.AddEntry('hist_5', "238U_glyptal", "f")
	leg.AddEntry('hist_6', "238U_epoxy", "f")
	leg.AddEntry('hist_comb', "combined", "f")
	leg.Draw()

	canvas.SetLogy()
	canvas.Update()
	canvas.Print(savename[0]+'.pdf')
	canvas.Clear()
	leg.Clear()

	############################################
	### only glyptal paint
	############################################
	hist_1.SetLineColor(2)
	hist_2.SetLineColor(2)
	hist_3.SetLineColor(3)
	hist_4.SetLineColor(3)
	hist_5.SetLineColor(4)
	hist_6.SetLineColor(4)

	print('\nCreating plot %s.pdf \n' %(savename[1]))
	hists_glyptal = ROOT.THStack("hists_glyptal", "")

	# summing histograms
	hist_comb_glyptal  = ROOT.TH1F('hist_comb_glyptal', '', bins, 0, E_range)
	hist_comb_glyptal.Add(hist_comb_glyptal, hist_1)
	hist_comb_glyptal.Add(hist_comb_glyptal, hist_3)
	hist_comb_glyptal.Add(hist_comb_glyptal, hist_5)
	hist_comb_glyptal.SetLineColor(1)
	hist_comb_glyptal.SetLineWidth(1)

	# adding all histograms
	hists_glyptal.Add(hist_1)
	hists_glyptal.Add(hist_3)
	hists_glyptal.Add(hist_5)
	hists_glyptal.Add(hist_comb_glyptal)
	hists_glyptal.Draw("nostack")

	# nice labeling
	hists_glyptal.GetXaxis().SetTitleSize(0.03)
	hists_glyptal.GetXaxis().SetTitle("E / keV")
	hists_glyptal.GetXaxis().SetTickLength(0.02)
	hists_glyptal.GetXaxis().SetTicks("-")
	hists_glyptal.GetXaxis().SetLabelOffset(-0.05)
	hists_glyptal.GetXaxis().SetTitleOffset(-1.4)

	hists_glyptal.GetYaxis().SetTitleSize(0.025)
	hists_glyptal.GetYaxis().SetTitle('counts/kg keV yr')
	hists_glyptal.GetYaxis().SetTickLength(0.02)
	hists_glyptal.GetYaxis().SetTicks("+")
	hists_glyptal.GetYaxis().SetLabelOffset(-0.025)
	hists_glyptal.GetYaxis().SetTitleOffset(-1.8)

	hists.SetMinimum(1e-3)			# Set minimum of y-axis
	hists.SetMaximum(1e4)			# Set maximum of y-axis

	# Legend
	leg = ROOT.TLegend(0.55, 0.8, 0.9, 0.9)
	leg.SetHeader("Isotopes from glyptal")
	leg.SetNColumns(4)
	leg.SetTextSize(0.025)

	leg.AddEntry('hist_1', "40K ", "f")
	leg.AddEntry('hist_3', "232Th ", "f")
	leg.AddEntry('hist_5', "238U ", "f")
	leg.AddEntry('hist_comb_glyptal', "combined", "f")
	leg.Draw()

	canvas.SetLogy()
	canvas.Update()
	canvas.Print(savename[1]+'.pdf')
	canvas.Clear()
	leg.Clear()


	############################################
	### only epoxy paint
	############################################
	print('\nCreating plot %s.pdf \n' %(savename[2]))
	hists_epoxy = ROOT.THStack("hists_epoxy", "")
	
	# summing histograms
	hist_comb_epoxy  = ROOT.TH1F('hist_comb_epoxy', '', bins, 0, E_range)
	hist_comb_epoxy.Add(hist_comb_epoxy, hist_2)
	hist_comb_epoxy.Add(hist_comb_epoxy, hist_4)
	hist_comb_epoxy.Add(hist_comb_epoxy, hist_6)
	hist_comb_epoxy.SetLineColor(1)
	hist_comb_epoxy.SetLineWidth(1)

	# adding all histograms
	hists_epoxy.Add(hist_2)
	hists_epoxy.Add(hist_4)
	hists_epoxy.Add(hist_6)
	hists_epoxy.Add(hist_comb_epoxy)
	hists_epoxy.Draw("nostack")

	# nice labeling
	hists_epoxy.GetXaxis().SetTitleSize(0.03)
	hists_epoxy.GetXaxis().SetTitle("E / keV")
	hists_epoxy.GetXaxis().SetTickLength(0.02)
	hists_epoxy.GetXaxis().SetTicks("-")
	hists_epoxy.GetXaxis().SetLabelOffset(-0.05)
	hists_epoxy.GetXaxis().SetTitleOffset(-1.4)

	hists_epoxy.GetYaxis().SetTitleSize(0.025)
	hists_epoxy.GetYaxis().SetTitle('counts/kg keV yr')
	hists_epoxy.GetYaxis().SetTickLength(0.02)
	hists_epoxy.GetYaxis().SetTicks("+")
	hists_epoxy.GetYaxis().SetLabelOffset(-0.025)
	hists_epoxy.GetYaxis().SetTitleOffset(-1.8)

	hists.SetMinimum(1e-3)			# Set minimum of y-axis
	hists.SetMaximum(1e4)			# Set maximum of y-axis

	# Legend
	leg = ROOT.TLegend(0.55, 0.8, 0.9, 0.9)
	leg.SetHeader("Isotopes from epoxy")
	leg.SetNColumns(4)
	leg.SetTextSize(0.025)

	leg.AddEntry('hist_2', "40K", "f")
	leg.AddEntry('hist_4', "232Th", "f")
	leg.AddEntry('hist_6', "238U", "f")
	leg.AddEntry('hist_comb_epoxy', "combined", "f")
	leg.Draw()

	canvas.SetLogy()
	canvas.Update()
	canvas.Print(savename[2]+'.pdf')
	canvas.Clear()
	leg.Clear()


	############################################
	### paint background, no paint distinction
	############################################
	print('\nCreating plot %s.pdf \n' %(savename[3]))
	hists_both = ROOT.THStack("hists_both", "")

	hist_40K = ROOT.TH1F('hist_40K', '', bins, 0, E_range)
	hist_40K.Add(hist_40K, hist_1)
	hist_40K.Add(hist_40K, hist_2)
	hist_40K.SetLineColor(2)
	hist_40K.SetLineWidth(2)

	hist_232Th  = ROOT.TH1F('hist_232Th', '', bins, 0, E_range)
	hist_232Th.Add(hist_232Th, hist_3)
	hist_232Th.Add(hist_232Th, hist_4)
	hist_232Th.SetLineColor(3)
	hist_232Th.SetLineWidth(2)

	hist_238U  = ROOT.TH1F('hist_238U', '', bins, 0, E_range)
	hist_238U.Add(hist_238U, hist_5)
	hist_238U.Add(hist_238U, hist_6)
	hist_238U.SetLineColor(4)
	hist_238U.SetLineWidth(2)

	# summing histograms
	hist_comb_both  = ROOT.TH1F('hist_comb_both', '', bins, 0, E_range)
	hist_comb_both.Add(hist_comb_both, hist_40K)
	hist_comb_both.Add(hist_comb_both, hist_232Th)
	hist_comb_both.Add(hist_comb_both, hist_238U)
	hist_comb_both.SetLineColor(1)
	hist_comb_both.SetLineWidth(1)

	# adding all histograms
	hists_both.Add(hist_40K)
	hists_both.Add(hist_232Th)
	hists_both.Add(hist_238U)
	hists_both.Add(hist_comb_both)
	hists_both.Draw("nostack")

	# nice labeling
	hists_both.GetXaxis().SetTitleSize(0.03)
	hists_both.GetXaxis().SetTitle("E / keV")
	hists_both.GetXaxis().SetTickLength(0.02)
	hists_both.GetXaxis().SetTicks("-")
	hists_both.GetXaxis().SetLabelOffset(-0.05)
	hists_both.GetXaxis().SetTitleOffset(-1.4)

	hists_both.GetYaxis().SetTitleSize(0.025)
	hists_both.GetYaxis().SetTitle('counts/kg keV yr')
	hists_both.GetYaxis().SetTickLength(0.02)
	hists_both.GetYaxis().SetTicks("+")
	hists_both.GetYaxis().SetLabelOffset(-0.025)
	hists_both.GetYaxis().SetTitleOffset(-1.8)

	hists.SetMinimum(1e-3)			# Set minimum of y-axis
	hists.SetMaximum(1e4)			# Set maximum of y-axis

	# Legend
	leg = ROOT.TLegend(0.55, 0.8, 0.9, 0.9)
	leg.SetHeader("Isotopes of both paints")
	leg.SetNColumns(4)
	leg.SetTextSize(0.025)

	leg.AddEntry('hist_40K', "40K", "f")
	leg.AddEntry('hist_232Th', "232Th", "f")
	leg.AddEntry('hist_238U', "238U", "f")
	leg.AddEntry('hist_comb_both', "combined", "f")
	leg.Draw()

	canvas.SetLogy()
	canvas.Update()
	canvas.Print(savename[3]+'.pdf')
	canvas.Clear()
	leg.Clear()


	############################################
	### difference plot
	############################################
	# this is the difference from the paint backgrounds!!!
	print('\nCreating plot %s.pdf \n' %(savename[4]))
	hists_diff = ROOT.THStack("hists_diff", "Differences between detector paints: N_epoxy - N_glyptal")

	hist_40K_diff  = ROOT.TH1F('hist_40K_diff', '', bins, 0, E_range)
	hist_40K_diff.Add(hist_1)
	hist_40K_diff.Add(hist_2, -1)
	hist_40K_diff.SetLineColor(2)
	hist_40K_diff.SetLineWidth(2)

	hist_232Th_diff  = ROOT.TH1F('hist_232Th_diff', '', bins, 0, E_range)
	hist_232Th_diff.Add(hist_3)
	hist_232Th_diff.Add(hist_4, -1)
	hist_232Th_diff.SetLineColor(3)
	hist_232Th_diff.SetLineWidth(2)

	hist_238U_diff  = ROOT.TH1F('hist_238U_diff', '', bins, 0, E_range)
	hist_238U_diff.Add( hist_5)
	hist_238U_diff.Add(hist_6, -1)
	hist_238U_diff.SetLineColor(4)
	hist_238U_diff.SetLineWidth(2)

	# difference between combind glyptal events and combined glptal events
	hist_comb_diff  = ROOT.TH1F('hist_comb_diff', '', bins, 0, E_range)
	hist_comb_diff.Add(hist_40K_diff)
	hist_comb_diff.Add(hist_232Th_diff)
	hist_comb_diff.Add(hist_238U_diff)
	hist_comb_diff.SetLineColor(1)
	hist_comb_diff.SetLineWidth(1)
	#hist_comb_diff.SetLineStyle(3)

	# adding all histograms
	hists_diff.Add(hist_40K_diff)
	hists_diff.Add(hist_232Th_diff)
	hists_diff.Add(hist_238U_diff)
	hists_diff.Add(hist_comb_diff)
	hists_diff.Draw("nostack")

	# nice labeling
	hists_diff.GetXaxis().SetTitleSize(0.03)
	hists_diff.GetXaxis().SetTitle("E / keV")
	hists_diff.GetXaxis().SetTickLength(0.02)
	hists_diff.GetXaxis().SetTicks("-")
	hists_diff.GetXaxis().SetLabelOffset(-0.05)
	hists_diff.GetXaxis().SetTitleOffset(-1.4)

	hists_diff.GetYaxis().SetTitleSize(0.025)
	hists_diff.GetYaxis().SetTitle('counts/kg keV yr')
	hists_diff.GetYaxis().SetTickLength(0.02)
	hists_diff.GetYaxis().SetTicks("+")
	hists_diff.GetYaxis().SetLabelOffset(-0.025)
	hists_diff.GetYaxis().SetTitleOffset(-1.8)

	hists_diff.SetMinimum(-20)			# Set minimum of y-axis
	hists_diff.SetMaximum(20)			# Set maximum of y-axis

	# Legend
	leg = ROOT.TLegend(0.55, 0.8, 0.9, 0.9)
	leg.SetHeader("Isotopes")
	leg.SetNColumns(4)
	leg.SetTextSize(0.025)

	leg.AddEntry('hist_40K_diff', "40K", "f")
	leg.AddEntry('hist_232Th_diff', "232Th", "f")
	leg.AddEntry('hist_238U_diff', "238U", "f")
	leg.AddEntry('hist_comb_diff', "combined", "f")
	leg.Draw()

	canvas.SetLogy(False)
	canvas.Update()
	canvas.Print(savename[4]+'.pdf')

	hists_diff.SetMinimum(-1)			# Set minimum of y-axis
	hists_diff.SetMaximum(1)			# Set maximum of y-axis

	canvas.Update()
	canvas.Print(savename[4]+'_2.pdf')
	canvas.Clear()
	leg.Clear()

	#################################################
	# differentiate between alpha, beta, gamma and other particles
	#################################################
	# particle_id meanings: "alpha"=E_range, "electron"=11 , "gamma"=22
	alpha_events = []
	beta_events = []
	gamma_events = []
	other_events = []
	temp_alpha = []
	temp_beta = []
	temp_gamma = []
	temp_other = []
	#print('len(particle_id): ', len(particle_id))
	#for i in range(len(particle_id)):
	#	print('len(particle_id[%i]): %i' %(i, len(particle_id[i])))

	for i in range(len(particle_id)):
		for j in range(len(particle_id[i])):
			temp = all_events[i][j]
			if (particle_id[i][j] == E_range):
				temp_alpha.append(temp)
			elif (particle_id[i][j] == 11):
				temp_beta.append(temp)
			elif (particle_id[i][j] == 22):
				temp_gamma.append(temp)
			else:
				temp_other.append(temp)
		alpha_events.append(temp_alpha)
		beta_events.append(temp_beta)
		gamma_events.append(temp_gamma)
		other_events.append(temp_other)
		temp_alpha = []
		temp_beta = []
		temp_gamma = []
		temp_other = []

	#print('alpha:')
	#for i in range(len(alpha_events)):
	#	print(len(alpha_events[i]))
	#print()
	#print('beta:')
	#for i in range(len(beta_events)):
	#	print(len(beta_events[i]))
	#print()
	#print('gamma:')
	#for i in range(len(gamma_events)):
	#	print(len(gamma_events[i]))
	#print()
	#print('other:')
	#for i in range(len(other_events)):
	#	print(len(other_events[i]))
	#print()

	hists_particle_sorts = ROOT.THStack("hists_particle_sorts", "")
	hists_particle_sorts_epox = ROOT.THStack("hists_particle_sorts_epox", "")
	hists_particle_sorts_glyp = ROOT.THStack("hists_particle_sorts_gylp", "")
	# need to create single histiograms and scaling each, then add to single particle sort

# alpha
	hist_alpha_40K_glyp = ROOT.TH1F('hist_alpha_40K_glyp', '', bins, 0, E_range)	
	hist_alpha_40K_epox = ROOT.TH1F('hist_alpha_40K_epox', '', bins, 0, E_range)	
	hist_alpha_232Th_glyp = ROOT.TH1F('hist_alpha_232Th_glyp', '', bins, 0, E_range)	
	hist_alpha_232Th_epox = ROOT.TH1F('hist_alpha_232Th_epox', '', bins, 0, E_range)	
	hist_alpha_238U_glyp = ROOT.TH1F('hist_alpha_238U_glyp', '', bins, 0, E_range)	
	hist_alpha_238U_epox = ROOT.TH1F('hist_alpha_238U_epox', '', bins, 0, E_range)
	for j in range(len(alpha_events[0])):
		hist_alpha_40K_glyp.Fill(alpha_events[0][j])
	stat_err_at116Cd_temp = np.sqrt(hist_alpha_40K_glyp.GetBinContent(atQ_116Cd))
	stat_err_at130Te_temp = np.sqrt(hist_alpha_40K_glyp.GetBinContent(atQ_130Te))
	#hist_4.SetFillColor(6)
	#hist_4.SetFillStyle(1001)
	hist_alpha_40K_glyp.SetLineColor(2)
	hist_alpha_40K_glyp.SetLineWidth(2)
	# Scaling
	scale = n_chain[0] * N_norm[0]/(1e6)		# for 1 million events # divide per bin width for per keV
	stat_err_at116Cd_temp = scale * stat_err_at116Cd_temp
	stat_err_at116Cd.append(stat_err_at116Cd_temp)
	stat_err_at130Te_temp = scale * stat_err_at130Te_temp
	stat_err_at130Te.append(stat_err_at130Te_temp)
	
	hist_alpha_40K_glyp.Scale(scale)

	for j in range(len(alpha_events[1])):
		hist_alpha_40K_epox.Fill(alpha_events[1][j])
	stat_err_at116Cd_temp = np.sqrt(hist_alpha_40K_epox.GetBinContent(atQ_116Cd))
	stat_err_at130Te_temp = np.sqrt(hist_alpha_40K_epox.GetBinContent(atQ_130Te))
	#hist_4.SetFillColor(6)
	#hist_4.SetFillStyle(1001)
	hist_alpha_40K_epox.SetLineColor(2)
	hist_alpha_40K_epox.SetLineWidth(2)
	# Scaling
	scale = n_chain[0] * N_norm[1]/(1e6)		# for 1 million events # divide per bin width for per keV
	stat_err_at116Cd_temp = scale * stat_err_at116Cd_temp
	stat_err_at116Cd.append(stat_err_at116Cd_temp)
	stat_err_at130Te_temp = scale * stat_err_at130Te_temp
	stat_err_at130Te.append(stat_err_at130Te_temp)
	
	hist_alpha_40K_epox.Scale(scale)

	for j in range(len(alpha_events[2])):
		hist_alpha_232Th_glyp.Fill(alpha_events[2][j])
	stat_err_at116Cd_temp = np.sqrt(hist_alpha_232Th_glyp.GetBinContent(atQ_116Cd))
	stat_err_at130Te_temp = np.sqrt(hist_alpha_232Th_glyp.GetBinContent(atQ_130Te))
	hist_alpha_232Th_glyp.SetLineColor(2)
	hist_alpha_232Th_glyp.SetLineWidth(2)
	# Scaling
	scale = n_chain[1] * N_norm[2]/(1e6)		# for 1 million events # divide per bin width for per keV
	stat_err_at116Cd_temp = scale * stat_err_at116Cd_temp
	stat_err_at116Cd.append(stat_err_at116Cd_temp)
	stat_err_at130Te_temp = scale * stat_err_at130Te_temp
	stat_err_at130Te.append(stat_err_at130Te_temp)
	
	hist_alpha_232Th_glyp.Scale(scale)

	for j in range(len(alpha_events[3])):
		hist_alpha_232Th_epox.Fill(alpha_events[3][j])
	stat_err_at116Cd_temp = np.sqrt(hist_alpha_232Th_epox.GetBinContent(atQ_116Cd))
	stat_err_at130Te_temp = np.sqrt(hist_alpha_232Th_epox.GetBinContent(atQ_130Te))
	hist_alpha_232Th_epox.SetLineColor(2)
	hist_alpha_232Th_epox.SetLineWidth(2)
	# Scaling
	scale = n_chain[1] * N_norm[3]/(1e6)		# for 1 million events # divide per bin width for per keV
	stat_err_at116Cd_temp = scale * stat_err_at116Cd_temp
	stat_err_at116Cd.append(stat_err_at116Cd_temp)
	stat_err_at130Te_temp = scale * stat_err_at130Te_temp
	stat_err_at130Te.append(stat_err_at130Te_temp)
	
	hist_alpha_232Th_epox.Scale(scale)

	for j in range(len(alpha_events[4])):
		hist_alpha_238U_glyp.Fill(alpha_events[4][j])
	stat_err_at116Cd_temp = np.sqrt(hist_alpha_238U_glyp.GetBinContent(atQ_116Cd))
	stat_err_at130Te_temp = np.sqrt(hist_alpha_238U_glyp.GetBinContent(atQ_130Te))
	hist_alpha_238U_glyp.SetLineColor(2)
	hist_alpha_238U_glyp.SetLineWidth(2)
	# Scaling
	scale = n_chain[2] * N_norm[4]/(1e6)		# for 1 million events # divide per bin width for per keV
	stat_err_at116Cd_temp = scale * stat_err_at116Cd_temp
	stat_err_at116Cd.append(stat_err_at116Cd_temp)
	stat_err_at130Te_temp = scale * stat_err_at130Te_temp
	stat_err_at130Te.append(stat_err_at130Te_temp)
	
	hist_alpha_238U_glyp.Scale(scale)

	for j in range(len(alpha_events[5])):
		hist_alpha_238U_epox.Fill(alpha_events[5][j])
	stat_err_at116Cd_temp = np.sqrt(hist_alpha_238U_epox.GetBinContent(atQ_116Cd))
	stat_err_at130Te_temp = np.sqrt(hist_alpha_238U_epox.GetBinContent(atQ_130Te))
	hist_alpha_238U_epox.SetLineColor(2)
	hist_alpha_238U_epox.SetLineWidth(2)
	# Scaling
	scale = n_chain[2] * N_norm[5]/(1e6)		# for 1 million events # divide per bin width for per keV
	stat_err_at116Cd_temp = scale * stat_err_at116Cd_temp
	stat_err_at116Cd.append(stat_err_at116Cd_temp)
	stat_err_at130Te_temp = scale * stat_err_at130Te_temp
	stat_err_at130Te.append(stat_err_at130Te_temp)
	
	hist_alpha_238U_epox.Scale(scale)


	# Errors for each coating
	stat_err_at116Cd_alpha_glyp = np.mean([stat_err_at116Cd[0], stat_err_at116Cd[2], stat_err_at116Cd[4]])
	stat_err_at116Cd_alpha_epox = np.mean([stat_err_at116Cd[1], stat_err_at116Cd[3], stat_err_at116Cd[5]])
	stat_err_at130Te_alpha_glyp = np.mean([stat_err_at130Te[0], stat_err_at130Te[2], stat_err_at130Te[4]])
	stat_err_at130Te_alpha_epox = np.mean([stat_err_at130Te[1], stat_err_at130Te[3], stat_err_at130Te[5]])

	stat_err_at116Cd = []
	stat_err_at130Te = []



	# Adding to one histogram
	# all
	hist_alpha = ROOT.TH1F('hist_alpha', '', bins, 0, E_range)
	hist_alpha.Add(hist_alpha_40K_glyp)
	hist_alpha.Add(hist_alpha_40K_epox)
	hist_alpha.Add(hist_alpha_232Th_glyp)
	hist_alpha.Add(hist_alpha_232Th_epox)
	hist_alpha.Add(hist_alpha_238U_glyp)
	hist_alpha.Add(hist_alpha_238U_epox)
	hist_alpha.SetLineColor(7)
	hist_alpha.SetLineWidth(2)

	# epoxy
	hist_alpha_epoxy = ROOT.TH1F('hist_alpha_epoxy', '', bins, 0, E_range)
	hist_alpha_epoxy.Add(hist_alpha_40K_epox)
	hist_alpha_epoxy.Add(hist_alpha_232Th_epox)
	hist_alpha_epoxy.Add(hist_alpha_238U_epox)
	hist_alpha_epoxy.SetLineColor(7)
	hist_alpha_epoxy.SetLineWidth(2)

	# glyptal
	hist_alpha_glyptal = ROOT.TH1F('hist_alpha_glyptal', '', bins, 0, E_range)
	hist_alpha_glyptal.Add(hist_alpha_40K_glyp)
	hist_alpha_glyptal.Add(hist_alpha_232Th_glyp)
	hist_alpha_glyptal.Add(hist_alpha_238U_glyp)
	hist_alpha_glyptal.SetLineColor(7)
	hist_alpha_glyptal.SetLineWidth(2)


# beta
	hist_beta_40K_glyp = ROOT.TH1F('hist_beta_40K_glyp', '', bins, 0, E_range)	
	hist_beta_40K_epox = ROOT.TH1F('hist_beta_40K_epox', '', bins, 0, E_range)	
	hist_beta_232Th_glyp = ROOT.TH1F('hist_beta_232Th_glyp', '', bins, 0, E_range)	
	hist_beta_232Th_epox = ROOT.TH1F('hist_beta_232Th_epox', '', bins, 0, E_range)	
	hist_beta_238U_glyp = ROOT.TH1F('hist_beta_238U_glyp', '', bins, 0, E_range)	
	hist_beta_238U_epox = ROOT.TH1F('hist_beta_238U_epox', '', bins, 0, E_range)
	for j in range(len(beta_events[0])):
		hist_beta_40K_glyp.Fill(beta_events[0][j])
	stat_err_at116Cd_temp = np.sqrt(hist_beta_40K_glyp.GetBinContent(atQ_116Cd))
	stat_err_at130Te_temp = np.sqrt(hist_beta_40K_glyp.GetBinContent(atQ_130Te))
	#hist_4.SetFillColor(6)
	#hist_4.SetFillStyle(1001)
	hist_beta_40K_glyp.SetLineColor(3)
	hist_beta_40K_glyp.SetLineWidth(2)
	# Scaling
	scale = n_chain[0] * N_norm[0]/(1e6)		# for 1 million events # divide per bin width for per keV
	stat_err_at116Cd_temp = scale * stat_err_at116Cd_temp
	stat_err_at116Cd.append(stat_err_at116Cd_temp)
	stat_err_at130Te_temp = scale * stat_err_at130Te_temp
	stat_err_at130Te.append(stat_err_at130Te_temp)
	
	hist_beta_40K_glyp.Scale(scale)

	for j in range(len(beta_events[1])):
		hist_beta_40K_epox.Fill(beta_events[1][j])
	stat_err_at116Cd_temp = np.sqrt(hist_beta_40K_epox.GetBinContent(atQ_116Cd))
	stat_err_at130Te_temp = np.sqrt(hist_beta_40K_epox.GetBinContent(atQ_130Te))
	#hist_4.SetFillColor(6)
	#hist_4.SetFillStyle(1001)
	hist_beta_40K_epox.SetLineColor(3)
	hist_beta_40K_epox.SetLineWidth(2)
	# Scaling
	scale = n_chain[0] * N_norm[1]/(1e6)		# for 1 million events # divide per bin width for per keV
	stat_err_at116Cd_temp = scale * stat_err_at116Cd_temp
	stat_err_at116Cd.append(stat_err_at116Cd_temp)
	stat_err_at130Te_temp = scale * stat_err_at130Te_temp
	stat_err_at130Te.append(stat_err_at130Te_temp)
	
	hist_beta_40K_epox.Scale(scale)

	for j in range(len(beta_events[2])):
		hist_beta_232Th_glyp.Fill(beta_events[2][j])
	stat_err_at116Cd_temp = np.sqrt(hist_beta_232Th_glyp.GetBinContent(atQ_116Cd))
	stat_err_at130Te_temp = np.sqrt(hist_beta_232Th_glyp.GetBinContent(atQ_130Te))
	hist_beta_232Th_glyp.SetLineColor(3)
	hist_beta_232Th_glyp.SetLineWidth(2)
	# Scaling
	scale = n_chain[1] * N_norm[2]/(1e6)		# for 1 million events # divide per bin width for per keV
	stat_err_at116Cd_temp = scale * stat_err_at116Cd_temp
	stat_err_at116Cd.append(stat_err_at116Cd_temp)
	stat_err_at130Te_temp = scale * stat_err_at130Te_temp
	stat_err_at130Te.append(stat_err_at130Te_temp)
	
	hist_beta_232Th_glyp.Scale(scale)

	for j in range(len(beta_events[3])):
		hist_beta_232Th_epox.Fill(beta_events[3][j])
	stat_err_at116Cd_temp = np.sqrt(hist_beta_232Th_epox.GetBinContent(atQ_116Cd))
	stat_err_at130Te_temp = np.sqrt(hist_beta_232Th_epox.GetBinContent(atQ_130Te))
	hist_beta_232Th_epox.SetLineColor(3)
	hist_beta_232Th_epox.SetLineWidth(2)
	# Scaling
	scale = n_chain[1] * N_norm[3]/(1e6)		# for 1 million events # divide per bin width for per keV
	stat_err_at116Cd_temp = scale * stat_err_at116Cd_temp
	stat_err_at116Cd.append(stat_err_at116Cd_temp)
	stat_err_at130Te_temp = scale * stat_err_at130Te_temp
	stat_err_at130Te.append(stat_err_at130Te_temp)
	
	hist_beta_232Th_epox.Scale(scale)

	for j in range(len(beta_events[4])):
		hist_beta_238U_glyp.Fill(beta_events[4][j])
	stat_err_at116Cd_temp = np.sqrt(hist_beta_238U_glyp.GetBinContent(atQ_116Cd))
	stat_err_at130Te_temp = np.sqrt(hist_beta_238U_glyp.GetBinContent(atQ_130Te))
	hist_beta_238U_glyp.SetLineColor(3)
	hist_beta_238U_glyp.SetLineWidth(2)
	# Scaling
	scale = n_chain[2] * N_norm[4]/(1e6)		# for 1 million events # divide per bin width for per keV
	stat_err_at116Cd_temp = scale * stat_err_at116Cd_temp
	stat_err_at116Cd.append(stat_err_at116Cd_temp)
	stat_err_at130Te_temp = scale * stat_err_at130Te_temp
	stat_err_at130Te.append(stat_err_at130Te_temp)
	
	hist_beta_238U_glyp.Scale(scale)

	for j in range(len(beta_events[5])):
		hist_beta_238U_epox.Fill(beta_events[5][j])
	stat_err_at116Cd_temp = np.sqrt(hist_beta_238U_glyp.GetBinContent(atQ_116Cd))
	stat_err_at130Te_temp = np.sqrt(hist_beta_238U_glyp.GetBinContent(atQ_130Te))
	hist_beta_238U_epox.SetLineColor(3)
	hist_beta_238U_epox.SetLineWidth(2)
	# Scaling
	scale = n_chain[2] * N_norm[5]/(1e6)		# for 1 million events # divide per bin width for per keV
	stat_err_at116Cd_temp = scale * stat_err_at116Cd_temp
	stat_err_at116Cd.append(stat_err_at116Cd_temp)
	stat_err_at130Te_temp = scale * stat_err_at130Te_temp
	stat_err_at130Te.append(stat_err_at130Te_temp)
	
	hist_beta_238U_epox.Scale(scale)

	# Errors for each coating
	stat_err_at116Cd_beta_glyp = np.mean([stat_err_at116Cd[0], stat_err_at116Cd[2], stat_err_at116Cd[4]])
	stat_err_at116Cd_beta_epox = np.mean([stat_err_at116Cd[1], stat_err_at116Cd[3], stat_err_at116Cd[5]])
	stat_err_at130Te_beta_glyp = np.mean([stat_err_at130Te[0], stat_err_at130Te[2], stat_err_at130Te[4]])
	stat_err_at130Te_beta_epox = np.mean([stat_err_at130Te[1], stat_err_at130Te[3], stat_err_at130Te[5]])

	stat_err_at116Cd = []
	stat_err_at130Te = []


	# Adding to one histogram
	hist_beta = ROOT.TH1F('hist_beta', '', bins, 0, E_range)
	hist_beta.Add(hist_beta_40K_glyp)
	hist_beta.Add(hist_beta_40K_epox)
	hist_beta.Add(hist_beta_232Th_glyp)
	hist_beta.Add(hist_beta_232Th_epox)
	hist_beta.Add(hist_beta_238U_glyp)
	hist_beta.Add(hist_beta_238U_epox)
	hist_beta.SetLineColor(6)
	hist_beta.SetLineWidth(2)

	# epoxy
	hist_beta_epoxy = ROOT.TH1F('hist_beta_epoxy', '', bins, 0, E_range)
	hist_beta_epoxy.Add(hist_beta_40K_epox)
	hist_beta_epoxy.Add(hist_beta_232Th_epox)
	hist_beta_epoxy.Add(hist_beta_238U_epox)
	hist_beta_epoxy.SetLineColor(6)
	hist_beta_epoxy.SetLineWidth(2)

	# glyptal
	hist_beta_glyptal = ROOT.TH1F('hist_beta_glyptal', '', bins, 0, E_range)
	hist_beta_glyptal.Add(hist_beta_40K_glyp)
	hist_beta_glyptal.Add(hist_beta_232Th_glyp)
	hist_beta_glyptal.Add(hist_beta_238U_glyp)
	hist_beta_glyptal.SetLineColor(6)
	hist_beta_glyptal.SetLineWidth(2)


# gamma
	hist_gamma_40K_glyp = ROOT.TH1F('hist_gamma_40K_glyp', '', bins, 0, E_range)	
	hist_gamma_40K_epox = ROOT.TH1F('hist_gamma_40K_epox', '', bins, 0, E_range)	
	hist_gamma_232Th_glyp = ROOT.TH1F('hist_gamma_232Th_glyp', '', bins, 0, E_range)	
	hist_gamma_232Th_epox = ROOT.TH1F('hist_gamma_232Th_epox', '', bins, 0, E_range)	
	hist_gamma_238U_glyp = ROOT.TH1F('hist_gamma_238U_glyp', '', bins, 0, E_range)	
	hist_gamma_238U_epox = ROOT.TH1F('hist_gamma_238U_epox', '', bins, 0, E_range)
	for j in range(len(gamma_events[0])):
		hist_gamma_40K_glyp.Fill(gamma_events[0][j])
	stat_err_at116Cd_temp = np.sqrt(hist_gamma_40K_glyp.GetBinContent(atQ_116Cd))
	stat_err_at130Te_temp = np.sqrt(hist_gamma_40K_glyp.GetBinContent(atQ_130Te))
	#hist_4.SetFillColor(6)
	#hist_4.SetFillStyle(1001)
	hist_gamma_40K_glyp.SetLineColor(4)
	hist_gamma_40K_glyp.SetLineWidth(2)
	# Scaling
	scale = n_chain[0] * N_norm[0]/(1e6)		# for 1 million events # divide per bin width for per keV
	stat_err_at116Cd_temp = scale * stat_err_at116Cd_temp
	stat_err_at116Cd.append(stat_err_at116Cd_temp)
	stat_err_at130Te_temp = scale * stat_err_at130Te_temp
	stat_err_at130Te.append(stat_err_at130Te_temp)
	
	hist_gamma_40K_glyp.Scale(scale)

	for j in range(len(gamma_events[1])):
		hist_gamma_40K_epox.Fill(gamma_events[1][j])
	stat_err_at116Cd_temp = np.sqrt(hist_gamma_40K_epox.GetBinContent(atQ_116Cd))
	stat_err_at130Te_temp = np.sqrt(hist_gamma_40K_epox.GetBinContent(atQ_130Te))
	#hist_4.SetFillColor(6)
	#hist_4.SetFillStyle(1001)
	hist_gamma_40K_epox.SetLineColor(4)
	hist_gamma_40K_epox.SetLineWidth(2)
	# Scaling
	scale = n_chain[0] * N_norm[1]/(1e6)		# for 1 million events # divide per bin width for per keV
	stat_err_at116Cd_temp = scale * stat_err_at116Cd_temp
	stat_err_at116Cd.append(stat_err_at116Cd_temp)
	stat_err_at130Te_temp = scale * stat_err_at130Te_temp
	stat_err_at130Te.append(stat_err_at130Te_temp)
	
	hist_gamma_40K_epox.Scale(scale)

	for j in range(len(gamma_events[2])):
		hist_gamma_232Th_glyp.Fill(gamma_events[2][j])
	stat_err_at116Cd_temp = np.sqrt(hist_gamma_232Th_glyp.GetBinContent(atQ_116Cd))
	stat_err_at130Te_temp = np.sqrt(hist_gamma_232Th_glyp.GetBinContent(atQ_130Te))
	hist_gamma_232Th_glyp.SetLineColor(4)
	hist_gamma_232Th_glyp.SetLineWidth(2)
	# Scaling
	scale = n_chain[1] * N_norm[2]/(1e6)		# for 1 million events # divide per bin width for per keV
	stat_err_at116Cd_temp = scale * stat_err_at116Cd_temp
	stat_err_at116Cd.append(stat_err_at116Cd_temp)
	stat_err_at130Te_temp = scale * stat_err_at130Te_temp
	stat_err_at130Te.append(stat_err_at130Te_temp)
	
	hist_gamma_232Th_glyp.Scale(scale)

	for j in range(len(gamma_events[3])):
		hist_gamma_232Th_epox.Fill(gamma_events[3][j])
	stat_err_at116Cd_temp = np.sqrt(hist_gamma_232Th_epox.GetBinContent(atQ_116Cd))
	stat_err_at130Te_temp = np.sqrt(hist_gamma_232Th_epox.GetBinContent(atQ_130Te))
	hist_gamma_232Th_epox.SetLineColor(4)
	hist_gamma_232Th_epox.SetLineWidth(2)
	# Scaling
	scale = n_chain[1] * N_norm[3]/(1e6)		# for 1 million events # divide per bin width for per keV
	stat_err_at116Cd_temp = scale * stat_err_at116Cd_temp
	stat_err_at116Cd.append(stat_err_at116Cd_temp)
	stat_err_at130Te_temp = scale * stat_err_at130Te_temp
	stat_err_at130Te.append(stat_err_at130Te_temp)
	
	hist_gamma_232Th_epox.Scale(scale)

	for j in range(len(gamma_events[4])):
		hist_gamma_238U_glyp.Fill(gamma_events[4][j])
	stat_err_at116Cd_temp = np.sqrt(hist_gamma_238U_glyp.GetBinContent(atQ_116Cd))
	stat_err_at130Te_temp = np.sqrt(hist_gamma_238U_glyp.GetBinContent(atQ_130Te))
	hist_gamma_238U_glyp.SetLineColor(4)
	hist_gamma_238U_glyp.SetLineWidth(2)
	# Scaling
	scale = n_chain[2] * N_norm[4]/(1e6)		# for 1 million events # divide per bin width for per keV
	stat_err_at116Cd_temp = scale * stat_err_at116Cd_temp
	stat_err_at116Cd.append(stat_err_at116Cd_temp)
	stat_err_at130Te_temp = scale * stat_err_at130Te_temp
	stat_err_at130Te.append(stat_err_at130Te_temp)
	
	hist_gamma_238U_glyp.Scale(scale)

	for j in range(len(gamma_events[5])):
		hist_gamma_238U_epox.Fill(gamma_events[5][j])
	stat_err_at116Cd_temp = np.sqrt(hist_gamma_238U_epox.GetBinContent(atQ_116Cd))
	stat_err_at130Te_temp = np.sqrt(hist_gamma_238U_epox.GetBinContent(atQ_130Te))
	hist_gamma_238U_epox.SetLineColor(4)
	hist_gamma_238U_epox.SetLineWidth(2)
	# Scaling
	scale = n_chain[2] * N_norm[5]/(1e6)		# for 1 million events # divide per bin width for per keV
	stat_err_at116Cd_temp = scale * stat_err_at116Cd_temp
	stat_err_at116Cd.append(stat_err_at116Cd_temp)
	stat_err_at130Te_temp = scale * stat_err_at130Te_temp
	stat_err_at130Te.append(stat_err_at130Te_temp)

	hist_gamma_238U_epox.Scale(scale)

	# Saving errors for each coating
	stat_err_at116Cd_gamma_glyp = np.mean([stat_err_at116Cd[0], stat_err_at116Cd[2], stat_err_at116Cd[4]])
	stat_err_at116Cd_gamma_epox = np.mean([stat_err_at116Cd[1], stat_err_at116Cd[3], stat_err_at116Cd[5]])
	stat_err_at130Te_gamma_glyp = np.mean([stat_err_at130Te[0], stat_err_at130Te[2], stat_err_at130Te[4]])
	stat_err_at130Te_gamma_epox = np.mean([stat_err_at130Te[1], stat_err_at130Te[3], stat_err_at130Te[5]])


	print()
	print('Checkpoint')
	print(stat_err_at116Cd_gamma_glyp)
	print(stat_err_at116Cd_gamma_epox)
	print(stat_err_at130Te_gamma_glyp)
	print(stat_err_at130Te_gamma_epox)

	stat_err_at116Cd = []
	stat_err_at130Te = []

	# Adding to one histogram
	hist_gamma = ROOT.TH1F('hist_gamma', '', bins, 0, E_range)
	hist_gamma.Add(hist_gamma_40K_glyp)
	hist_gamma.Add(hist_gamma_40K_epox)
	hist_gamma.Add(hist_gamma_232Th_glyp)
	hist_gamma.Add(hist_gamma_232Th_epox)
	hist_gamma.Add(hist_gamma_238U_glyp)
	hist_gamma.Add(hist_gamma_238U_epox)
	hist_gamma.SetLineColor(5)
	hist_gamma.SetLineWidth(2)

	# epoxy
	hist_gamma_epoxy = ROOT.TH1F('hist_gamma_epoxy', '', bins, 0, E_range)
	hist_gamma_epoxy.Add(hist_gamma_40K_epox)
	hist_gamma_epoxy.Add(hist_gamma_232Th_epox)
	hist_gamma_epoxy.Add(hist_gamma_238U_epox)
	hist_gamma_epoxy.SetLineColor(5)
	hist_gamma_epoxy.SetLineWidth(2)

	# glyptal
	hist_gamma_glyptal = ROOT.TH1F('hist_gamma_glyptal', '', bins, 0, E_range)
	hist_gamma_glyptal.Add(hist_gamma_40K_glyp)
	hist_gamma_glyptal.Add(hist_gamma_232Th_glyp)
	hist_gamma_glyptal.Add(hist_gamma_238U_glyp)
	hist_gamma_glyptal.SetLineColor(5)
	hist_gamma_glyptal.SetLineWidth(2)

	# Saving statistical errors for each particle
	stat_err_at116Cd_particles = [stat_err_at116Cd_alpha_glyp, stat_err_at116Cd_alpha_epox, stat_err_at116Cd_beta_glyp, stat_err_at116Cd_beta_epox, stat_err_at116Cd_gamma_glyp, stat_err_at116Cd_gamma_epox]
	stat_err_at130Te_particles = [stat_err_at130Te_alpha_glyp, stat_err_at130Te_alpha_epox, stat_err_at130Te_beta_glyp, stat_err_at130Te_beta_epox, stat_err_at130Te_gamma_glyp, stat_err_at130Te_gamma_epox]

	#print(stat_err_at116Cd_particles)
	#print(stat_err_at130Te_particles)


# create combined histogram
	hist_comb_particles  = ROOT.TH1F('hist_comb_particles', '', bins, 0, E_range)
	hist_comb_particles.Add(hist_alpha)
	hist_comb_particles.Add(hist_beta)
	hist_comb_particles.Add(hist_gamma)
	hist_comb_particles.SetLineColor(1)
	hist_comb_particles.SetLineWidth(1)

	# adding all histograms
	hists_particle_sorts.Add(hist_alpha)
	hists_particle_sorts.Add(hist_beta)
	hists_particle_sorts.Add(hist_gamma)
	hists_particle_sorts.Add(hist_comb_particles)
	hists_particle_sorts.Draw("nostack")

	# nice labeling
	hists_particle_sorts.GetXaxis().SetTitleSize(0.03)
	hists_particle_sorts.GetXaxis().SetTitle("E / keV")
	hists_particle_sorts.GetXaxis().SetTickLength(0.02)
	hists_particle_sorts.GetXaxis().SetTicks("-")
	hists_particle_sorts.GetXaxis().SetLabelOffset(-0.05)
	hists_particle_sorts.GetXaxis().SetTitleOffset(-1.4)

	hists_particle_sorts.GetYaxis().SetTitleSize(0.025)
	hists_particle_sorts.GetYaxis().SetTitle('counts/kg keV yr')
	hists_particle_sorts.GetYaxis().SetTickLength(0.02)
	hists_particle_sorts.GetYaxis().SetTicks("+")
	hists_particle_sorts.GetYaxis().SetLabelOffset(-0.025)
	hists_particle_sorts.GetYaxis().SetTitleOffset(-1.8)

	hists.SetMinimum(1e-3)			# Set minimum of y-axis
	hists.SetMaximum(1e4)			# Set maximum of y-axis

	# Legend
	leg = ROOT.TLegend(0.55, 0.8, 0.9, 0.9)
	leg.SetHeader("Particles from both paints")
	leg.SetNColumns(4)
	leg.SetTextSize(0.025)

	leg.AddEntry('hist_alpha', "alphas", "f")
	leg.AddEntry('hist_beta', "betas", "f")
	leg.AddEntry('hist_gamma', "gammas", "f")
	leg.AddEntry('hist_comb_particles', "combined", "f")
	leg.Draw()

	canvas.SetLogy()
	canvas.Update()
	canvas.Print(savename[5]+'.pdf')
	canvas.Clear()
	leg.Clear()

	# Getting entries at Q-Values
	N_at116Cd_particles.append(hist_alpha_glyptal.GetBinContent(atQ_116Cd))
	N_at116Cd_particles.append(hist_alpha_epoxy.GetBinContent(atQ_116Cd))
	N_at116Cd_particles.append(hist_beta_glyptal.GetBinContent(atQ_116Cd))
	N_at116Cd_particles.append(hist_beta_epoxy.GetBinContent(atQ_116Cd))
	N_at116Cd_particles.append(hist_gamma_glyptal.GetBinContent(atQ_116Cd))
	N_at116Cd_particles.append(hist_gamma_epoxy.GetBinContent(atQ_116Cd))
	N_at116Cd_particles = unp.uarray(N_at116Cd_particles, stat_err_at116Cd_particles)

	N_at130Te_particles.append(hist_alpha_glyptal.GetBinContent(atQ_130Te))
	N_at130Te_particles.append(hist_alpha_epoxy.GetBinContent(atQ_130Te))
	N_at130Te_particles.append(hist_beta_glyptal.GetBinContent(atQ_130Te))
	N_at130Te_particles.append(hist_beta_epoxy.GetBinContent(atQ_130Te))
	N_at130Te_particles.append(hist_gamma_glyptal.GetBinContent(atQ_130Te))
	N_at130Te_particles.append(hist_gamma_epoxy.GetBinContent(atQ_130Te))
	N_at130Te_particles = unp.uarray(N_at130Te_particles, stat_err_at130Te_particles)

# for epoxy
# create combined histogram
	hist_comb_particles_epox  = ROOT.TH1F('hist_comb_particles_epox', '', bins, 0, E_range)
	hist_comb_particles_epox.Add(hist_alpha_epoxy)
	hist_comb_particles_epox.Add(hist_beta_epoxy)
	hist_comb_particles_epox.Add(hist_gamma_epoxy)
	hist_comb_particles_epox.SetLineColor(1)
	hist_comb_particles_epox.SetLineWidth(1)

	# adding all histograms
	hists_particle_sorts_epox.Add(hist_alpha_epoxy)
	hists_particle_sorts_epox.Add(hist_beta_epoxy)
	hists_particle_sorts_epox.Add(hist_gamma_epoxy)
	hists_particle_sorts_epox.Add(hist_comb_particles_epox)
	hists_particle_sorts_epox.Draw("nostack")

	# nice labeling
	hists_particle_sorts_epox.GetXaxis().SetTitleSize(0.03)
	hists_particle_sorts_epox.GetXaxis().SetTitle("E / keV")
	hists_particle_sorts_epox.GetXaxis().SetTickLength(0.02)
	hists_particle_sorts_epox.GetXaxis().SetTicks("-")
	hists_particle_sorts_epox.GetXaxis().SetLabelOffset(-0.05)
	hists_particle_sorts_epox.GetXaxis().SetTitleOffset(-1.4)

	hists_particle_sorts_epox.GetYaxis().SetTitleSize(0.025)
	hists_particle_sorts_epox.GetYaxis().SetTitle('counts/kg keV yr')
	hists_particle_sorts_epox.GetYaxis().SetTickLength(0.02)
	hists_particle_sorts_epox.GetYaxis().SetTicks("+")
	hists_particle_sorts_epox.GetYaxis().SetLabelOffset(-0.025)
	hists_particle_sorts_epox.GetYaxis().SetTitleOffset(-1.8)

	hists.SetMinimum(1e-3)			# Set minimum of y-axis
	hists.SetMaximum(1e4)			# Set maximum of y-axis

	# Legend
	leg = ROOT.TLegend(0.55, 0.8, 0.9, 0.9)
	leg.SetHeader("Particles from epoxy paints")
	leg.SetNColumns(4)
	leg.SetTextSize(0.025)

	leg.AddEntry('hist_alpha_epoxy', "alphas", "f")
	leg.AddEntry('hist_beta_epoxy', "betas", "f")
	leg.AddEntry('hist_gamma_epoxy', "gammas", "f")
	leg.AddEntry('hist_comb_particles_epox', "combined", "f")
	leg.Draw()

	canvas.SetLogy()
	canvas.Update()
	canvas.Print(savename[7]+'.pdf')
	canvas.Clear()
	leg.Clear()



# for glyptal
# create combined histogram
	hist_comb_particles_glyp  = ROOT.TH1F('hist_comb_particles_glyp', '', bins, 0, E_range)
	hist_comb_particles_glyp.Add(hist_alpha_glyptal)
	hist_comb_particles_glyp.Add(hist_beta_glyptal)
	hist_comb_particles_glyp.Add(hist_gamma_glyptal)
	hist_comb_particles_glyp.SetLineColor(1)
	hist_comb_particles_glyp.SetLineWidth(1)

	# adding all histograms
	hists_particle_sorts_glyp.Add(hist_alpha_glyptal)
	hists_particle_sorts_glyp.Add(hist_beta_glyptal)
	hists_particle_sorts_glyp.Add(hist_gamma_glyptal)
	hists_particle_sorts_glyp.Add(hist_comb_particles_glyp)
	hists_particle_sorts_glyp.Draw("nostack")

	# nice labeling
	hists_particle_sorts_glyp.GetXaxis().SetTitleSize(0.03)
	hists_particle_sorts_glyp.GetXaxis().SetTitle("E / keV")
	hists_particle_sorts_glyp.GetXaxis().SetTickLength(0.02)
	hists_particle_sorts_glyp.GetXaxis().SetTicks("-")
	hists_particle_sorts_glyp.GetXaxis().SetLabelOffset(-0.05)
	hists_particle_sorts_glyp.GetXaxis().SetTitleOffset(-1.4)

	hists_particle_sorts_glyp.GetYaxis().SetTitleSize(0.025)
	hists_particle_sorts_glyp.GetYaxis().SetTitle('counts/kg keV yr')
	hists_particle_sorts_glyp.GetYaxis().SetTickLength(0.02)
	hists_particle_sorts_glyp.GetYaxis().SetTicks("+")
	hists_particle_sorts_glyp.GetYaxis().SetLabelOffset(-0.025)
	hists_particle_sorts_glyp.GetYaxis().SetTitleOffset(-1.8)

	hists.SetMinimum(1e-3)			# Set minimum of y-axis
	hists.SetMaximum(1e4)			# Set maximum of y-axis

	# Legend
	leg = ROOT.TLegend(0.55, 0.8, 0.9, 0.9)
	leg.SetHeader("Particles from glyptal paints")
	leg.SetNColumns(4)
	leg.SetTextSize(0.025)

	leg.AddEntry('hist_alpha_glyptal', "alphas", "f")
	leg.AddEntry('hist_beta_glyptal', "betas", "f")
	leg.AddEntry('hist_gamma_glyptal', "gammas", "f")
	leg.AddEntry('hist_comb_particles_glyp', "combined", "f")
	leg.Draw()

	canvas.SetLogy()
	canvas.Update()
	canvas.Print(savename[6]+'.pdf')
	canvas.Clear()
	leg.Clear()


	# print backgroundcontributions:
	# 1: sorted by glyptal-epoxy isotopes (40K_glyptal, 40K_epoxy, 232Th_glyptal,...)
	# 1: sorted by glyptal-epoxy particle sort (alpha_glyptal, ...)
	N_at116Cd = [N_at116Cd, N_at116Cd_particles]
	N_at130Te = [N_at130Te, N_at130Te_particles]

	#print("\n##################################")
	#print("# Background contributions")
	#print("##################################")
	print("116Cd:")
	for i in range(len(N_at116Cd)):
		print(N_at116Cd[i])
	print()
	print("130Te:")
	for i in range(len(N_at130Te)):
		print(N_at130Te[i])
	#### saving
	dir = 'Latex_tab_prep/'
	data = [[N_at116Cd[0][0], N_at116Cd[0][1]], [N_at116Cd[0][2], N_at116Cd[0][3]], [N_at116Cd[0][4], N_at116Cd[0][5]]]
	savename='at116Cd_isotopes.tex'
	transform2latex_tab_2(data, dir, savename)

	data = []
	data = [[N_at116Cd[1][0], N_at116Cd[1][1]], [N_at116Cd[1][2], N_at116Cd[1][3]], [N_at116Cd[1][4], N_at116Cd[1][5]]]
	savename='at116Cd_particles.tex'
	transform2latex_tab_2(data, dir, savename)

	data = []
	data = [[N_at130Te[0][0], N_at130Te[0][1]], [N_at130Te[0][2], N_at130Te[0][3]], [N_at130Te[0][4], N_at130Te[0][5]]]
	savename='at130Te_isotopes.tex'
	transform2latex_tab_2(data, dir, savename)

	data = []
	data = [[N_at130Te[1][0], N_at130Te[1][1]], [N_at130Te[1][2], N_at130Te[1][3]], [N_at130Te[1][4], N_at130Te[1][5]]]
	savename='at130Te_particles.tex'
	transform2latex_tab_2(data, dir, savename)


	

def background_czt(eventfile, datafile, savename):
	################################
	# Unstacked Histogram
	################################

	canvas = ROOT.TCanvas('canv', 'Histogramm')
	leg = ROOT.TLegend(0.6, 0.7, 0.9, 0.9)

	bins = 1000
	E_range = 3000

	data = read_File(datafile)
	print('Reading "%s" ...\n' %(eventfile))

	iso_list = data[:,0]							# get isotope
	N_norm = convert_str2num(data[:,4])				# get norming factors, number convertion necessary

	# entries at Qvalues:
	Q_116Cd = 2813
	Q_130Te = 2527
	N_at116Cd = []
	N_at130Te = []
	N_at116Cd_particles = []
	N_at130Te_particles = []
	stat_err_at116Cd = []
	stat_err_at130Te = []

	atQ_116Cd = np.int(np.around(Q_116Cd*bins/E_range, decimals=0))
	#print(atQ_116Cd)
	atQ_130Te = np.int(np.around(Q_130Te*bins/E_range, decimals=0))

	print('Data:\n%s \n' %(data))
	print('Isotopes:\n%s \n' %(iso_list))
	print('Number of events per kilogram and year: \n %s\n' %(N_norm))
	all_events = []
	particle_id = []

	print('Reading files:')
	for i in range(len(eventfile)):
		print('\n"%s"' %(eventfile[i]))
	# reading all eventdata and: MeV -> keV
		ev_data = np.genfromtxt(eventfile[i], unpack=True)
		events = 1000*ev_data[0] 				# keV
		particle_id_temp = ev_data[2]
		particle_id.append(particle_id_temp) 
		all_events.append(events)
	# for i in range(len(eventfile)):
	#	events = np.genfromtxt(eventfile[i])
	#	events = 1000*events 				# keV
	#	all_events.append(events)

	hists = ROOT.THStack("hists", "")
	hists_no_lim = ROOT.THStack("hists_no_lim", "")

	# ROOT dislikes numpy-objects so it's necessary to cast to python-objects
	# all_events[i]-entries: np.ndarray is casted to list, 
	# all_events[i][j]-entries: np.float64 is casted to float (automatically)
	for i in range(len(all_events)):
		all_events[i] = all_events[i].tolist()
	
	hist_1 = ROOT.TH1F('hist_1', '', bins, 0, E_range)
	bin_width = hist_1.GetBinWidth(1)
	for j in range(len(all_events[0])):
		hist_1.Fill(all_events[0][j])
	stat_err_at116Cd_temp = np.sqrt(hist_1.GetBinContent(atQ_116Cd))
	stat_err_at130Te_temp = np.sqrt(hist_1.GetBinContent(atQ_130Te))
	#hist_1.SetFillColor(2)
	#hist_1.SetFillStyle(1001)
	hist_1.SetLineColor(2)
	hist_1.SetLineWidth(2)

	#bin_width = hist_1.GetBinWidth(1)

	# Scaling
	scale = N_norm[0]/(1e6)		# for 1 million events # divide per bin width for per keV
	stat_err_at116Cd_temp = scale * stat_err_at116Cd_temp
	stat_err_at116Cd.append(stat_err_at116Cd_temp)
	stat_err_at130Te_temp = scale * stat_err_at130Te_temp
	stat_err_at130Te.append(stat_err_at130Te_temp)
	hist_1.Scale(scale)

	hist_2 = ROOT.TH1F('hist_2', '', bins, 0, E_range)
	for j in range(len(all_events[1])):
		hist_2.Fill(all_events[1][j])
	stat_err_at116Cd_temp = np.sqrt(hist_2.GetBinContent(atQ_116Cd))
	stat_err_at130Te_temp = np.sqrt(hist_2.GetBinContent(atQ_130Te))
	#hist_2.SetFillColor(3)
	#hist_2.SetFillStyle(1001)
	hist_2.SetLineColor(3)
	hist_2.SetLineWidth(2)
	# Scaling
	scale = N_norm[1]/(1e6)		# for 1 million events # divide per bin width for per keV
	stat_err_at116Cd_temp = scale * stat_err_at116Cd_temp
	stat_err_at116Cd.append(stat_err_at116Cd_temp)
	stat_err_at130Te_temp = scale * stat_err_at130Te_temp
	stat_err_at130Te.append(stat_err_at130Te_temp)
	hist_2.Scale(scale)

	hist_3 = ROOT.TH1F('hist_3', '', bins, 0, E_range)
	for j in range(len(all_events[2])):
		hist_3.Fill(all_events[2][j])
	stat_err_at116Cd_temp = np.sqrt(hist_3.GetBinContent(atQ_116Cd))
	stat_err_at130Te_temp = np.sqrt(hist_3.GetBinContent(atQ_130Te))
	#hist_3.SetFillColor(4)
	#hist_3.SetFillStyle(1001)
	hist_3.SetLineColor(4)
	hist_3.SetLineWidth(2)
	# Scaling
	scale = N_norm[2]/(1e6)		# for 1 million events # divide per bin width for per keV
	stat_err_at116Cd_temp = scale * stat_err_at116Cd_temp
	stat_err_at116Cd.append(stat_err_at116Cd_temp)
	stat_err_at130Te_temp = scale * stat_err_at130Te_temp
	stat_err_at130Te.append(stat_err_at130Te_temp)
	hist_3.Scale(scale)

	hist_4 = ROOT.TH1F('hist_4', '', bins, 0, E_range)
	for j in range(len(all_events[3])):
		hist_4.Fill(all_events[3][j])
	stat_err_at116Cd_temp = np.sqrt(hist_4.GetBinContent(atQ_116Cd))
	stat_err_at130Te_temp = np.sqrt(hist_4.GetBinContent(atQ_130Te))
	#hist_4.SetFillColor(6)
	#hist_4.SetFillStyle(1001)
	hist_4.SetLineColor(6)
	hist_4.SetLineWidth(2)
	# Scaling
	scale = N_norm[3]/(1e6)		# for 1 million events # divide per bin width for per keV
	stat_err_at116Cd_temp = scale * stat_err_at116Cd_temp
	stat_err_at116Cd.append(stat_err_at116Cd_temp)
	stat_err_at130Te_temp = scale * stat_err_at130Te_temp
	stat_err_at130Te.append(stat_err_at130Te_temp)
	hist_4.Scale(scale)

	hist_5 = ROOT.TH1F('hist_5', '', bins, 0, E_range)
	for j in range(len(all_events[4])):
		hist_5.Fill(all_events[4][j])
	stat_err_at116Cd_temp = np.sqrt(hist_5.GetBinContent(atQ_116Cd))
	stat_err_at130Te_temp = np.sqrt(hist_5.GetBinContent(atQ_130Te))
	#hist_5.SetFillColor(7)
	#hist_5.SetFillStyle(1001)
	hist_5.SetLineColor(7)
	hist_5.SetLineWidth(2)
	# Scaling
	scale = N_norm[4]/(1e6)		# for 1 million events # divide per bin width for per keV
	stat_err_at116Cd_temp = scale * stat_err_at116Cd_temp
	stat_err_at116Cd.append(stat_err_at116Cd_temp)
	stat_err_at130Te_temp = scale * stat_err_at130Te_temp
	stat_err_at130Te.append(stat_err_at130Te_temp)
	hist_5.Scale(scale)


	################################
	# include upper limits
	################################
	# combine all histograms
	hist_comb_1  = ROOT.TH1F('hist_comb_1', '', bins, 0, E_range)
	hist_comb_1.Add(hist_comb_1, hist_1)
	hist_comb_1.Add(hist_comb_1, hist_2)
	hist_comb_1.Add(hist_comb_1, hist_3)
	hist_comb_1.Add(hist_comb_1, hist_4)
	hist_comb_1.Add(hist_comb_1, hist_5)
	hist_comb_1.SetLineColor(1)
	hist_comb_1.SetLineWidth(1)

	
	hists.Add(hist_1)
	hists.Add(hist_2)
	hists.Add(hist_3)
	hists.Add(hist_4)
	hists.Add(hist_5)
	hists.Add(hist_comb_1)

	hists.Draw("nostack")

	hists.GetXaxis().SetTitleSize(0.03)
	hists.GetXaxis().SetTitle("E / keV")
	hists.GetXaxis().SetTickLength(0.02)
	hists.GetXaxis().SetTicks("-")
	hists.GetXaxis().SetLabelOffset(-0.05)
	hists.GetXaxis().SetTitleOffset(-1.4)

	hists.GetYaxis().SetTitleSize(0.025)
	hists.GetYaxis().SetTitle('counts/kg keV yr')
	hists.GetYaxis().SetTickLength(0.02)
	hists.GetYaxis().SetTicks("+")
	hists.GetYaxis().SetLabelOffset(-0.025)
	hists.GetYaxis().SetTitleOffset(-1.8)

	hists.SetMinimum(1e-4)			# Set minimum of y-axis
	hists.SetMaximum(1e6)			# Set minimum of y-axis


	leg.SetHeader("Isotopes")
	leg.SetNColumns(2)
	leg.SetTextSize(0.025)
	for i in range(0, 5, 1):
		leg.AddEntry('hist_%i' %(i+1), iso_list[i], "f")
	leg.AddEntry('hist_comb_1', "combined", "f")
	leg.Draw()

	canvas.SetLogy()
	canvas.Update()
	canvas.Print(savename+'_unstacked_root_with_upper_limits.pdf')
	canvas.Clear()
	leg.Clear()


	################################
	# without limits
	################################
	# combine all histograms
	hist_comb_2  = ROOT.TH1F('hist_comb_2', '', bins, 0, E_range)
	hist_comb_2.Add(hist_comb_2, hist_2)
	hist_comb_2.Add(hist_comb_2, hist_4)
	hist_comb_2.Add(hist_comb_2, hist_5)
	hist_comb_2.SetLineColor(1)
	hist_comb_2.SetLineWidth(1)
	

	hists_no_lim.Add(hist_2)
	hists_no_lim.Add(hist_4)
	hists_no_lim.Add(hist_5)
	hists_no_lim.Add(hist_comb_2)

	hists_no_lim.Draw("nostack")

	hists_no_lim.GetXaxis().SetTitleSize(0.03)
	hists_no_lim.GetXaxis().SetTitle("E / keV")
	hists_no_lim.GetXaxis().SetTickLength(0.02)
	hists_no_lim.GetXaxis().SetTicks("-")
	hists_no_lim.GetXaxis().SetLabelOffset(-0.05)
	hists_no_lim.GetXaxis().SetTitleOffset(-1.4)

	hists_no_lim.GetYaxis().SetTitleSize(0.025)
	hists_no_lim.GetYaxis().SetTitle('counts/kg keV yr')
	hists_no_lim.GetYaxis().SetTickLength(0.02)
	hists_no_lim.GetYaxis().SetTicks("+")
	hists_no_lim.GetYaxis().SetLabelOffset(-0.025)
	hists_no_lim.GetYaxis().SetTitleOffset(-1.8)

	hists_no_lim.SetMinimum(1e-4)			# Set minimum of y-axis
	hists_no_lim.SetMaximum(1e6)			# Set minimum of y-axis


	leg.SetHeader("Isotopes")
	leg.SetNColumns(2)
	leg.SetTextSize(0.025)
	leg.AddEntry('hist_2', iso_list[1], "f")
	leg.AddEntry('hist_4', iso_list[3], "f")
	leg.AddEntry('hist_5', iso_list[4], "f")
	leg.AddEntry('hist_comb_2', "combined", "f")
	leg.Draw()

	canvas.SetLogy()
	canvas.Update()
	canvas.Print(savename+'_unstacked_root_no_limits.pdf')

	# Getting number of entries at Q-Values
	N_at116Cd.append(hist_1.GetBinContent(atQ_116Cd))
	N_at116Cd.append(hist_2.GetBinContent(atQ_116Cd))
	N_at116Cd.append(hist_3.GetBinContent(atQ_116Cd))
	N_at116Cd.append(hist_4.GetBinContent(atQ_116Cd))
	N_at116Cd.append(hist_5.GetBinContent(atQ_116Cd))
	N_at116Cd = unp.uarray(N_at116Cd, stat_err_at116Cd)
	N_at130Te.append(hist_1.GetBinContent(atQ_130Te))
	N_at130Te.append(hist_2.GetBinContent(atQ_130Te))
	N_at130Te.append(hist_3.GetBinContent(atQ_130Te))
	N_at130Te.append(hist_4.GetBinContent(atQ_130Te))
	N_at130Te.append(hist_5.GetBinContent(atQ_130Te))
	N_at130Te = unp.uarray(N_at130Te, stat_err_at130Te)

	#### saving
	dir = 'Latex_tab_prep/'
	data = [N_at116Cd, N_at130Te]
	savename='doublebetaminus_contribution_at_qval.tex'
	transform2latex_tab_2(data, dir, savename)



################
# run program
################
# background_coating
################

datafile = './Lists/activities_detpaint.txt'
eventfile = ['../rootfiles/edep_entries/entries_40K_glyptal.txt', 
				'../rootfiles/edep_entries/entries_40K_epoxy.txt', 
				'../rootfiles/edep_entries/entries_232Th_glyptal.txt',
				'../rootfiles/edep_entries/entries_232Th_epoxy.txt',
				'../rootfiles/edep_entries/entries_238U_glyptal.txt',
				'../rootfiles/edep_entries/entries_238U_epoxy.txt'
				]


if not os.path.exists('Plots'):
	os.makedirs('Plots')
savename = ['./Plots/paint_background_single',
				'./Plots/paint_background_glyptal',
				'./Plots/paint_background_epoxy',
				'./Plots/paint_background_isotopes',
				'./Plots/paint_background_difference',
				'./Plots/paint_background_particle_sorts',
				'./Plots/paint_background_particle_glyptal',
				'./Plots/paint_background_particle_epoxy'
				]

background_coating(eventfile, datafile, savename)


################
# background CZT
################
datafile = './calc_solutions/calculated_events.txt'
eventfile = ['../rootfiles/edep_entries/entries_114Cd.txt', 
				'../rootfiles/edep_entries/entries_116Cd.txt', 
				'../rootfiles/edep_entries/entries_70Zn.txt',
				'../rootfiles/edep_entries/entries_128Te.txt',
				'../rootfiles/edep_entries/entries_130Te.txt']

savename = './Plots/CZT_spectra'

#background_czt(eventfile, datafile, savename)