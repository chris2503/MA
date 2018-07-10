import numpy as np
import uncertainties
import uncertainties.unumpy as unp
from uncertainties.unumpy import(nominal_values as noms, std_devs as stds)
import scipy.constants as const
import re
import os
import sys
from time import sleep
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt

import ROOT

sys.path.append(os.path.abspath('./Methods/'))
#print (sys.path)
from Methods.my_methods import get_varname, read_File, convert_str2num, transform2latex_tab_2

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

def read_evData(file):
	print('\nReading file: \t"%s" \n' %(file))
	data = np.genfromtxt(file, unpack=True)
	return data

def create_sector_hists(ev_data, scale):
	################################
	# creates a histogram for every single sector of all nine detectors
	################################
	bins = 1000
	E_range = 10000

	# entries at Qvalues:
	#Q_116Cd = 2813
	#Q_130Te = 2527
	
	all_hists = []
	# reading all eventdata and: MeV -> keV
	events = 1000*ev_data[0] 				# keV
	crystal_id = ev_data[1]
	sector_id = ev_data[2]
	particle_id = ev_data[4]	
	all_events = events

	subdets = 4
	dets = 9

	for i in range(dets):			# for all 9 detectors
		i_h = i+1
		for j in range(subdets):	# for all 4 sectors
			# hist_11, hist_12, hist_13, hist_14, hist_21, ..., hist_94
			j_h = j+1

			locals()['hist_%i%i' %(i_h,j_h)] =  ROOT.TH1F(('hist_%i%i' %(i_h,j_h)), ('Detector %i, Sector %i' %(i_h,j_h)), bins, 0, E_range)
			for i_entry in range(len(all_events)):
				if crystal_id[i_entry] == i_h and sector_id[i_entry] == j_h:
					locals()['hist_%i%i' %(i_h,j_h)].Fill(all_events[i_entry])
			locals()['hist_%i%i' %(i_h,j_h)].Scale(scale)
			locals()['hist_%i%i' %(i_h,j_h)].SetStats(False)
			all_hists.append(locals()['hist_%i%i' %(i_h,j_h)])
	return all_hists


def create_det_hists(hists):
	################################
	# creates a histogram for all nine detectors
	################################
	all_hists = []
	subdets = 4
	dets = 9
	bins = 1000
	E_range = 10000

	for i in range(dets):			# for all 9 detectors
		i_h = i+1
		locals()['hist_%i' %(i_h)] = ROOT.TH1F(("hist_%i" %(i_h)), "Detector %i" %(i_h), bins, 0, E_range)
		for j in range(subdets):	# for all 4 sectors
			# hist_11, hist_12, hist_13, hist_14, hist_21, ..., hist_94
			j_h = j+1
			
			locals()['hist_%i%i' %(i_h,j_h)] = get_sectorHist(hists, i, j)
			locals()['hist_%i%i' %(i_h, j_h)].SetStats(False)
			locals()['hist_%i%i' %(i_h,j_h)].SetLineWidth(2)
			locals()['hist_%i%i' %(i_h,j_h)].SetLineColor(j+2)
			locals()['hist_%i' %(i_h)].Add(locals()['hist_%i%i' %(i_h,j_h)])
		all_hists.append(locals()['hist_%i' %(i_h)])
	return (all_hists)


def create_iso_hist(hists, eventfile):
	dets = 9
	################################
	# creates a histogram for all nine detectors
	################################
	sim = eventfile.split('/')
	sim = sim[len(sim)-1]
	sim = sim.split('entries_')
	sim = sim[len(sim)-1]
	sim = sim.split('.')
	sim = sim[0]

	bins = 1000
	E_range = 10000
	hist = []

	locals()['hist_%s' %(sim)] = ROOT.TH1F("hist_%s" %(sim), "", bins, 0, E_range)
	for j in range(dets):	# for all 4 sectors
		# hist_11, hist_12, hist_13, hist_14, hist_21, ..., hist_94
		j_h = j+1
		
		locals()['hist_%i' %(j_h)] = get_detHist(hists, j)
		locals()['hist_%i' %(j_h)].SetStats(False)
		locals()['hist_%i' %(j_h)].SetLineWidth(2)
		locals()['hist_%i' %(j_h)].SetLineColor(j+2)
		locals()['hist_%s' %(sim)].Add(locals()['hist_%i' %(j_h)])
		hist.append(locals()['hist_%s' %(sim)])
	return hist


# get the histogram for a sector
def get_sectorHist(hists, i, j):
	return hists[4*i+j]

def get_detHist(hists, i):
	return hists[i]

def get_isoHist(hist):
	return hist[0]


# delete all histograms of each sector
def delete_all_sectorHists(hists):
	subdets = 4
	dets = 9
	for i in range(dets):			# for all 9 detectors
		i_h = i+1
		for j in range(subdets):	# for all 4 sectors
			j_h = j+1
			locals()['hist_%i%i' %(i_h,j_h)] = get_sectorHist(hists, i, j)
			locals()['hist_%i%i' %(i_h,j_h)].Delete()

def delete_all_detHists(hists):
	dets = 9
	for i in range(dets):			# for all 9 detectors
		i_h = i+1
		locals()['hist_%i' %(i_h)] = get_detHist(hists, i)
		locals()['hist_%i' %(i_h)].Delete()

def delete_iso_Hist(hist, eventfile):
	sim =  eventfile.split('/')
	sim = sim[len(sim)-1]
	sim = sim.split('entries_')
	sim = sim[len(sim)-1]
	sim = sim.split('.')
	sim = sim[0]
	
	locals()['hist_%s' %(sim)] = get_isoHist(hist)
	locals()['hist_%s' %(sim)].Delete()


# save histogram of every single sectors
def save_single_sector_hists(hists, eventfile, x_range):
	subdets = 4
	dets = 9
	print('Plotting sector histograms ...')
	# creating part of the path for the directory
	sim =  eventfile.split('/')
	sim = sim[len(sim)-1]
	sim = sim.split('entries_')
	sim = sim[len(sim)-1]
	sim = sim.split('.')
	sim = sim[0]

	dir = './Plots/%s/sector/' %(sim)

	canvas = ROOT.TCanvas('canv', 'Histogramm')
	for i in range(dets):			# for all 9 detectors
		i_h = i+1
		for j in range(subdets):	# for all 4 sectors
			# hist_11, hist_12, hist_13, hist_14, hist_21, ..., hist_94
			j_h = j+1
			#print('j_h ', j_h)
			locals()['hist_%i%i' %(i_h,j_h)] = get_sectorHist(hists, i, j)
			locals()['hist_%i%i' %(i_h,j_h)].SetLineColor(2)
			locals()['hist_%i%i' %(i_h,j_h)].SetLineWidth(2)
			
			# nice labeling
			locals()['hist_%i%i' %(i_h,j_h)].GetXaxis().SetTitleSize(0.03)
			locals()['hist_%i%i' %(i_h,j_h)].GetXaxis().SetTitle("E in keV")
			locals()['hist_%i%i' %(i_h,j_h)].GetXaxis().SetTickLength(0.02)
			locals()['hist_%i%i' %(i_h,j_h)].GetXaxis().SetTicks("-")
			locals()['hist_%i%i' %(i_h,j_h)].GetXaxis().SetLabelOffset(-0.05)
			locals()['hist_%i%i' %(i_h,j_h)].GetXaxis().SetTitleOffset(-1.4)

			locals()['hist_%i%i' %(i_h,j_h)].GetYaxis().SetTitleSize(0.025)
			locals()['hist_%i%i' %(i_h,j_h)].GetYaxis().SetTitle('N in #frac{counts}{kg keV yr}')
			locals()['hist_%i%i' %(i_h,j_h)].GetYaxis().SetTickLength(0.02)
			locals()['hist_%i%i' %(i_h,j_h)].GetYaxis().SetTicks("+")
			locals()['hist_%i%i' %(i_h,j_h)].GetYaxis().SetLabelOffset(-0.025)
			locals()['hist_%i%i' %(i_h,j_h)].GetYaxis().SetTitleOffset(-1.8)
			locals()['hist_%i%i' %(i_h,j_h)].SetStats(False)

			locals()['hist_%i%i' %(i_h,j_h)].SetAxisRange(0., x_range, "X")	# Set maximum of x-axis
			locals()['hist_%i%i' %(i_h,j_h)].SetAxisRange(1e-4, 1e3, "Y")	# Set range of y-axis

			locals()['hist_%i%i' %(i_h,j_h)].Draw()


			canvas.SetLogy(True)
			canvas.Update()
			
			if not os.path.exists(dir):
				os.makedirs(dir)
			canvas.Print('%shist_%i%i.pdf' %(dir,i_h,j_h))
			canvas.Clear()
	print('Plotting successful :) \n')


def save_single_iso_hists(hist, eventfile, x_range):
	subdets = 4
	dets = 9
	print('Plotting detector histograms ...')
	sim =  eventfile.split('/')
	sim = sim[len(sim)-1]
	sim = sim.split('entries_')
	sim = sim[len(sim)-1]
	sim = sim.split('.')
	sim = sim[0]

	dir = './Plots/%s/' %(sim)

	canvas = ROOT.TCanvas('canv', 'Histogramm')
	locals()['hist_%s' %(sim)] = get_isoHist(hist)
			
	# nice labeling
	locals()['hist_%s' %(sim)].Draw()
	locals()['hist_%s' %(sim)].SetStats(False)
	locals()['hist_%s' %(sim)].SetLineColor(3)
	locals()['hist_%s' %(sim)].SetLineWidth(2)
	locals()['hist_%s' %(sim)].GetXaxis().SetTitleSize(0.03)
	locals()['hist_%s' %(sim)].GetXaxis().SetTitle("E in keV")
	locals()['hist_%s' %(sim)].GetXaxis().SetTickLength(0.02)
	locals()['hist_%s' %(sim)].GetXaxis().SetTicks("-")
	locals()['hist_%s' %(sim)].GetXaxis().SetLabelOffset(-0.05)
	locals()['hist_%s' %(sim)].GetXaxis().SetTitleOffset(-1.4)

	locals()['hist_%s' %(sim)].GetYaxis().SetTitleSize(0.025)
	locals()['hist_%s' %(sim)].GetYaxis().SetTitle('N in #frac{counts}{kg keV yr}')
	locals()['hist_%s' %(sim)].GetYaxis().SetTickLength(0.02)
	locals()['hist_%s' %(sim)].GetYaxis().SetTicks("+")
	locals()['hist_%s' %(sim)].GetYaxis().SetLabelOffset(-0.025)
	locals()['hist_%s' %(sim)].GetYaxis().SetTitleOffset(-1.8)

	locals()['hist_%s' %(sim)].SetAxisRange(0., x_range, "X")	# Set maximum of x-axis
	locals()['hist_%s' %(sim)].SetAxisRange(1e-4, 1e3, "Y")	# Set range of y-axis

	# Legend
	leg = ROOT.TLegend(0.55, 0.8, 0.9, 0.9)
	leg.SetHeader("Detector")
	leg.SetNColumns(4)
	leg.SetTextSize(0.025)

	for j in range(dets):	# for all 4 sectors
		j_h = j+1
		leg.AddEntry('hist_%i' %(j_h), "%i" %(j_h), "f")
	leg.Draw()

	canvas.SetLogy(True)
	canvas.Update()
			
	if not os.path.exists(dir):
		os.makedirs(dir)
	canvas.Print('%shist_%s.pdf' %(dir, sim))
	canvas.Clear()


def save_single_det_hists(hists, eventfile, x_range):
	dets = 9
	subdets = 4
	print('Plotting detector histograms ...')
	sim =  eventfile.split('/')
	sim = sim[len(sim)-1]
	sim = sim.split('entries_')
	sim = sim[len(sim)-1]
	sim = sim.split('.')
	sim = sim[0]

	dir = './Plots/%s/single_det/' %(sim)

	canvas = ROOT.TCanvas('canv', 'Histogramm')
	for i in range(dets):			# for all 9 detectors
		i_h = i+1
		locals()['hist_%i' %(i_h)] = get_detHist(hists, i)
			
		# nice labeling
		locals()['hist_%i' %(i_h)].Draw()
		locals()['hist_%i' %(i_h)].SetLineColor(3)
		locals()['hist_%i' %(i_h)].SetLineWidth(1)
		locals()['hist_%i' %(i_h)].SetStats(False)
		locals()['hist_%i' %(i_h)].GetXaxis().SetTitleSize(0.03)
		locals()['hist_%i' %(i_h)].GetXaxis().SetTitle("E in keV")
		locals()['hist_%i' %(i_h)].GetXaxis().SetTickLength(0.02)
		locals()['hist_%i' %(i_h)].GetXaxis().SetTicks("-")
		locals()['hist_%i' %(i_h)].GetXaxis().SetLabelOffset(-0.05)
		locals()['hist_%i' %(i_h)].GetXaxis().SetTitleOffset(-1.4)

		locals()['hist_%i' %(i_h)].GetYaxis().SetTitleSize(0.025)
		locals()['hist_%i' %(i_h)].GetYaxis().SetTitle('N in #frac{counts}{kg keV yr}')
		locals()['hist_%i' %(i_h)].GetYaxis().SetTickLength(0.02)
		locals()['hist_%i' %(i_h)].GetYaxis().SetTicks("+")
		locals()['hist_%i' %(i_h)].GetYaxis().SetLabelOffset(-0.025)
		locals()['hist_%i' %(i_h)].GetYaxis().SetTitleOffset(-1.8)

		locals()['hist_%i' %(i_h)].SetAxisRange(0., x_range, "X")	# Set maximum of x-axis
		locals()['hist_%i' %(i_h)].SetAxisRange(1e-4, 1e3, "Y")	# Set range of y-axis

		# Legend
		leg = ROOT.TLegend(0.55, 0.8, 0.9, 0.9)
		leg.SetHeader("Sector")
		leg.SetNColumns(4)
		leg.SetTextSize(0.025)

		for j in range(subdets):	# for all 4 sectors
			j_h = j+1
			leg.AddEntry('hist_%i%i' %(i_h,j_h), "%i" %(j_h), "f")
		leg.Draw()

		canvas.SetLogy(True)
		canvas.Update()
				
		if not os.path.exists(dir):
			os.makedirs(dir)
		canvas.Print('%shist_%i.pdf' %(dir,i_h))
		canvas.Clear()

#########################################################
# histograms of how many events are seen by the detectors
def create_dep_secHist(ev_data, scale):
	bins = 4
	E_range = 4

	# entries at Qvalues:
	#Q_116Cd = 2813
	#Q_130Te = 2527
	
	all_hists = []
	# reading all eventdata and: MeV -> keV
	events = 1000*ev_data[0] 				# keV
	crystal_id = ev_data[1]
	sector_id = ev_data[2]
	particle_id = ev_data[4]	
	all_events = events

	subdets = 4
	dets = 9

	for i in range(dets):			# for all 9 detectors
		i_h = i+1
		locals()['hist_%i' %(i_h)] =  ROOT.TH1F(('hist_%i' %(i_h)), ('Detector %i' %(i_h)), bins, 0, E_range)
		for i_entry in range(len(all_events)):
			if crystal_id[i_entry] == i_h:
				locals()['hist_%i' %(i_h)].Fill(sector_id[i_entry]-1)
		for i_bin in range(bins):
			bin_temp = locals()['hist_%i' %(i_h)].GetBinContent(i_bin+1)
			locals()['hist_%i' %(i_h)].SetBinError(i_bin+1, 1/np.sqrt(bin_temp))
			#print(bin_temp)
			#print(locals()['hist_%i' %(i_h)].GetBinError(i_bin+1))
		locals()['hist_%i' %(i_h)].Scale(scale)
		all_hists.append(locals()['hist_%i' %(i_h)])
	return all_hists


# save the histograms above
def save_dep_sec_hists(hists, eventfile):
	dets = 9
	print('\nPlotting detector histograms ...\n')
	sim =  eventfile.split('/')
	sim = sim[len(sim)-1]
	sim = sim.split('entries_')
	sim = sim[len(sim)-1]
	sim = sim.split('.')
	sim = sim[0]
	y_max = 0

	dir = './Plots/%s/dep/single_det/' %(sim)
	for i in range(dets):
		i_h = i+1
		locals()['hist_%i' %(i_h)] = get_detHist(hists, i)
		y_tmax = locals()['hist_%i' %(i_h)].GetMaximum()
		if y_tmax > y_max:
			y_max = y_tmax
	y_max = 1.2 * y_max


	canvas = ROOT.TCanvas('canv', 'Histogramm')
	for i in range(dets):			# for all 9 detectors
		i_h= i+1

		# nice labeling
		locals()['hist_%i' %(i_h)].Draw()
		locals()['hist_%i' %(i_h)].SetLineColor(3)
		locals()['hist_%i' %(i_h)].SetLineWidth(1)
		locals()['hist_%i' %(i_h)].SetStats(False)
		locals()['hist_%i' %(i_h)].GetXaxis().SetTitleSize(0.03)
		locals()['hist_%i' %(i_h)].GetXaxis().SetTitle("Sector")
		locals()['hist_%i' %(i_h)].GetXaxis().SetTickLength(0.02)
		locals()['hist_%i' %(i_h)].GetXaxis().SetTicks("-")
		locals()['hist_%i' %(i_h)].GetXaxis().SetLabelOffset(-0.05)
		locals()['hist_%i' %(i_h)].GetXaxis().SetTitleOffset(-1.4)
		locals()['hist_%i' %(i_h)].GetXaxis().SetNdivisions(4)

		locals()['hist_%i' %(i_h)].GetYaxis().SetTitleSize(0.025)
		locals()['hist_%i' %(i_h)].GetYaxis().SetTitle('N in #frac{counts}{kg yr}')
		locals()['hist_%i' %(i_h)].GetYaxis().SetTickLength(0.02)
		locals()['hist_%i' %(i_h)].GetYaxis().SetTicks("+")
		locals()['hist_%i' %(i_h)].GetYaxis().SetLabelOffset(-0.025)
		locals()['hist_%i' %(i_h)].GetYaxis().SetTitleOffset(-1.8)

		locals()['hist_%i' %(i_h)].SetAxisRange(0, y_max, "Y")			# Set range of y-axis
		canvas.SetLogy(False)
		canvas.Update()

		if not os.path.exists(dir):
			os.makedirs(dir)
		canvas.Print('%shist_%i.pdf' %(dir,i_h))
	canvas.Clear()
	canvas.Close()

def save_dep_sec_hists_2(hists, eventfile=None, background=None):
	dets = 9
	print('\nPlotting detector histograms ...\n')
	if eventfile:
		sim =  eventfile.split('/')
		sim = sim[len(sim)-1]
		sim = sim.split('entries_')
		sim = sim[len(sim)-1]
		sim = sim.split('.')
		sim = sim[0]
		dir = './Plots/%s/dep/single_det/' %(sim)
	else:
		dir = './Plots/all_dep/'

	canvas = ROOT.TCanvas('canv', 'Histogramm', 500, 500)
	canvas.Divide(3, 3)
	y_max = 0
	for i in range(dets):
		i_h = i+1
		
		locals()['hist_%i' %(i_h)] = get_detHist(hists, i)
		y_tmax = locals()['hist_%i' %(i_h)].GetMaximum()
		if y_tmax > y_max:
			y_max = y_tmax
	y_max = 1.2 * y_max

	for i in range(dets):			# for all 9 detectors
		i_h = i+1
		
		canvas.cd(i_h)

		# nice labeling
		locals()['hist_%i' %(i_h)].Draw()
		locals()['hist_%i' %(i_h)].SetLineColor(3)
		locals()['hist_%i' %(i_h)].SetLineWidth(1)
		locals()['hist_%i' %(i_h)].SetStats(False)
		locals()['hist_%i' %(i_h)].GetXaxis().SetTitleSize(0.03)
		locals()['hist_%i' %(i_h)].GetXaxis().SetTitle("Sector")
		locals()['hist_%i' %(i_h)].GetXaxis().SetTickLength(0.02)
		locals()['hist_%i' %(i_h)].GetXaxis().SetTicks("-")
		locals()['hist_%i' %(i_h)].GetXaxis().SetLabelOffset(-0.05)
		locals()['hist_%i' %(i_h)].GetXaxis().SetTitleOffset(-1.4)
		locals()['hist_%i' %(i_h)].GetXaxis().SetNdivisions(4)

		locals()['hist_%i' %(i_h)].GetYaxis().SetTitleSize(0.025)
		locals()['hist_%i' %(i_h)].GetYaxis().SetTitle('N in #frac{counts}{kg yr}')
		locals()['hist_%i' %(i_h)].GetYaxis().SetTickLength(0.02)
		locals()['hist_%i' %(i_h)].GetYaxis().SetTicks("+")
		locals()['hist_%i' %(i_h)].GetYaxis().SetLabelOffset(-0.025)
		locals()['hist_%i' %(i_h)].GetYaxis().SetTitleOffset(-2.0)
		#locals()['hist_%i' %(i_h)].LabelsOption("u", "Y")


	for i in range(dets):
		i_h = i+1
		locals()['hist_%i' %(i_h)].SetAxisRange(0, y_max, "Y")			# Set range of y-axis
		canvas.cd(i_h)
		canvas.SetLogy(True)
		canvas.Update()

	if not os.path.exists(dir):
		os.makedirs(dir)
	if eventfile:
		canvas.Print('%shist_allSects_%s.pdf' %(dir, sim))
	if background:
		canvas.Print('%shist_allSects_%s.pdf' %(dir, background))
	else:
		print('This case has not been implemented yet!')
	canvas.Clear()
	canvas.Close()


# histograms of how many events are seen by the detectors
def create_dep_detHist(ev_data, eventfile, scale):
	# entries at Qvalues:
	#Q_116Cd = 2813
	#Q_130Te = 2527
	
	hist = []
	# reading all eventdata and: MeV -> keV
	events = 1000*ev_data[0] 				# keV
	crystal_id = ev_data[1]
	sector_id = ev_data[2]
	particle_id = ev_data[4]	
	all_events = events

	sim = eventfile.split('/')
	sim = sim[len(sim)-1]
	sim = sim.split('entries_')
	sim = sim[len(sim)-1]
	sim = sim.split('.')
	sim = sim[0]

	bins = 9
	E_range = 9

	locals()['hist_%s' %(sim)] = ROOT.TH1F("hist_%s" %(sim), "", bins, 0, E_range)
	for i_entry in range(len(all_events)):
		locals()['hist_%s' %(sim)].Fill(crystal_id[i_entry]-1)
	for i_bin in range(bins):
		bin_temp = locals()['hist_%s' %(sim)].GetBinContent(i_bin+1)
		locals()['hist_%s' %(sim)].SetBinError(i_bin+1, 1/np.sqrt(bin_temp))
	locals()['hist_%s' %(sim)].Scale(scale)
	hist.append(locals()['hist_%s' %(sim)])
	return hist


# save the histograms above
def save_dep_det_hists(hist, eventfile):
	print('\nPlotting detector histograms ...\n')
	sim =  eventfile.split('/')
	sim = sim[len(sim)-1]
	sim = sim.split('entries_')
	sim = sim[len(sim)-1]
	sim = sim.split('.')
	sim = sim[0]

	dir = './Plots/%s/dep/single_det/' %(sim)

	canvas = ROOT.TCanvas('canv', 'Histogramm')
	locals()['hist_%s' %(sim)] = get_isoHist(hist)

	# nice labeling
	locals()['hist_%s' %(sim)].Draw()
	locals()['hist_%s' %(sim)].SetLineColor(7)
	locals()['hist_%s' %(sim)].SetLineWidth(1)
	locals()['hist_%s' %(sim)].SetStats(False)
	locals()['hist_%s' %(sim)].GetXaxis().SetTitleSize(0.03)
	locals()['hist_%s' %(sim)].GetXaxis().SetTitle("Detector")
	locals()['hist_%s' %(sim)].GetXaxis().SetTickLength(0.02)
	locals()['hist_%s' %(sim)].GetXaxis().SetTicks("-")
	locals()['hist_%s' %(sim)].GetXaxis().SetLabelOffset(-0.05)
	locals()['hist_%s' %(sim)].GetXaxis().SetTitleOffset(-1.4)
	locals()['hist_%s' %(sim)].GetXaxis().SetNdivisions(9)

	locals()['hist_%s' %(sim)].GetYaxis().SetTitleSize(0.025)
	locals()['hist_%s' %(sim)].GetYaxis().SetTitle('N in #frac{counts}{kg yr}')
	locals()['hist_%s' %(sim)].GetYaxis().SetTickLength(0.02)
	locals()['hist_%s' %(sim)].GetYaxis().SetTicks("+")
	locals()['hist_%s' %(sim)].GetYaxis().SetLabelOffset(-0.025)
	locals()['hist_%s' %(sim)].GetYaxis().SetTitleOffset(-1.8)

	y_max = locals()['hist_%s' %(sim)].GetMaximum()
	y_max = 1.2 * y_max
	locals()['hist_%s' %(sim)].SetAxisRange(0, y_max, "Y")			# Set range of y-axis
	
	canvas.SetLogy(False)
	canvas.Update()

	if not os.path.exists(dir):
		os.makedirs(dir)
	canvas.Print('%shist_%s.pdf' %(dir, sim))
	canvas.Clear()
	canvas.Close()


def save_dep_heatmap(hists, eventfile=None, background=None):
	dets = 9
	subdets = 4
	if eventfile:
		print('\nPlotting heat map ...\n')
		sim =  eventfile.split('/')
		sim = sim[len(sim)-1]
		sim = sim.split('entries_')
		sim = sim[len(sim)-1]
		sim = sim.split('.')
		sim = sim[0]
		dir = './Plots/%s/' %(sim)
	else: 
		dir = './Plots/'

	y_min = 0
	entries = []
	x = []

	for i in range(dets):
		i_h = i+1
		locals()['hist_%i' %(i_h)] = get_detHist(hists, i)
		for j in range(subdets):
			j_h = j+1
			entries.append(locals()['hist_%i' %(i_h)].GetBinContent(j_h))
		x.append(entries)
		entries=[]


	x_transf = np.zeros((9, 2, 2))
	for i in range(9):
		x_transf[i][0][0] = x[i][0]
		x_transf[i][0][1] = x[i][1]
		x_transf[i][1][0] = x[i][2]
		x_transf[i][1][1] = x[i][3]


	sns.set()
	asp = x_transf.shape[0]/float(x_transf.shape[1])
	figw = 5
	figh = 5
	cmap= plt.cm.inferno
	norm = matplotlib.colors.Normalize(vmin= 0, vmax= x_transf.max())

	gridspec_kw = {"height_ratios":[2,2,2], "width_ratios" : [2,2,2]}
	heatmapkws = dict(square=False, cbar=False, linewidths=0.0, cmap=cmap, vmin=0, vmax= x_transf.max() ) 
	tickskw =  dict(xticklabels=False, yticklabels=False)

	left = 0.1; right = 0.8
	bottom = 0.1; top = 0.9
	fig, axes = plt.subplots(ncols=3, nrows=3, figsize=(figw, figh), gridspec_kw=gridspec_kw)
	plt.subplots_adjust(left=left, right=right,bottom=bottom, top=top, wspace=0.1, hspace=0.1 )

	k=0
	for i in range(3):
		for j in range(3):
			sns.heatmap(x_transf[k], ax=axes[i,j], xticklabels=False, yticklabels=False, **heatmapkws)
			k=k+1
	cax = fig.add_axes([0.85,0.1,0.04,0.8])
	axes1 = fig.add_axes([0.09,0.101,0.72,0.8], frameon=False)
	axes1.grid(None)
	axes1.yaxis.set_major_locator(plt.NullLocator())
	axes1.xaxis.set_major_formatter(plt.NullFormatter())
	axes1.set_xlabel('x')
	axes1.set_ylabel('y')
	sm = matplotlib.cm.ScalarMappable(cmap=cmap, norm=norm)
	sm.set_array([])
	clb = fig.colorbar(sm, cax=cax)
	clb.set_label(r'N in $\frac{\mathrm{counts}}{\mathrm{kg\,yr}}$', labelpad=-20, y=1.075, rotation=0)
	if eventfile:
		plt.savefig('%sheatmap_%s.pdf' %(dir, sim))
	if background:
		plt.savefig('%sheatmap_%s.pdf' %(dir, background))
	else:
		print('This case has not been implemented yet!')


##############################################
# main
##############################################
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

#### Coating
m_det_or = np.array([34.21, 35.15, 34.63, 33.91, 35.50, 35.50, 35.50, 35.50, 35.50]) * 1e-3			#kg
m_detpaint_glyptal_or = np.array([0.06, 0.07, 0.06, 0.07]) * 1e-3									#kg
m_detpaint_epoxy_or = np.array([0.14, 0.13, 0.14, 0.12, 0.16]) * 1e-3								#kg

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


#### CZT
datafile = './calc_solutions/calculated_events.txt'
data = read_File(datafile)

iso_list = data[:,0]							# get isotope
N_norm_czt = convert_str2num(data[:,4])			# get norming factors, number convertion necessary
N_simEv = 1e6									# 1 Mio simulated Events

for i in range(len(N_norm_czt)):									
	scale_czt.append(1/4 * N_norm_czt[i] / N_simEv)


def my_main(eventfile, scale, x_range, background):
	# creating plots
	# 1) single plots
	for i_file in range(len(eventfile)):
		data = read_evData(eventfile[i_file])
		my_secHists = create_sector_hists(data, scale[i_file])
		save_single_sector_hists(my_secHists, eventfile[i_file], x_range)
		delete_all_sectorHists(my_secHists)

		my_secHists = create_sector_hists(data, scale[i_file])
		my_detHists = create_det_hists(my_secHists)
		save_single_det_hists(my_detHists, eventfile[i_file], x_range)
		delete_all_detHists(my_detHists)
		delete_all_sectorHists(my_secHists)

		my_secHists = create_sector_hists(data, scale[i_file])
		my_detHists = create_det_hists(my_secHists)

		iso_hist = create_iso_hist(my_detHists, eventfile[i_file])
		save_single_iso_hists(iso_hist, eventfile[i_file], x_range)
		delete_iso_Hist(iso_hist, eventfile[i])
		delete_all_detHists(my_detHists)
		delete_all_sectorHists(my_secHists)

	for i_file in range(len(eventfile)):
		data = read_evData(eventfile[i_file])
		my_depHists = create_dep_secHist(data, scale[i_file])
		save_dep_sec_hists(my_depHists, eventfile[i_file])
		delete_all_detHists(my_depHists)

		my_depHist = create_dep_detHist(data, eventfile[i_file], scale[i_file])
		save_dep_det_hists(my_depHist, eventfile[i_file])
		delete_iso_Hist(my_depHist, eventfile[i])

		my_depHists = create_dep_secHist(data, scale[i_file])
		save_dep_sec_hists_2(my_depHists, eventfile[i_file], None)
		delete_all_detHists(my_depHists)

	for i_file in range(len(eventfile)):
		data = read_evData(eventfile[i_file])
		my_depHists_heat = create_dep_secHist(data, scale[i_file])
		save_dep_heatmap(my_depHists_heat, eventfile[i_file])
		delete_all_detHists(my_depHists_heat)

	# 2) combined plots
	data = []
	for i_file in range(len(eventfile)):
		temp_data = read_evData(eventfile[i_file])
		data.append(temp_data)
		temp_data = []


background = ['coating', 'czt']

my_main(eventfile_coating, scale_coating, x_range_coating, background[0])
my_main(eventfile_czt, scale_czt, x_range_czt, background[1])
