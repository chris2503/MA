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

from ROOT import kOrange
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

# the optional parameter k is optional if you would like to view the whole contributions e.g. for czt background
def create_sector_hists(ev_data, scale, k=None, Q_val_returns=None):
	################################
	# creates a histogram for every single sector of all nine detectors
	################################
	bins = 1000
	E_range = 10000

	# entries at Qvalues:
	if Q_val_returns:
		Q_116Cd = 2813
		Q_130Te = 2527
		contrib_at116Cd = []
		contrib_at116Cd_err = []
		contrib_at130Te = []
		contrib_at130Te_err = []


		atQ_116Cd = np.int(np.around(Q_116Cd*bins/E_range, decimals=0))
		atQ_130Te = np.int(np.around(Q_130Te*bins/E_range, decimals=0))

	all_hists = []
	# reading all eventdata and: MeV -> keV
	events = 1000*ev_data[0] 				# keV
	crystal_id = ev_data[1]
	sector_id = ev_data[2]
	particle_id = ev_data[4]
	all_events = events

	subdets = 4
	dets = 9

	if k:
		for i in range(dets):			# for all 9 detectors
			i_h = i+1
			t_d_contrib_at116Cd =[]
			t_d_contrib_at130Te = []
			t_d_contrib_at116Cd_err = []
			t_d_contrib_at130Te_err = []

			for j in range(subdets):	# for all 4 sectors
				# hist_11, hist_12, hist_13, hist_14, hist_21, ..., hist_94
				j_h = j+1
				t_sd_contrib_at116Cd = []
				t_sd_contrib_at130Te = []
				t_sd_contrib_at116Cd_err = []
				t_sd_contrib_at130Te_err = []

				locals()['hist_%s_%i%i' %(k,i_h,j_h)] =  ROOT.TH1F(('hist_%s_%i%i' %(k,i_h,j_h)), ('Detector %i, Sector %i' %(i_h,j_h)), bins, 0, E_range)
				for i_entry in range(len(all_events)):
					if crystal_id[i_entry] == i_h and sector_id[i_entry] == j_h:
						locals()['hist_%s_%i%i' %(k,i_h,j_h)].Fill(all_events[i_entry])
				if Q_val_returns:
					temp_contrib_at116Cd = locals()['hist_%s_%i%i' %(k,i_h,j_h)].GetBinContent(Q_116Cd)
					temp_contrib_at130Te = locals()['hist_%s_%i%i' %(k,i_h,j_h)].GetBinContent(Q_130Te)
					temp_contrib_at116Cd_err = 1/np.sqrt(contrib_at116Cd) * scale
					temp_contrib_at130Te_err = 1/np.sqrt(contrib_at130Cd) * scale

					t_sd_contrib_at116Cd.append(temp_contrib_at116Cd)
					t_sd_contrib_at130Te.append(temp_contrib_at130Te)
					t_sd_contrib_at116Cd_err.append(temp_contrib_at116Cd_err)
					t_sd_contrib_at130Te_err.append(temp_contrib_at130Te_err)

				locals()['hist_%s_%i%i' %(k,i_h,j_h)].Scale(scale)
				locals()['hist_%s_%i%i' %(k,i_h,j_h)].SetStats(False)
				all_hists.append(locals()['hist_%s_%i%i' %(k,i_h,j_h)])
			if Q_val_returns:
				t_d_contrib_at116Cd.append(t_sd_contrib_at116Cd)
				t_d_contrib_at130Te.append(t_sd_contrib_at130Te)
				t_d_contrib_at116Cd_err.append(t_sd_contrib_at116Cd_err)
				t_d_contrib_at130Te_err.append(t_sd_contrib_at130Te_err)
		if Q_val_returns:
			contrib_at116Cd.append(t_d_contrib_at116Cd)
			contrib_at130Te.append(t_d_contrib_at130Te)
			contrib_at116Cd_err.append(t_d_contrib_at116Cd_err)
			contrib_at130Te_err.append(t_d_contrib_at130Te_err)


	else:
		for i in range(dets):			# for all 9 detectors
			i_h = i+1
			for j in range(subdets):	# for all 4 sectors
				# hist_11, hist_12, hist_13, hist_14, hist_21, ..., hist_94
				j_h = j+1

				locals()['hist_%i%i' %(i_h,j_h)] =  ROOT.TH1F(('hist_%i%i' %(i_h,j_h)), ('Detector %i, Sector %i' %(i_h,j_h)), bins, 0, E_range)
				for i_entry in range(len(all_events)):
					if crystal_id[i_entry] == i_h and sector_id[i_entry] == j_h:
						locals()['hist_%i%i' %(i_h,j_h)].Fill(all_events[i_entry])
				if Q_val_returns:
					temp_contrib_at116Cd = locals()['hist_%i%i' %(i_h,j_h)].GetBinContent(Q_116Cd)
					temp_contrib_at130Te = locals()['hist_%i%i' %(i_h,j_h)].GetBinContent(Q_130Te)
					temp_contrib_at116Cd_err = 1/np.sqrt(contrib_at116Cd)
					temp_contrib_at130Te_err = 1/np.sqrt(contrib_at130Te)
				locals()['hist_%i%i' %(i_h,j_h)].Scale(scale)
				locals()['hist_%i%i' %(i_h,j_h)].SetStats(False)
				all_hists.append(locals()['hist_%i%i' %(i_h,j_h)])
	if Q_val_returns:
		#print(Q_116Cd[10])
		return all_hists, contrib_at116Cd, contrib_at116Cd_err, contrib_at130Te, contrib_at130Te_err
	else:
		return all_hists


def create_sumsecHist(hists, hcolor=None, k=None):
	dets = 9
	subdets = 4
	bins = 1000
	E_range = 10000
	all_sumhists = []

	if k:
		for i in range(dets):			# for all 9 detectors
			i_h = i+1
			locals()['sumhist_%s_%i' %(k,i_h)] = ROOT.TH1F("sumhist_%s_%i" %(k,i_h), "", bins, 0, E_range)

			for j in range(subdets):	# for all 4 sectors
				# hist_11, hist_12, hist_13, hist_14, hist_21, ..., hist_94
				j_h = j+1
				locals()['hist_%s_%i%i' %(k, i_h,j_h)] = get_sectorHist(hists, i, j)
				locals()['sumhist_%s_%i' %(k, i_h)].Add(locals()['hist_%s_%i%i' %(k,i_h,j_h)])
				locals()['sumhist_%s_%i' %(k, i_h)].SetLineWidth(1)
				if hcolor:
					if i_h == 9:
						locals()['sumhist_%s_%i' %(k,i_h)].SetLineColor(kOrange-3)
					else:
						locals()['sumhist_%s_%i' %(k,i_h)].SetLineColor(i_h+1)

				else:
					locals()['sumhist_%s_%i' %(k, i_h)].SetLineColor(1)
			all_sumhists.append(locals()['sumhist_%s_%i' %(k,i_h)])

	else:
		for i in range(dets):			# for all 9 detectors
			i_h = i+1
			locals()['sumhist_%i' %(i_h)] = ROOT.TH1F("sumhist_%i" %(i_h), "", bins, 0, E_range)

			for j in range(subdets):	# for all 4 sectors
				# hist_11, hist_12, hist_13, hist_14, hist_21, ..., hist_94
				j_h = j+1
				locals()['hist_%i%i' %(i_h,j_h)] = get_sectorHist(hists, i, j)
				locals()['sumhist_%i' %(i_h)].Add(locals()['hist_%i%i' %(i_h,j_h)])
				locals()['sumhist_%i' %(i_h)].SetLineWidth(1)
				if hcolor:
					if i_h == 9:
						locals()['sumhist_%i' %(i_h)].SetLineColor(kOrange-3)
					else:
						locals()['sumhist_%i' %(i_h)].SetLineColor(i_h+1)

				else:
					locals()['sumhist_%i' %(i_h)].SetLineColor(1)
			all_sumhists.append(locals()['sumhist_%i' %(i_h)])
	return all_sumhists


def create_sumdetHist(hists, hcolor=None, k=None):
	dets = 9
	subdets = 4
	bins = 1000
	E_range = 10000
	all_sumhists = []

	if k:
		locals()['sumdethist_%s' %(k)] = ROOT.TH1F("sumdethist_%s" %(k) , "", bins, 0, E_range)
		locals()['sumdethist_%s' %(k)].SetLineWidth(2)
		for i in range(dets):			# for all 9 detectors
			i_h = i+1

			locals()['sumhist_%s_%i' %(k,i_h)] = get_detHist(hists, i)
			locals()['sumhist_%s_%i' %(k,i_h)].SetLineWidth(1)
			locals()['sumdethist_%s' %(k)].Add(locals()['sumhist_%s_%i' %(k, i_h)])
			if i_h == 9:
				locals()['sumhist_%s_%i' %(k,i_h)].SetLineColor(kOrange-3)
			else:
				locals()['sumhist_%s_%i' %(k,i_h)].SetLineColor(i_h+1)
			locals()['sumdethist_%s' %(k)].SetLineColor(1)
		all_sumhists.append(locals()['sumdethist_%s' %(k)])

	else:
		locals()['sumdethist'] = ROOT.TH1F("sumdethist" , "", bins, 0, E_range)
		locals()['sumdethist'].SetLineWidth(1)
		for i in range(dets):			# for all 9 detectors
			i_h = i+1

			locals()['sumhist_%i' %(i_h)] = get_detHist(hists, i)
			locals()['sumhist_%i' %(i_h)].SetLineWidth(1)
			locals()['sumdethist'].Add(locals()['sumhist_%i' %(i_h)])
			if i_h == 9:
				locals()['sumhist_%i' %(i_h)].SetLineColor(kOrange-3)
			else:
				locals()['sumhist_%i' %(i_h)].SetLineColor(i_h+1)
			locals()['sumdethist'].SetLineColor(1)
		all_sumhists.append(locals()['sumdethist'])
	return all_sumhists


def create_det_hists(hists, k=None):
	################################
	# creates a histogram for all nine detectors
	################################
	all_hists = []
	subdets = 4
	dets = 9
	bins = 1000
	E_range = 10000

	if k:
		for i in range(dets):			# for all 9 detectors
			i_h = i+1
			locals()['hist_%s_%i' %(k,i_h)] = ROOT.TH1F(("hist_%s_%i" %(k,i_h)), "Detector %i" %(i_h), bins, 0, E_range)
			for j in range(subdets):	# for all 4 sectors
				# hist_11, hist_12, hist_13, hist_14, hist_21, ..., hist_94
				j_h = j+1

				locals()['hist_%s_%i%i' %(k,i_h,j_h)] = get_sectorHist(hists, i, j)
				locals()['hist_%s_%i%i' %(k,i_h,j_h)].SetStats(False)
				locals()['hist_%s_%i%i' %(k,i_h,j_h)].SetLineWidth(1)
				locals()['hist_%s_%i%i' %(k,i_h,j_h)].SetLineColor(j+2)
				locals()['hist_%s_%i%i' %(k,i_h)].Add(locals()['hist_%s_%i%i' %(k,i_h,j_h)])
			all_hists.append(locals()['hist_%s_%i' %(k,i_h)])

	else:
		for i in range(dets):			# for all 9 detectors
			i_h = i+1
			locals()['hist_%i' %(i_h)] = ROOT.THStack(("hist_%i" %(i_h)), "Detector %i" %(i_h))
			for j in range(subdets):	# for all 4 sectors
				# hist_11, hist_12, hist_13, hist_14, hist_21, ..., hist_94
				j_h = j+1

				locals()['hist_%i%i' %(i_h,j_h)] = get_sectorHist(hists, i, j)
				locals()['hist_%i%i' %(i_h,j_h)].SetStats(False)
				locals()['hist_%i%i' %(i_h,j_h)].SetLineWidth(1)
				locals()['hist_%i%i' %(i_h,j_h)].SetLineColor(j+2)
				locals()['hist_%i' %(i_h)].Add(locals()['hist_%i%i' %(i_h,j_h)])
			all_hists.append(locals()['hist_%i' %(i_h)])
	return all_hists


def create_iso_hist(hists, eventfile, k=None):
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

	if k:
		locals()['hist_%s_%s' %(k,sim)] = ROOT.TH1F("hist_%s_%s" %(k,sim), "", bins, 0, E_range)
		for j in range(dets):	# for all 4 sectors
			# hist_11, hist_12, hist_13, hist_14, hist_21, ..., hist_94
			j_h = j+1

			locals()['hist_%s_%i' %(k,j_h)] = get_detHist(hists, j)
			locals()['hist_%s_%i' %(k,j_h)].SetStats(False)
			locals()['hist_%s_%i' %(k,j_h)].SetLineWidth(1)
			locals()['hist_%s_%i' %(k,j_h)].SetLineColor(j+2)
			locals()['hist_%s_%s' %(k,sim)].Add(locals()['hist_%s_%i' %(k,j_h)])
			hist.append(locals()['hist_%s_%s' %(k,sim)])

	else:
		locals()['hist_%s' %(sim)] = ROOT.THStack("hist_%s" %(sim), "")
		for j in range(dets):	# for all 4 sectors
			# hist_11, hist_12, hist_13, hist_14, hist_21, ..., hist_94
			j_h = j+1

			locals()['hist_%i' %(j_h)] = get_detHist(hists, j)
			locals()['hist_%i' %(j_h)].SetLineWidth(1)
			if j == 9:
				locals()['hist_%i' %(j_h)].SetLineColor(kOrange-3)
			else:
				locals()['hist_%i' %(j_h)].SetLineColor(j_h+1)
			locals()['hist_%s' %(sim)].Add(locals()['hist_%i' %(j_h)])
			hist.append(locals()['hist_%s' %(sim)])
	return hist

def create_material_hist(hist, b):
#	sim = eventfile.split('/')
#	sim = sim[len(sim)-1]
#	sim = sim.split('entries_')
#	sim = sim[len(sim)-1]
#	sim = sim.split('.')
#	sim = sim[0]
	mathist= []
	locals()['summathist_%s' %(b)] = ROOT.THStack("summathist_%s" %(b), "")
	for i in range(len(hist)):
		locals()['sumdethist_%i_%s' %(i,b)] = get_isoHist(hist[i])
		locals()['sumdethist_%i_%s' %(i,b)].SetLineWidth(1)
		locals()['sumdethist_%i_%s' %(i,b)].SetLineColor(i+2)
		locals()['summathist_%s' %(b)].Add(locals()['sumdethist_%i_%s' %(i,b)])
	mathist.append(locals()['summathist_%s' %(b)])
	return mathist

def create_summaterial_hist(hist, b):
	bins = 1000
	E_range = 10000

	summathist= []
	locals()['allmathist_%s' %(b)] = ROOT.TH1F("allmathist_%s" %(b), "", bins, 0, E_range)
	for i in range(len(hist)):
		locals()['sumdethist_%i_%s' %(i,b)] = get_isoHist(hist[i])
		locals()['sumdethist_%i_%s' %(i,b)].SetLineWidth(1)
		locals()['allmathist_%s' %(b)].Add(locals()['sumdethist_%i_%s' %(i,b)])
	summathist.append(locals()['allmathist_%s' %(b)])
	return summathist

def create_sumHist(hists, b):
	allhist=[]
	locals()['all_hists'] = ROOT.THStack("all_hists", "")
	for i in range(len(hists)):
		locals()['summathist_%s' %(b)] = get_isoHist(hists[i])
		locals()['summathist_%s' %(b)].SetLineWidth(1)
		locals()['summathist_%s' %(b)].SetLineColor(i+2)
		locals()['all_hists'].Add(locals()['summathist_%s' %(b)])
	allhist.append(locals()['all_hists'])
	return allhist

def create_allsumHist(hists, b):
	bins = 1000
	E_range = 10000

	sumhist= []
	locals()['all_sumhists'] = ROOT.TH1F("all_sumhists", "", bins, 0, E_range)
	for i in range(len(hists)):
		locals()['summathist_%s' %(b)] = get_isoHist(hists[i])
		locals()['summathist_%s' %(b)].SetLineWidth(1)
		locals()['summathist_%s' %(b)].SetLineWidth(1)
		locals()['all_sumhists'].Add(locals()['summathist_%s' %(b)])
	sumhist.append(locals()['all_sumhists'])
	return sumhist


# get the histogram for a sector
def get_sectorHist(hists, i, j):
	return hists[4*i+j]

def get_detHist(hists, i):
	return hists[i]

def get_isoHist(hist):
	return hist[0]


# delete all histograms of each sector
def delete_all_sectorHists(hists, k=None):
	subdets = 4
	dets = 9
	for i in range(dets):			# for all 9 detectors
		i_h = i+1
		for j in range(subdets):	# for all 4 sectors
			j_h = j+1
			if k:
				locals()['hist_%s_%i%i' %(k,i_h,j_h)] = get_sectorHist(hists, i, j)
				locals()['hist_%s_%i%i' %(k,i_h,j_h)].Delete()
			else:
				locals()['hist_%i%i' %(i_h,j_h)] = get_sectorHist(hists, i, j)
				locals()['hist_%i%i' %(i_h,j_h)].Delete()


def delete_all_detHists(hists, k=None):
	dets = 9
	for i in range(dets):			# for all 9 detectors
		i_h = i+1
		if k:
			locals()['hist_%s_%i' %(k,i_h)] = get_detHist(hists, i)
			locals()['hist_%s_%i' %(k,i_h)].Delete()
		else:
			locals()['hist_%i' %(i_h)] = get_detHist(hists, i)
			locals()['hist_%i' %(i_h)].Delete()

def delete_all_sumdetHists(hists, k=None):
	dets = 9
	for i in range(dets):			# for all 9 detectors
		i_h = i+1
		if k:
			locals()['sumhist_%s_%i' %(k,i_h)] = get_detHist(hists, i)
			locals()['sumhist_%s_%i' %(k,i_h)].Delete()
		else:
			locals()['sumhist_%i' %(i_h)] = get_detHist(hists, i)
			locals()['sumhist_%i' %(i_h)].Delete()

def delete_iso_Hist(hist, eventfile, k=None):
	sim = eventfile.split('/')
	sim = sim[len(sim)-1]
	sim = sim.split('entries_')
	sim = sim[len(sim)-1]
	sim = sim.split('.')
	sim = sim[0]

	if k:
		locals()['hist_%s_%s' %(k,sim)] = get_isoHist(hist)
		locals()['hist_%s_%s' %(k,sim)].Delete()

	else:
		locals()['hist_%s' %(sim)] = get_isoHist(hist)
		locals()['hist_%s' %(sim)].Delete()


def delete_iso_sumHist(hist, k=None):
	if k:
		locals()['sumdethist_%s' %(k)] = get_isoHist(hist)
		locals()['sumdethist_%s' %(k)].Delete()

	else:
		locals()['sumdethist'] = get_isoHist(hist)
		locals()['sumdethist'].Delete()

def save_sumHist(hists, sumhist, b):
	canvas = ROOT.TCanvas('canv', 'Histogramm')
	labels=[]
	for i_lab in range(len(b)):
		if b[i_lab] == 'czt':
			labels.append('CdZnTe')
		if b[i_lab] == 'czt_nolim':
			labels.append('CdZnTe')
		elif b[i_lab] == 'coating':
			labels.append('Detector Coatings')


	print('Plotting contributions of all background contributions ...')

	dir = './Plots/'

	locals()['all_hists'] = get_isoHist(hists)
	locals()['all_sumhists'] = get_isoHist(sumhist)
	locals()['all_sumhists'].SetLineColor(1)
	locals()['all_hists'].Add(locals()['all_sumhists'])
	# nice labeling
	locals()['all_hists'].Draw('nostack')
	locals()['all_hists'].GetXaxis().SetTitleSize(0.03)
	locals()['all_hists'].GetXaxis().SetTitle("E in keV")
	locals()['all_hists'].GetXaxis().SetTickLength(0.02)
	locals()['all_hists'].GetXaxis().SetTicks("-")
	locals()['all_hists'].GetXaxis().SetLabelOffset(-0.05)
	locals()['all_hists'].GetXaxis().SetTitleOffset(-1.4)

	locals()['all_hists'].GetYaxis().SetTitleSize(0.025)
	locals()['all_hists'].GetYaxis().SetTitle('N in #frac{counts}{kg keV yr}')
	locals()['all_hists'].GetYaxis().SetTickLength(0.02)
	locals()['all_hists'].GetYaxis().SetTicks("+")
	locals()['all_hists'].GetYaxis().SetLabelOffset(-0.025)
	locals()['all_hists'].GetYaxis().SetTitleOffset(-1.8)

	locals()['all_hists'].GetHistogram().SetStats(False)
	locals()['all_hists'].GetHistogram().SetAxisRange(0., 10000, "X")	# Set maximum of x-axis
	locals()['all_hists'].GetHistogram().SetAxisRange(1e-4, 1e3, "Y")	# Set range of y-axis

	# Legend
	leg = ROOT.TLegend(0.55, 0.8, 0.9, 0.9)
	leg.SetHeader("Isotope contributions")
	leg.SetNColumns(np.int(np.ceil(len(labels)/2)))
	leg.SetTextSize(0.025)

	for i in range(locals()['all_hists'].GetNhists()-1):	# isotopes
		if i == 0:
			locals()['all_hists'].GetHists().First().SetLineColor(2)
			leg.AddEntry(locals()['all_hists'].GetHists().First(), "%s" %(labels[i]), "f")
		else:
			locals()['all_hists'].GetHists().At(i).SetLineColor(i+2)
			leg.AddEntry(locals()['all_hists'].GetHists().At(i), "%s" %(labels[i]), "f")
	leg.AddEntry(locals()['all_sumhists'], "combined", "f")
	leg.Draw()

	canvas.SetLogy(True)
	canvas.Update()

	if not os.path.exists(dir):
		os.makedirs(dir)
	canvas.Print(dir+'whole_contributions_%s.pdf' %(b))
	canvas.Clear()
	print('Plotting successful :) \n')


def save_material_hists(hists, sumhist, b):
	canvas = ROOT.TCanvas('canv', 'Histogramm')
	if b == 'czt':
		labels = ['114Cd', '116Cd', '70Zn', '128Te', '130Te']

	if b == 'czt_nolim':
		labels = ['116Cd', '128Te', '130Te']

	if b == 'coating':
		labels = ['40K_glyp', '40K_epox', '232Th_glyp', '232Th_epox', '238U_glyp', '238U_epox']

	print('Plotting contributions of %s ...\n' %(b))

	dir = './Plots/'

	locals()['summathist_%s' %(b)] = get_isoHist(hists)
	locals()['allmathist_%s' %(b)] = get_isoHist(sumhist)
	locals()['allmathist_%s' %(b)].SetLineColor(1)
	locals()['summathist_%s' %(b)].Add(locals()['allmathist_%s' %(b)])
	# nice labeling
	locals()['summathist_%s' %(b)].Draw('nostack')
	locals()['summathist_%s' %(b)].GetXaxis().SetTitleSize(0.03)
	locals()['summathist_%s' %(b)].GetXaxis().SetTitle("E in keV")
	locals()['summathist_%s' %(b)].GetXaxis().SetTickLength(0.02)
	locals()['summathist_%s' %(b)].GetXaxis().SetTicks("-")
	locals()['summathist_%s' %(b)].GetXaxis().SetLabelOffset(-0.05)
	locals()['summathist_%s' %(b)].GetXaxis().SetTitleOffset(-1.4)

	locals()['summathist_%s' %(b)].GetYaxis().SetTitleSize(0.025)
	locals()['summathist_%s' %(b)].GetYaxis().SetTitle('N in #frac{counts}{kg keV yr}')
	locals()['summathist_%s' %(b)].GetYaxis().SetTickLength(0.02)
	locals()['summathist_%s' %(b)].GetYaxis().SetTicks("+")
	locals()['summathist_%s' %(b)].GetYaxis().SetLabelOffset(-0.025)
	locals()['summathist_%s' %(b)].GetYaxis().SetTitleOffset(-1.8)

	locals()['summathist_%s' %(b)].GetHistogram().SetStats(False)
	locals()['summathist_%s' %(b)].GetHistogram().SetAxisRange(0., 10000, "X")	# Set maximum of x-axis
	locals()['summathist_%s' %(b)].GetHistogram().SetAxisRange(1e-4, 1e3, "Y")	# Set range of y-axis

	# Legend
	leg = ROOT.TLegend(0.55, 0.8, 0.9, 0.9)
	leg.SetHeader("Isotope contributions")
	leg.SetNColumns(np.int(np.ceil(len(labels)/2)))
	leg.SetTextSize(0.025)

	for i in range(locals()['summathist_%s' %(b)].GetNhists()-1):	# isotopes
		if i == 0:
			locals()['summathist_%s' %(b)].GetHists().First().SetLineColor(2)
			leg.AddEntry(locals()['summathist_%s' %(b)].GetHists().First(), "%s" %(labels[i]), "f")
		else:
			locals()['summathist_%s' %(b)].GetHists().At(i).SetLineColor(i+2)
			leg.AddEntry(locals()['summathist_%s' %(b)].GetHists().At(i), "%s" %(labels[i]), "f")
	leg.AddEntry(locals()['allmathist_%s' %(b)], "combined", "f")
	leg.Draw()

	canvas.SetLogy(True)
	canvas.Update()

	if not os.path.exists(dir):
		os.makedirs(dir)
	canvas.Print(dir+'%s.pdf' %(b))
	canvas.Clear()
	print('Plotting successful :) \n')


# save histogram of every single sectors
def save_single_sector_hists(hists, eventfile, x_range, k=None):
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
			if k:
				locals()['hist_%i_%i%i' %(k,i_h,j_h)] = get_sectorHist(hists, i, j)
				locals()['hist_%i_%i%i' %(k,i_h,j_h)].SetLineColor(2)
				locals()['hist_%i_%i%i' %(k,i_h,j_h)].SetLineWidth(1)

				# nice labeling
				locals()['hist_%i_%i%i' %(k,i_h,j_h)].GetXaxis().SetTitleSize(0.03)
				locals()['hist_%i_%i%i' %(k,i_h,j_h)].GetXaxis().SetTitle("E in keV")
				locals()['hist_%i_%i%i' %(k,i_h,j_h)].GetXaxis().SetTickLength(0.02)
				locals()['hist_%i_%i%i' %(k,i_h,j_h)].GetXaxis().SetTicks("-")
				locals()['hist_%i_%i%i' %(k,i_h,j_h)].GetXaxis().SetLabelOffset(-0.05)
				locals()['hist_%i_%i%i' %(k,i_h,j_h)].GetXaxis().SetTitleOffset(-1.4)

				locals()['hist_%i_%i%i' %(k,i_h,j_h)].GetYaxis().SetTitleSize(0.025)
				locals()['hist_%i_%i%i' %(k,i_h,j_h)].GetYaxis().SetTitle('N in #frac{counts}{kg keV yr}')
				locals()['hist_%i_%i%i' %(k,i_h,j_h)].GetYaxis().SetTickLength(0.02)
				locals()['hist_%i_%i%i' %(k,i_h,j_h)].GetYaxis().SetTicks("+")
				locals()['hist_%i_%i%i' %(k,i_h,j_h)].GetYaxis().SetLabelOffset(-0.025)
				locals()['hist_%i_%i%i' %(k,i_h,j_h)].GetYaxis().SetTitleOffset(-1.8)
				locals()['hist_%i_%i%i' %(k,i_h,j_h)].SetStats(False)

				locals()['hist_%i_%i%i' %(k,i_h,j_h)].SetAxisRange(0., x_range, "X")	# Set maximum of x-axis
				locals()['hist_%i_%i%i' %(k,i_h,j_h)].SetAxisRange(1e-4, 1e3, "Y")	# Set range of y-axis
				locals()['hist_%i_%i%i' %(k,i_h,j_h)].Draw()

			else:
				locals()['hist_%i%i' %(i_h,j_h)] = get_sectorHist(hists, i, j)
				locals()['hist_%i%i' %(i_h,j_h)].SetLineColor(2)
				locals()['hist_%i%i' %(i_h,j_h)].SetLineWidth(1)

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
				locals()['hist_%i%i' %(i_h,j_h)].Draw()


			canvas.SetLogy(True)
			canvas.Update()

			if not os.path.exists(dir):
				os.makedirs(dir)
			canvas.Print('%shist_%i%i.pdf' %(dir,i_h,j_h))
			canvas.Clear()
	print('Plotting successful :) \n')


def save_single_det_hists(hists, mysumhist, eventfile, x_range, k=None):
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
	if k:
		for i in range(dets):			# for all 9 detectors
			i_h = i+1
			locals()['hist_%i_%i' %(k,i_h)] = get_detHist(hists, i)

			# nice labeling
			locals()['sumhist_%i_%i' %(k,i_h)] = get_detHist(mysumhist, i)
			locals()['hist_%i_%i' %(k,i_h)].Add(locals()['sumhist_%i' %(i_h)])
			locals()['hist_%i_%i' %(k,i_h)].Draw('nostack')
			locals()['hist_%i_%i' %(k,i_h)].SetLineColor(3)
			locals()['hist_%i_%i' %(k,i_h)].SetLineWidth(1)
			locals()['hist_%i_%i' %(k,i_h)].SetStats(False)
			locals()['hist_%i_%i' %(k,i_h)].GetXaxis().SetTitleSize(0.03)
			locals()['hist_%i_%i' %(k,i_h)].GetXaxis().SetTitle("E in keV")
			locals()['hist_%i_%i' %(k,i_h)].GetXaxis().SetTickLength(0.02)
			locals()['hist_%i_%i' %(k,i_h)].GetXaxis().SetTicks("-")
			locals()['hist_%i_%i' %(k,i_h)].GetXaxis().SetLabelOffset(-0.05)
			locals()['hist_%i_%i' %(k,i_h)].GetXaxis().SetTitleOffset(-1.4)

			locals()['hist_%i_%i' %(k,i_h)].GetYaxis().SetTitleSize(0.025)
			locals()['hist_%i_%i' %(k,i_h)].GetYaxis().SetTitle('N in #frac{counts}{kg keV yr}')
			locals()['hist_%i_%i' %(k,i_h)].GetYaxis().SetTickLength(0.02)
			locals()['hist_%i_%i' %(k,i_h)].GetYaxis().SetTicks("+")
			locals()['hist_%i_%i' %(k,i_h)].GetYaxis().SetLabelOffset(-0.025)
			locals()['hist_%i_%i' %(k,i_h)].GetYaxis().SetTitleOffset(-1.8)
			locals()['hist_%i_%i' %(k,i_h)].SetAxisRange(0., x_range, "X")	# Set maximum of x-axis
			locals()['hist_%i_%i' %(k,i_h)].SetAxisRange(1e-4, 1e3, "Y")	# Set range of y-axis


	else:
		for i in range(dets):			# for all 9 detectors
			i_h = i+1
			locals()['hist_%i' %(i_h)] = get_detHist(hists, i)
			# nice labeling
			locals()['sumhist_%i' %(i_h)] = get_detHist(mysumhist, i)
			locals()['hist_%i' %(i_h)].Add(locals()['sumhist_%i' %(i_h)])
			locals()['hist_%i' %(i_h)].Draw('nostack')
			#locals()['hist_%i' %(i_h)].SetLineColor(3)
			#locals()['hist_%i' %(i_h)].SetLineWidth(1)
			locals()['hist_%i' %(i_h)].GetHistogram().SetStats(False)
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
			locals()['hist_%i' %(i_h)].GetHistogram().SetAxisRange(0., x_range, "X")	# Set maximum of x-axis
			locals()['hist_%i' %(i_h)].GetHistogram().SetAxisRange(1e-4, 1e3, "Y")	# Set range of y-axis

			# Legend
			leg = ROOT.TLegend(0.55, 0.8, 0.9, 0.9)
			leg.SetHeader("Sector")
			leg.SetNColumns(5)
			leg.SetTextSize(0.025)

			if k:
				for j in range(subdets):	# for all 4 sectors
					j_h = j+1
					leg.AddEntry('hist_%i_%i%i' %(k,i_h,j_h), "%i" %(j_h), "f")
				leg.AddEntry(locals()['sumhist_%i_%i' %(k,i_h)], "combined", "f")
			else:
				for j in range(subdets):	# for all 4 sectors
					j_h = j+1
					leg.AddEntry('hist_%i%i' %(i_h,j_h), "%i" %(j_h), "f")
				leg.AddEntry(locals()['sumhist_%i' %(i_h)], "combined", "f")
			leg.Draw()

			canvas.SetLogy(True)
			canvas.Update()

			if not os.path.exists(dir):
				os.makedirs(dir)
			canvas.Print('%shist_%i.pdf' %(dir,i_h))
			canvas.Clear()

def save_single_iso_hists(hist, sumhist, eventfile, x_range, k=None):
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

	if k:
		locals()['hist_%i_%s' %(k,sim)] = get_isoHist(hist)
		locals()['sumdethist_%i' %(k)] = get_isoHist(sumhist)
		locals()['sumdethist_%i' %(k)].SetLineWidth(1)
		locals()['hist_%i_%s' %(k,sim)].Add(locals()['sumdethist'])
		# nice labeling
		locals()['hist_%i_%s' %(k,sim)].Draw('nostack')
		locals()['hist_%i_%s' %(k,sim)].SetStats(False)
		locals()['hist_%i_%s' %(k,sim)].SetLineColor(3)
		locals()['hist_%i_%s' %(k,sim)].SetLineWidth(2)
		locals()['hist_%i_%s' %(k,sim)].GetXaxis().SetTitleSize(0.03)
		locals()['hist_%i_%s' %(k,sim)].GetXaxis().SetTitle("E in keV")
		locals()['hist_%i_%s' %(k,sim)].GetXaxis().SetTickLength(0.02)
		locals()['hist_%i_%s' %(k,sim)].GetXaxis().SetTicks("-")
		locals()['hist_%i_%s' %(k,sim)].GetXaxis().SetLabelOffset(-0.05)
		locals()['hist_%i_%s' %(k,sim)].GetXaxis().SetTitleOffset(-1.4)

		locals()['hist_%i_%s' %(k,sim)].GetYaxis().SetTitleSize(0.025)
		locals()['hist_%i_%s' %(k,sim)].GetYaxis().SetTitle('N in #frac{counts}{kg keV yr}')
		locals()['hist_%i_%s' %(k,sim)].GetYaxis().SetTickLength(0.02)
		locals()['hist_%i_%s' %(k,sim)].GetYaxis().SetTicks("+")
		locals()['hist_%i_%s' %(k,sim)].GetYaxis().SetLabelOffset(-0.025)
		locals()['hist_%i_%s' %(k,sim)].GetYaxis().SetTitleOffset(-1.8)
		locals()['hist_%i_%s' %(k,sim)].SetAxisRange(0., x_range, "X")	# Set maximum of x-axis
		locals()['hist_%i_%s' %(k,sim)].SetAxisRange(1e-4, 1e3, "Y")	# Set range of y-axis

	else:
		locals()['hist_%s' %(sim)] = get_isoHist(hist)
		locals()['sumdethist'] = get_isoHist(sumhist)
		locals()['sumdethist'].SetLineWidth(1)
		locals()['hist_%s' %(sim)].Add(locals()['sumdethist'])
		# nice labeling
		locals()['hist_%s' %(sim)].Draw('nostack')
		#locals()['hist_%s' %(sim)].SetStats(False)
		#locals()['hist_%s' %(sim)].SetLineColor(3)
		#locals()['hist_%s' %(sim)].SetLineWidth(2)
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
		locals()['hist_%s' %(sim)].GetHistogram().SetAxisRange(0., x_range, "X")	# Set maximum of x-axis
		locals()['hist_%s' %(sim)].GetHistogram().SetAxisRange(1e-4, 1e3, "Y")	# Set range of y-axis

	# Legend
	leg = ROOT.TLegend(0.55, 0.8, 0.9, 0.9)
	leg.SetHeader("Detector")
	leg.SetNColumns(5)
	leg.SetTextSize(0.025)

	if k:
		for j in range(dets):	# for all 4 sectors
			j_h = j+1
			leg.AddEntry('sumhist_%i_%i' %(k,j_h), "%i" %(j_h), "f")
		leg.AddEntry('sumdethist_%i' %(k), "combined", "f")
	else:
		for j in range(dets):	# for all 4 sectors
			j_h = j+1
			leg.AddEntry('sumhist_%i' %(j_h), "%i" %(j_h), "f")
		leg.AddEntry('sumdethist', "combined", "f")
	leg.Draw()

	canvas.SetLogy(True)
	canvas.Update()

	if not os.path.exists(dir):
		os.makedirs(dir)
	canvas.Print('%shist_%s.pdf' %(dir, sim))
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
				locals()['hist_%i' %(i_h)].Fill(sector_id[i_entry]-0.5)
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
		locals()['hist_%i' %(i_h)].GetXaxis().Set(4, 0.5, 4.5)

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
		locals()['hist_%i' %(i_h)].GetXaxis().Set(4, 0.5,4.5)

		locals()['hist_%i' %(i_h)].GetYaxis().SetTitleSize(0.025)
		locals()['hist_%i' %(i_h)].GetYaxis().SetTitle('N in #frac{counts}{kg yr}')
		locals()['hist_%i' %(i_h)].GetYaxis().SetTickLength(0.02)
		locals()['hist_%i' %(i_h)].GetYaxis().SetTicks("+")
		locals()['hist_%i' %(i_h)].GetYaxis().SetLabelOffset(-0.025)
		locals()['hist_%i' %(i_h)].GetYaxis().SetTitleOffset(-2.0)


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
	elif background:
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
	locals()['hist_%s' %(sim)].GetXaxis().Set(9, 0.5,9.5)
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
	if background:
		plt.savefig('%sheatmap_%s.pdf' %(dir, background))
	elif eventfile:
		plt.savefig('%sheatmap_%s.pdf' %(dir, sim))
	else:
		print('This case has not been implemented yet!')
