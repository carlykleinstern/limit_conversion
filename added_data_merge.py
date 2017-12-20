#CALL SEQUENCE: added_data_merge.py
#PURPOSE: Convert limits (DD to LHC), UPDATE LUX DATA (no PICO), make plots using matplotlib
#DATE: 20 Dec 2017
#WRITTEN BY: Carly KleinStern, with help from Bjoern Penning and Ryan Wang

#necessary imports
from ROOT import TCanvas, TPad, TFormula, TF1, TPaveLabel, TH1F, TFile, TH2F, TF3, TGraph2D, TRandom, TLegend
from ROOT import gROOT, gBenchmark
from array import array
import sys
import math
import ROOT 
import os
import numpy as np
import matplotlib.pyplot as plt

#define constants
gDM = 1.0
gSM = 0.
mn=0.938
conv_units = 2.568*pow(10.0,27.0) 
Delta_d, Delta_u, Delta_s, = -0.42, 0.85, -0.08
v=246; 
fup, fdp = 0.0208, 0.0411; 
fsp=0.043; 
fTG=1-fup-fdp-fsp

#specify the paths to new datasets
path_n = "/home/vagrant/test/limit_conversion/data/LUX_2016_SD_n.txt"
path_p = "/home/vagrant/test/limit_conversion/data/LUX_2016_SD_p.txt"

data = open("/home/vagrant/notebooks/limit_conversion/axial/input/CMS_monojet_July2017_AXIAL_3.txt", "r")

mMed_lhc_l = []
mDM_lhc_l = []
for line in data:
	
	elems = line.split();
	mMed_lhc = float(elems[0])
	mDM_lhc = float(elems[1])
	mMed_lhc_l.append(float(mMed_lhc))
	mDM_lhc_l.append(float(mDM_lhc))
	print str(mMed_lhc)+" "+str(mDM_lhc)
	
mMed_lhc_a = array('d', mMed_lhc_l)
mDM_lhc_a = array('d', mDM_lhc_l)

#this function yields the name of the dataset (without path or extension)
def base(path):
	base=os.path.basename(path)
	split = os.path.splitext(base)
	name_only = split[0]
	return name_only
	
path1 = base(path_n)
path2 = base(path_p)

#this function does the conversion
def conversion_axial(path):
	mmed_l = []
	mdm_l = [] 
	dataset=open(path)

	for line in dataset:
		elems = line.split();
		mDM = float(elems[0])
		sigma = float(elems[1])*conv_units 
		mu_nDM=mn*mDM/(mn+mDM)
		
		if base(path)=="LUX_2016_SD_n": #change if-else
			gSM = .25
		elif base(path)=="LUX_2016_SD_p":
			gSM = 1.0
		else:
			gSM  = .0
				
		#print base(path)
		f=abs(gDM*(gSM*Delta_u+gSM*Delta_d+gSM*Delta_s)) 
		mMed=pow(f*mu_nDM,0.5)/pow(math.pi*sigma/3.,0.25);
		mmed_l.append(float(mMed))
		mdm_l.append(float(mDM))
		print str(mMed)+" "+str(mDM)

	mmed_a = array('d', mmed_l)
	mdm_a = array('d', mdm_l)
	
	#find range of x and y-axes for plotting
	#print np.amax(mmed_a) #max x axis
	#print np.amin(mmed_a) #min x axis
	#print np.amax(mdm_a) #max y axis
	#print np.amin(mdm_a) #min y axis
	
	#print mmed_a
	#print mdm_a
	
	return mmed_a, mdm_a	
	

array1, array2 = conversion_axial(path_n)
array3, array4 = conversion_axial(path_p)

plt.title("AXIAL VECTOR; LUXSD_n, LUXSD_p & CMS")
plt.plot(array1, array2, 'k-', color='b', label="LUX_2016_SD_n")
plt.plot(array3, array4, 'k-', color='r', label="LUX_2016_SD_p")
plt.plot(mMed_lhc_a, mDM_lhc_a, 'k-', color='g', label="CMS_axial")
plt.yscale("log")
plt.ylabel("mDM")
plt.xlabel("mMed")

plt.grid(True)
plt.legend(loc=4, ncol=1, borderaxespad=0.)

plt.savefig("Axial_LUX_PICO_CMS_update.pdf")


