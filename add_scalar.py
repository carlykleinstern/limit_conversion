#CALL SEQUENCE: add_vector.py
#PURPOSE: Convert limits (DD to LHC), ADD SCALAR INTERACTION, 
#create separate plots for axial_vector and vector interactions
#update LUX vector to 2016 data
#WRITTEN BY: Carly KleinStern, with help from Bjoern Penning

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
gSM = 0.0
#all the same now, but could perhaps change these values in the future
gu=gd=gs=1.0
mn=0.938
conv_units = 2.568*pow(10.0,27.0) 
#all the same now, but could perhaps change these values in the future
Delta_d_p, Delta_d_n, Delta_u_p, Delta_u_n, Delta_s_p, Delta_s_n, = -0.42, -0.42, 0.85, 0.85, -0.08, -0.08
v=246; 
fup, fdp = 0.0208, 0.0411; 
fsp=0.043; 
fTG=1-fup-fdp-fsp

#specify the paths to new datasets (relative paths)
path_n = "./DATA/LUX_2016_SD_n.txt"
path_p = "./DATA/LUX_2016_SD_p.txt"
path_v = "./DATA/LUX_2016_SI.txt"
path_s = "./DATA/LUX_2016_SI_scalar.txt"

#get CMS data
data1 = open("./DATA/CMS_monojet_July2017_AXIAL_3.txt", "r")

#vector AND scalar
data2 = open("./DATA/CMS_monojet_July2017_VECTOR.txt", "r")

#this function extracts the CMS data, and puts it into arrays (to be graphed later)
def get_cms_data(data):
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
	
	return mMed_lhc_a, mDM_lhc_a
	
cms_array1, cms_array2 = get_cms_data(data1)
cms_array3, cms_array4 = get_cms_data(data2)

#this function yields the name of the dataset (without path or extension)
def base(path):
	base=os.path.basename(path)
	split = os.path.splitext(base)
	name_only = split[0]
	return name_only
	
path1 = base(path_n)
path2 = base(path_p)
path3 = base(path_v)
path4 = base(path_s)

#this function does the conversion (from DD to LHC)
def conversion(path):
	mmed_l = []
	mdm_l = [] 
	dataset=open(path)

	for line in dataset:
		elems = line.split();
		mDM = float(elems[0])
		sigma = float(elems[1])*conv_units 
		mu_nDM=mn*mDM/(mn+mDM)
		
		if base(path)=="LUX_2016_SD_p" or "LUX_2016_SD_p":
			gSM = 1.0
			
		if base(path)=="LUX_2016_SI" or "LUX1_scalar" or "LUX1_2016_SI_scalar":
			gSM = 1.0
				
		#print base(path)
		if base(path) == "LUX_2016_SD_n" or "LUX_2016_SD_p":
			f=abs(gDM*(gSM*Delta_u_p+gSM*Delta_d_p+gSM*Delta_s_p))
			fsq=f*f #we need to square this quantity
			mMed=pow(fsq*mu_nDM,0.5)/pow(math.pi*sigma/3.,0.25);
			
		if base(path) == "LUX_2016_SI":
			mMed=pow((2*gu+gd)*gDM*mu_nDM,0.5)/pow((math.pi*sigma),0.25);
		
		if base(path) == "LUX_2016_SI_scalar":
			fmMed2=(mn/v)*gSM*gDM*(fup+fdp+fsp+2./27.*fTG*3.);
			mMed=pow(fmMed2*mu_nDM,0.5)/pow(math.pi*sigma,0.25);

		mmed_l.append(float(mMed))
		mdm_l.append(float(mDM))
		#print str(mMed)+" "+str(mDM)

	mmed_a = array('d', mmed_l)
	mdm_a = array('d', mdm_l)
	
	#find range of x and y-axes for plotting
	#print np.amax(mmed_a) #max x axis
	#print np.amin(mmed_a) #min x axis
	#print np.amax(mdm_a) #max y axis
	#print np.amin(mdm_a) #min y axis
	
	print mmed_a
	print mdm_a
	
	return mmed_a, mdm_a	
	
#axial_n arrays
array1, array2 = conversion(path_n)

#axial_p arrays
array3, array4 = conversion(path_p)

#vector arrays
array5, array6 = conversion(path_v)

#scalar arrays
array7, array8 = conversion(path_s)

plt.title("AXIAL VECTOR: LUXSD_n, LUXSD_p & CMS")
plt.plot(array1, array2, 'k-', color='b', label="LUX_2016_SD_n")
plt.plot(array3, array4, 'k-', color='r', label="LUX_2016_SD_p")
plt.plot(cms_array1, cms_array2, 'k-', color='g', label="CMS_2017_axial")
plt.yscale("log")
plt.xscale("log")
plt.ylim((-1000,10000))
plt.xlim((0,2500))
plt.ylabel("mDM [GeV]")
plt.xlabel("mMed [GeV]")
plt.grid(True)
plt.legend(loc=1, ncol=1, borderaxespad=0.0)
#plt.savefig("Axial_LUX_CMS_update2.pdf")
plt.close()

plt.title("VECTOR: LUX & CMS")
plt.plot(array5, array6, 'k-', color="#0165fc", label="LUX_2016_vector")
plt.plot(cms_array3, cms_array4, 'k-', color="#f97306", label="CMS_2017_vector")
plt.ylabel("mDM [GeV]")
plt.xlabel("mMed [GeV]")
plt.yscale("log")
plt.xscale("log")
plt.ylim((0,10E5))
plt.xlim((-10,6.5E4))
plt.grid(True)
plt.legend(loc=1, ncol=1, borderaxespad=0.0)
#plt.savefig("Vector_LUX_CMS_update2.pdf")
plt.close()

plt.title("SCALAR: LUX_2016")
plt.plot(array7, array8, 'k-', color="#f10c45", label="LUX_2016_scalar")
#add LHC_scalar data here
#plt.plot(array, array, 'k-', color="#21fc0d", label="LHC_scalar")
plt.ylabel("mDM [GeV]")
plt.xlabel("mMed [GeV]")
plt.yscale("log")
plt.xscale("log")
plt.ylim((0,10E4))
plt.xlim((-10,3000))
plt.grid(True)
plt.legend(loc=1, ncol=1, borderaxespad=0.0)
plt.savefig("Scalar_LUX_update2016.pdf")




