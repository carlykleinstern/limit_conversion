#CALL SEQUENCE: merge_test8.py
#PURPOSE: Convert limits (DD to LHC), compare plots. Merges axial_p and axial_n interaction limit conversions.
#DATE: 18- Dec 2017
#WRITTEN BY: Carly KleinStern, with help from Bjoern Penning and Ryan Wang

#necessary imports
from ROOT import TCanvas, TPad, TFormula, TF1, TPaveLabel, TH1F, TFile, TH2F, TF3, TGraph2D, TRandom
from ROOT import gROOT, gBenchmark
from array import array
import sys
import math
import ROOT 
import os
import numpy as np

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

#specify the paths of the datasets
path_n = "/home/vagrant/notebooks/limit_conversion/axial/input/LUXSD_n.dat"
path_p = "/home/vagrant/notebooks/limit_conversion/axial/input/PICOSD_p.dat"

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
		
		if base(path)=="LUXSD_n": #change if-else
			gSM = .25
		elif base(path)=="PICOSD_p":
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
	print np.amax(mmed_a) #max x axis
	print np.amin(mmed_a) #min x axis
	print np.amax(mdm_a) #max y axis
	print np.amin(mdm_a) #min y axis
	
	dt = ROOT.TGraph(len(mmed_a), mmed_a, mdm_a)
	
	#graph object
	return dt	

graph1 = conversion_axial(path_n)
graph2 = conversion_axial(path_p)
 
#print graph1
#print graph2

#to make individual plots
#def makeplot(path, graph):
	
	#c1 = TCanvas( 'c1', 'c', 800, 700)	
	
	#if base(path)=="PICOSD_p":
		#c1.SetLogy()
	#else:
		#pass
		
	#graph.SetTitle("DD2LHC " + base(path) +" m_{Med};m_{DM}");
	#graph.Draw("apl");
	#c1.Update()
	#c1.SaveAs("merge_" + base(path) + ".pdf")	

#makeplot(path_n, graph1)
#makeplot(path_p, graph2)

#to compare plots
def superimpose(path_a, path_b, graph_a, graph_b):
	c1 = TCanvas( 'c1', 'c', 800, 700, 1000, 1000)
	c1.SetLogy()
	graph_b.GetYaxis().SetRangeUser(5, 20000)
	graph_b.GetYaxis().SetTitleOffset(1.1)
	graph_b.GetYaxis().SetTitle("mDM")
	graph_b.SetLineColor(1)
	graph_b.SetLineWidth(3)
	graph_b.GetXaxis().SetRangeUser(50, 520)
	graph_b.GetXaxis().SetTitle("mMed")
	graph_b.SetTitle("DD2LHC "+ base(path_a) + " & " + base(path_b));
	graph_b.Draw("apl")
	graph_a.SetLineColor(2)
	graph_a.SetLineWidth(3)
	graph_a.Draw("pl")
	c1.Update()
	c1.SaveAs("merge_" + base(path_a) + " & " + base(path_b) + ".pdf")
	
superimpose(path_n, path_p, graph1, graph2)


