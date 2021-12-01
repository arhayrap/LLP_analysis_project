import ROOT as rt
# import root_numpy as rtnp
import csv
import re
import sys
import collections
import os
from collections import OrderedDict
import uproot
import pandas as pd
import scipy
import awkward
import numpy as np
import time
import numba
from numba import jit
from matplotlib import pyplot as plt
sys.path.append('/storage/af/user/christiw/gpu/christiw/llp/delayed_jet_analyzer/lib/')
from histo_utilities import create_TH1D, create_TH2D, std_color_list, create_TGraph, make_ratio_plot
import CMS_lumi, tdrstyle
import importlib

a = tdrstyle.setTDRStyle()
CMS_lumi.writeExtraText = 0

wH = 1
Z_MASS = 91.2

# donotdelete = []
print(sys.version)

cut_based = True
cut_based_version='v4'
r_denom = {}
z_denom = {}
r_nom = {}
z_nom = {}
r_nRechits = {}
z_nRechits = {}
weight_nom = {}
weight_denom = {}
weight_nRechits = {}
llp_r = {}
llp_z = {}
weight_llp = {}
decay = 'bbbb'
lumi = [ 35.9, 41.5, 59.7 ]
l = 2

#paths =['/eos/user/a/arhayrap/OUTPUT_DATA/WplusHToSS_STodd_ms15_pl1000.root']
#paths =['/eos/user/a/arhayrap/OUTPUT_DATA/WplusHToSS_SToTauTau_ms15_pl1000.root']

# paths =['/eos/user/a/arhayrap/OUTPUT_DATA/GMSB_STau_mass247_ctau1000.root']
paths =['/eos/user/a/arhayrap/OUTPUT_DATA/GMSB_STau_mass247_ctau3000.root']

path = paths[0]

cut_based = True
cut_based_version='v4'
r_denom = {}
z_denom = {}
r_nom = {}
z_nom = {}
r_nRechits = {}
z_nRechits = {}
weight_nom = {}
weight_denom = {}
weight_nRechits = {}
llp_r = {}
llp_z = {}
weight_llp = {}
decay = 'bbbb'
lumi = [ 35.9, 41.5, 59.7 ]
for m in [ '15', '40','55']:
    r_denom[m] = []
    z_denom[m] = []
    r_nom[m] = []
    z_nom[m] = []
    r_nRechits[m] = []
    z_nRechits[m] = []
    weight_nom[m] = []
    weight_denom[m] = []
    weight_nRechits[m] = []
    llp_r[m] = []
    llp_z[m] =[]
    weight_llp[m] = []
    for i, year in enumerate(["MC_Summer16", "MC_Fall17", "MC_Fall18"]):
        analyzer_version = "/v1/v108/"
        # path = "/storage/cms/store/group/phys_exotica/delayedjets/displacedJetMuonAnalyzer/csc/V1p17/"+year+analyzer_version+"normalized/"
        # path = "/eos/user/a/arhayrap/"
        # root_dir =uproot.open(path+"GMSB_STau_mass247_ctau3000.root")
        root_dir =uproot.open(path)
        T = root_dir['MuonSystem']
        #print("*************/////---->  ", T.array('jetPt'))
        sel_jet = np.logical_and(T.array('jetPt') > 50, np.abs(T.array('jetEta')) < 2.4 )
        #print("#############/////---->  ", sel_jet)
        
        ########### SELECTION: EVENTS ############
        hlt = T['HLTDecision'].array()

        #sel_llp = np.abs(T.array('gLLP_eta'))<3
        #sel_llp = np.logical_and(sel_llp, T.array('gLLP_decay_vertex_r')<800)
        #sel_llp = np.logical_and(sel_llp, np.abs(T.array('gLLP_decay_vertex_z'))<1200)
        #sel_llp = np.logical_and(sel_llp, np.abs(T.array('gLLP_decay_vertex_z'))>200)
        #sel_ev = np.sum(sel_llp,axis = 1) == 1

        sel_ev = np.sum(T.array('gLLP_csc'),axis = 1) == 1
        #sel_ev = np.logical_and(sel_ev, T.array('METNoMuTrigger'))
        #sel_ev = np.logical_and(sel_ev, T.array('category') == 0)
        #sel_ev = np.logical_and(sel_ev, T.array('nLeptons') == 0)
        #sel_ev = np.logical_and(sel_ev, T.array('Flag2_all'))
        #sel_ev = np.logical_and(sel_ev, T.array('metEENoise') >= 200)
        sel_ev = np.logical_and(sel_ev, sel_jet.sum()>=1)

        sel_llp = T.array('gLLP_csc')[sel_ev] == 1

        llp_r[m] += list(np.abs(T.array('gLLP_decay_vertex_r'))[sel_ev][sel_llp])
        llp_z[m] += list(np.abs(T.array('gLLP_decay_vertex_z'))[sel_ev][sel_llp])
	#print("-----> ", len(llp_r[m]), len(llp_z[m]))
	ones = np.ones(len(T.array('higgsPtWeight'))) #(T.array('higgsPtWeight')*T.array('weight')*T.array('pileupWeight')*T.array('metSF')
        weight_llp[m] += list(ones[sel_ev]*lumi[i])
        cluster_index = '3'

        ########### SELECTION: CLUSTERS ############
        sel_rechitcluster = T.array('cscRechitCluster'+cluster_index+'_match_gLLP_csc')
        sel_ev = np.logical_and(sel_ev,sel_rechitcluster.sum()==1)
        r_denom[m] += list(np.abs(T.array('cscRechitCluster' + cluster_index + '_match_gLLP_decay_r'))[sel_rechitcluster][sel_ev].flatten())
        z_denom[m] += list(np.abs(T.array('cscRechitCluster' + cluster_index + '_match_gLLP_decay_z'))[sel_rechitcluster][sel_ev].flatten())
	# print(T.array('cscRechitCluster' + cluster_index + '_match_gLLP_decay_r'), sel_rechitcluster.sum().sum(), sel_ev.sum())

        weight_denom[m] += list(ones[sel_ev]*lumi[i])
	print("start ", sum(sel_rechitcluster.sum()))
        #sel_rechitcluster = np.logical_and(sel_rechitcluster, np.logical_and(T.array('cscRechitCluster'+cluster_index+'Time') > -5.0, T.array('cscRechitCluster'+cluster_index+'Time') < 12.5))
	#sel_rechitcluster = np.logical_and(sel_rechitcluster, np.abs(T.array('cscRechitCluster'+cluster_index+'Eta')) < 2.1)
        #sel_rechitcluster = np.logical_and(sel_rechitcluster, T.array('cscRechitCluster'+cluster_index+'NRechitChamberPlus11') <= 0)
        #sel_rechitcluster = np.logical_and(sel_rechitcluster, T.array('cscRechitCluster'+cluster_index+'NRechitChamberPlus12') <= 0)
        #sel_rechitcluster = np.logical_and(sel_rechitcluster, T.array('cscRechitCluster'+cluster_index+'NRechitChamberMinus11') <= 0)
        #sel_rechitcluster = np.logical_and(sel_rechitcluster, T.array('cscRechitCluster'+cluster_index+'NRechitChamberMinus12') <= 0)

        #sel_rechitcluster = np.logical_and( sel_rechitcluster, T.array('cscRechitCluster'+cluster_index+'JetVetoPt') < 10)
        #sel_rechitcluster = np.logical_and( sel_rechitcluster, T.array('cscRechitCluster'+cluster_index+'MuonVetoPt') < 20)
        #sel_rechitcluster = np.logical_and(sel_rechitcluster, T.array('cscRechitCluster'+cluster_index+'TimeSpread') < 20)
        #sel_rechitcluster = np.logical_and(sel_rechitcluster, T.array('cscRechitCluster'+cluster_index+'_match_MB1_0p4') <= 0)
        #sel_rechitcluster = np.logical_and(sel_rechitcluster, T.array('cscRechitCluster'+cluster_index+'_match_RB1_0p4') <= 0)
        #sel_rechitcluster = np.logical_and(sel_rechitcluster, T.array('cscRechitCluster'+cluster_index+'_match_RE12_0p4') <= 0)
        print("end ", sum(sel_rechitcluster.sum()))
        print()
        sel_ev = np.logical_and(sel_ev,sel_rechitcluster.sum()==1)

        ##### bdt variables ####

        cscRechitClusterAvgStation = T.array('cscRechitCluster' + cluster_index + 'AvgStation5')[sel_rechitcluster][sel_ev].flatten()
        cscRechitClusterNStation = T.array('cscRechitCluster' + cluster_index + 'NStation5')[sel_rechitcluster][sel_ev].flatten()
        cscRechitClusterEta = T.array('cscRechitCluster' + cluster_index + 'Eta')[sel_rechitcluster][sel_ev].flatten()
        cscRechitClusterPhi = T.array('cscRechitCluster' + cluster_index + 'Phi')[sel_rechitcluster][sel_ev].flatten()
        cscRechitClusterAvgStation10 = T.array('cscRechitCluster' + cluster_index + 'AvgStation10')[sel_rechitcluster][sel_ev].flatten()
        cscRechitClusterNStation10 = T.array('cscRechitCluster' + cluster_index + 'NStation10')[sel_rechitcluster][sel_ev].flatten()

        if cut_based_version == 'v4':
            cond2 = np.logical_and(np.abs(cscRechitClusterAvgStation10)==2, np.abs(cscRechitClusterEta) < 1.6)
            cond3 = np.logical_and(np.abs(cscRechitClusterAvgStation10)==3, np.abs(cscRechitClusterEta) < 1.6)
            cond4 = np.logical_and(np.abs(cscRechitClusterAvgStation10)==4, np.abs(cscRechitClusterEta) < 1.8)
            cond1 = np.logical_and(cscRechitClusterNStation10==1, np.logical_or(np.logical_or(np.abs(cscRechitClusterAvgStation10)==1, cond2), np.logical_or(cond3, cond4)))
            cond2 = np.logical_and(cscRechitClusterNStation10 > 1, np.abs(cscRechitClusterEta) < 1.9)
            bdt_sel = np.logical_or(np.logical_or(cond1, cond2), np.logical_or(cond3, cond4))
            # print("bdt_sel-> ", bdt_sel)
        # sel_rechitcluster = np.logical_and(sel_rechitcluster, bdt_score >= BDT_CUT)
        # sel_ev = np.logical_and(sel_ev,sel_rechitcluster.sum()>=1)
        else:
            print("ERROR")

        r_nom[m] += list(np.abs(T.array('cscRechitCluster' + cluster_index + '_match_gLLP_decay_r'))[sel_rechitcluster][sel_ev].flatten()[bdt_sel])
        z_nom[m] += list(np.abs(T.array('cscRechitCluster' + cluster_index + '_match_gLLP_decay_z'))[sel_rechitcluster][sel_ev].flatten()[bdt_sel])
        weight_nom[m] += list(ones[sel_ev][bdt_sel]*lumi[i])
        nRechitscut = T.array('cscRechitCluster' + cluster_index + 'Size')[sel_rechitcluster][sel_ev].flatten()[bdt_sel]
        r_nRechits[m] += list(np.abs(T.array('cscRechitCluster' + cluster_index + '_match_gLLP_decay_r'))[sel_rechitcluster][sel_ev].flatten()[bdt_sel][nRechitscut>=130])
        z_nRechits[m] += list(np.abs(T.array('cscRechitCluster' + cluster_index + '_match_gLLP_decay_z'))[sel_rechitcluster][sel_ev].flatten()[bdt_sel][nRechitscut>=130])
        weight_nRechits[m] += list(ones[sel_ev][bdt_sel][nRechitscut>=130]*lumi[i])
    #print(sum(weight_llp[m]))

#print(weight_nom[m])
#print(weight_denom[m])
#print('higgsPtWeight', 'weight', 'pileupWeight', 'metSF')
#print(sum(T.array('higgsPtWeight')), sum(T.array('weight')), sum(T.array('pileupWeight')), sum(T.array('metSF')))

# importlib.reload(sys.modules['CMS_lumi'])

names = ['cluster_efficiency', 'sig_eff', 'nRechit_efficiency']
m = '40'
# print(weight_denom)
for name in names:
    for r in [0,1]:
        if r == 1:continue
        if not name == 'cluster_efficiency':continue
        c = rt.TCanvas('c','c', 800, 600)
        if r:
            bins = [20,0,700]
            xaxis_title = 'LLP decay R [cm]'
        else:
            bins = [35, 550, 1075]
            bins = [35, 400, 1075]
            xaxis_title = 'z decay position [cm]'

            #hm = create_TH1D(cscRechitCluster_match_gLLP_r[eta.flatten()==True].flatten(), 'hm1', axis_title = [xaxis_title,'Signal Efficiency'], binning=bins)
        #hb = create_TH1D(cscRechitCluster_match_gLLP_r.flatten(), 'hb1', axis_title = [xaxis_title,'Signal Efficiency'], binning=bins)
        if r:
            if name == 'cluster_efficiency':
                hm = create_TH1D(r_denom[m], 'hb1', axis_title = [xaxis_title,'Cluster Efficiency'], binning=bins, weights = weight_denom[m])
                hb = create_TH1D(llp_r[m], 'hb1', axis_title = [xaxis_title,'Cluster Efficiency'], binning=bins, weights = weight_llp[m])
            elif name == 'sig_eff':
                hm = create_TH1D(r_nom[m], 'hm1', axis_title = [xaxis_title,'Signal Efficiency'], binning=bins, weights = weight_nom[m])
                hb = create_TH1D(llp_r[m], 'hb1', axis_title = [xaxis_title,'Signal Efficiency'], binning=bins, weights = weight_llp[m])
            else:
                hm = create_TH1D(r_nRechits[m], 'hm1', axis_title = [xaxis_title,'N_{rechits} Cut Efficiency'], binning=bins, weights = weight_nRechits[m])
                hb = create_TH1D(r_nom[m], 'hb1', axis_title = [xaxis_title,'N_{rechits} Cut Efficiency'], binning=bins, weights = weight_nom[m])
        else:
        #     if name == 'cluster_efficiency':hm = create_TH1D(z_denom, 'hm1', axis_title = [xaxis_title,'Signal Efficiency'], binning=bins)
        #     else: hm = create_TH1D(z_nom, 'hb1', axis_title = [xaxis_title,'Signal Efficiency'], binning=bins)
        #     hb = create_TH1D(llp_z, 'hb1', axis_title = [xaxis_title,'Signal Efficiency'], binning=bins)
            if name == 'cluster_efficiency':
                hm = create_TH1D(z_denom[m], 'hb1', axis_title = [xaxis_title,'Cluster Efficiency'], binning=bins, weights = weight_denom[m])
                hb = create_TH1D(llp_z[m], 'hb1', axis_title = [xaxis_title,'Cluster Efficiency'], binning=bins, weights = weight_llp[m])
            elif name == 'sig_eff':
                hm = create_TH1D(z_nom[m], 'hm1', axis_title = [xaxis_title,'Signal Efficiency'], binning=bins, weights = weight_nom[m])
                hb = create_TH1D(llp_z[m], 'hb1', axis_title = [xaxis_title,'Signal Efficiency'], binning=bins, weights = weight_llp[m])
            else:
                hm = create_TH1D(z_nRechits[m], 'hm1', axis_title = [xaxis_title,'N_{rechits} Cut Efficiency'], binning=bins, weights = weight_nRechits[m])
                hb = create_TH1D(z_nom[m], 'hb1', axis_title = [xaxis_title,'N_{rechits} Cut Efficiency'], binning=bins, weights = weight_nom[m])

        #print(len(z_denom), len(llp_z))
        #print(len(r_denom), len(llp_r))

        pEff1 = rt.TEfficiency(hm,hb)
        pEff1.SetLineColor(std_color_list[2])
        pEff1.SetLineWidth(3)
        pEff1.Draw()
	print(hm.GetEntries(), "    ", hb.GetEntries())
#        tdrstyle.setTDRStyle()
        
        c.Draw()
        pEff1.GetPaintedGraph().GetHistogram().GetYaxis().SetTitleOffset(0.0);

        ymax = pEff1.GetPaintedGraph().GetHistogram().GetMaximum()
        if name == 'nRechit_efficiency' and r == 1:ymax*=1.05
        ymin = pEff1.GetPaintedGraph().GetHistogram().GetMinimum()
        xmin = pEff1.GetPaintedGraph().GetHistogram().GetXaxis().GetXmin()
        xmax = pEff1.GetPaintedGraph().GetHistogram().GetXaxis().GetXmax()
        pEff1.GetPaintedGraph().GetHistogram().SetMaximum(ymax)
        # print(xmax, xmin, ymin)
        if r:
            l = rt.TLine(350,ymin, 350, ymax)
            l.SetLineWidth(2)
            l.SetLineStyle(2)
            l.Draw()
            text = rt.TLatex()
            text.SetTextSize(0.05)
            text.DrawLatex(max(50,xmin+25), ymax*0.92, "Inner ring")
            text.DrawLatex(400, ymax*0.92, "Outer ring")
        else:
            boxes = []

            boxes.append(rt.TBox(xmin,ymin,568,ymax)) #in front of ME11
            boxes.append(rt.TBox(632,ymin,671,ymax)) #between ME11 and ME12
            boxes.append(rt.TBox(724,ymin,789,ymax)) #between ME12 and station2
            boxes.append(rt.TBox(849,ymin,911,ymax)) #between station2 and station3
            boxes.append(rt.TBox(970,ymin,1002,ymax)) #between station3 and station4
            boxes.append(rt.TBox(1073,ymin,xmax,ymax)) #beyond CMS
            for b in boxes:
                b.SetFillColor(15)
                b.SetFillStyle(3001)
                b.Draw('same')
            
            l = rt.TLatex()
            l.SetTextSize(0.08)
            l.SetTextColor(12)
            l.SetTextAngle(90)
#             l.DrawLatex(xmin+80, ymax*0.6, "Steel")

            l.DrawLatex(xmin+80, ymax*0.55, "Steel")

            l2 = rt.TLatex()
            l2.SetTextSize(0.06)
            l2.SetTextColor(13)
            l2.SetTextAngle(90)
            l2.DrawLatex(1110, ymax*0.5, "Beyond CMS")
            text = rt.TLatex()
            text.SetTextSize(0.04)
            text.DrawLatex(570, ymax*1.01, "ME1/1")
            text.DrawLatex(660, ymax*1.01, "ME1/2-3")
            text.DrawLatex(795, ymax*1.01, "ME2")
            text.DrawLatex(920, ymax*1.01, "ME3")
            text.DrawLatex(1015, ymax*1.01, "ME4")

        CMS_lumi.cmsText = "CMS"
        CMS_lumi.iPos = 11
        CMS_lumi.relPosX = 0.1
        CMS_lumi.relPosY = 0.05
        CMS_lumi.writeExtraText = 1
        CMS_lumi.extraText   = "Simulation Preliminary"
#        CMS_lumi.extraText   = "Simulation Supplementary"

        CMS_lumi.CMS_lumi(c, 0, 11)
        pEff1.Draw('same')

        if cut_based: bdt_name = "cut_based_"+cut_based_version
        # print(bdt_name)
        # outDir = "../../plots/MuonSystem_Analysis/"+bdt_name+analyzer_version
        outDir = "./"
        if not os.path.isdir(outDir):os.makedirs(outDir)
        fileName = outDir+name+"_"
        if r: fileName += "r_denomAllLLP_m"+m
        else:fileName += "z_denomAllLLP_m"+m
        if CMS_lumi.writeExtraText and CMS_lumi.extraText == "Simulation Preliminary": fileName += "_pas"
        elif CMS_lumi.writeExtraText and CMS_lumi.extraText == "Simulation Supplementary": fileName += "_paper"

        c.SaveAs(fileName + ".png")
#        c.SaveAs(fileName + ".pdf")
#        c.SaveAs(fileName + ".C")
print("")
import numpy.ma as ma
# importlib.reload(sys.modules['CMS_lumi'])

names = ['sig_eff']
for m in ['15', '40', '55']:
    if not m == '15':continue
    c = rt.TCanvas('c','c', 800, 600)
    
    # bins = [10,100,700, 20, 550, 1075]
    # title = ['Decay R [cm]', 'Decay Z [cm]']
    
    bins = [20, 400, 1100, 20,100,750]
    title = ['r decay position [cm]', '|z| decay position [cm]']
    #print("------------------>    ", sum(z_nom[m]), sum(r_nom[m]), sum(llp_z[m]), sum(llp_r[m]))
    hm = create_TH2D(np.column_stack((z_nom[m],r_nom[m])), 'hm1', axis_title = [title[1], title[0], 'Signal efficiency'], binning=bins, weights = weight_nom[m])
    hb = create_TH2D(np.column_stack((llp_z[m],llp_r[m])), 'hb1', axis_title = [title[1], title[0], 'Signal efficiency'], binning=bins, weights = weight_llp[m])

    #print("")
    #print("2D plot ", hm.Integral(), hb.Integral())
    #print("")

    #print(np.column_stack((z_nom[m],r_nom[m])).shape, np.column_stack((z_nom[m],r_nom[m])).sum())
    #print(len(z_denom) , len(llp_z))
    #print(len(r_denom) , len(llp_r))

    pEff1 = rt.TEfficiency(hm,hb)
    hist = hm.Clone()
    
    for y in range(1, hm.GetYaxis().GetNbins()+1):
        for x in range(1, hm.GetXaxis().GetNbins()+1):
            global_bin = pEff1.GetGlobalBin(x,y)
            if pEff1.GetEfficiency(global_bin)<0.01 and pEff1.GetEfficiency(global_bin)>0:
                print(hm.GetXaxis().GetBinCenter(x), hm.GetYaxis().GetBinCenter(y), pEff1.GetEfficiency(global_bin))
                pEff1.SetPassedEvents(global_bin, 0)
            if x==hm.GetXaxis().GetNbins():pEff1.SetPassedEvents(global_bin, 0)
            hist.SetBinContent(x,y,pEff1.GetEfficiency(global_bin))
    
#     pEff1.Draw('colz')
    hist.SetMaximum(1.0)
    hist.Draw('colz')

    c.SetTopMargin(0.07)

    c.SetRightMargin(0.2)
    c.Draw()
    
    CMS_lumi.cmsText = "CMS"
    iPos = 0
    CMS_lumi.writeExtraText = 1
    CMS_lumi.extraText = "Simulation"
#     CMS_lumi.extraText = "Simulation Preliminary"


    if( iPos==0 ): CMS_lumi.relPosX = 0.15
    # CMS_lumi.CMS_lumi(c, 4, 0)
    CMS_lumi.CMS_lumi(c, 0, iPos)

    boxes = []
    boxes.append(rt.TBox(400,402,661,449)) #MB1
    boxes.append(rt.TBox(400,490,661,533)) #MB2
    boxes.append(rt.TBox(400,597,661,636)) #MB3
    boxes.append(rt.TBox(400,700,661,738)) #MB4

    boxes.append(rt.TBox(400,295,661,380)) #solenoid

    boxes.append(rt.TBox(791,357,850,700)) #ME2/2
    boxes.append(rt.TBox(911,357,970,700)) #ME3/2
    boxes.append(rt.TBox(1002,357,1063,700)) #ME4/2

    boxes.append(rt.TBox(789,139,850,345)) #ME2/1
    boxes.append(rt.TBox(915,160,970,345)) #ME3/1
    boxes.append(rt.TBox(1002,178,1063,345)) #ME4/1

    boxes.append(rt.TBox(580,100,632,275)) #ME1/1
    boxes.append(rt.TBox(668,275,724,465)) #ME1/2
    boxes.append(rt.TBox(686,505,724,700)) #ME1/3
    boxes.append(create_TGraph([400.5,560,560,440,400.5,400.5],[100.5,100.5,270,270,190,100.5])) #hcal

    for b in boxes:
        b.SetFillColorAlpha(15,0)
        b.SetLineWidth(2)
        b.SetLineColor(1)
    #     b.SetFillStyle(3001)
        b.Draw('l same')

    steel = []
    steel.append(rt.TBox(400,449,661,490)) #MB1
    steel.append(rt.TBox(400,533,661,597)) #MB2
    steel.append(rt.TBox(400,636,661,700)) #MB3
    x=[633,633, 1100, 1100, 1075, 1075, 1063, 1063,1002,1002,975,975,912,912,849,849,789,789,724,724, 633]
    y=[275, 100, 100, 150, 150, 700,700,150,150,700,700,140,140,700,700,110,110,700,700,275,275]
    # print(len(x),len(y))
    steel.append(create_TGraph(x,y)) 

    for b in steel:
        b.SetFillColorAlpha(15,0.5)
        b.SetFillStyle(3001)
#         b.SetFillColorAlpha(15,0.9)
        b.Draw('fsame')

    text = rt.TLatex()
    text.SetTextSize(0.04)

    text.DrawLatex(665, 705, "ME1/3")
    text.DrawLatex(780, 705, "ME2/2")
    text.DrawLatex(900, 705, "ME3/2")
    text.DrawLatex(990, 705, "ME4/2")

    text.DrawLatex(410, 415, "MB1")
    text.DrawLatex(410, 500, "MB2")
    text.DrawLatex(410, 605, "MB3")
    text.DrawLatex(410, 708, "MB4")

    text.DrawLatex(680, 140, "Steel")
    text.DrawLatex(410, 660, "Steel")
    text.DrawLatex(410, 330, "Solenoid")
    text.DrawLatex(430, 150, "HCAL")

    text.SetTextAngle(90)
    text.DrawLatex(615, 110, "ME1/1")
    text.DrawLatex(695, 331, "ME1/2")
    text.DrawLatex(830, 145, "ME2/1")
    text.DrawLatex(955, 170, "ME3/1")
    text.DrawLatex(1040, 190, "ME4/1")
    
    if cut_based: bdt_name = "cut_based_"+cut_based_version
    # print(bdt_name)
    #outDir = "/storage/af/user/christiw/gpu/christiw/llp/delayed_jet_analyzer//plots/MuonSystem_Analysis/"+bdt_name+analyzer_version
    outDir = "./"
    if not os.path.isdir(outDir):os.makedirs(outDir)
    if CMS_lumi.extraText == 'Simulation Preliminary':
        #c.SaveAs(outDir+"2D_sig_eff_"+decay+"_m"+m+"_pas.C")
        #c.SaveAs(outDir+"2D_sig_eff_"+decay+"_m"+m+"_pas.pdf")
        c.SaveAs(outDir+"2D_sig_eff_"+decay+"_m"+m+"_pas.png")
    else:
        #c.SaveAs(outDir+"2D_sig_eff_"+decay+"_m"+m+".C")
        c.SaveAs(outDir+"2D_sig_eff_"+decay+"_m"+m+".png")
        #c.SaveAs(outDir+"2D_sig_eff_"+decay+"_m"+m+".pdf")
    