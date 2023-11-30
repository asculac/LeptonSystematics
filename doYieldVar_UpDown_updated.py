import math
import ROOT
ROOT.gROOT.SetBatch(True)
from ROOT import *

gRandom.SetSeed(101)

year = "2016"
#year = "2017"
year = "2018"

if (year=="2016"):
   lumi = 36.92
   #muSF=TFile.Open("/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/RunIILegacy/200205_CutBased/MuonSF/final_HZZ_SF_2016_legacy_mupogsysts.root")
   # muSF=TFile.Open("$CMSSW_BASE/src/ZZAnalysis/AnalysisStep/data/LeptonEffScaleFactors/ScaleFactors_mu_Moriond2017_v2.root")
   muSF=TFile.Open("/afs/cern.ch/user/t/tjavaid/public/UL_SFrms/final_HZZ_SF_2016UL_mupogsysts_newLoose_abseta3_rmsSystLooseSIPIso.root") #latest muons with "RMS"
elif(year=="2017"):
   lumi = 41.53
   #muSF=TFile.Open("/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/RunIILegacy/200205_CutBased/MuonSF/final_HZZ_SF_2017_rereco_mupogsysts_3010.root")
   # muSF=TFile.Open("$CMSSW_BASE/src/ZZAnalysis/AnalysisStep/data/LeptonEffScaleFactors/ScaleFactors_mu_Moriond2018_final.root")
   muSF=TFile.Open("/afs/cern.ch/user/t/tjavaid/public/UL_SFrms/final_HZZ_2017UL_SF_mupogSysts_newLoose_abseta3_rmsSystLooseSIPIso.root") #latest muons with "RMS"
   eleSF_ID_nogap=TFile.Open("/afs/cern.ch/work/a/asculac/REGULAR_CRZLL_TEST/CMSSW_10_6_26/src/ZZAnalysis/AnalysisStep/data/LeptonEffScaleFactors/ElectronSF_UL2017_nogap.root")
   eleSF_ID_gap=TFile.Open("/afs/cern.ch/work/a/asculac/REGULAR_CRZLL_TEST/CMSSW_10_6_26/src/ZZAnalysis/AnalysisStep/data/LeptonEffScaleFactors/ElectronSF_UL2017_gap.root")
   eleSF_RECO_high=TFile.Open("/afs/cern.ch/work/a/asculac/REGULAR_CRZLL_TEST/CMSSW_10_6_26/src/ZZAnalysis/AnalysisStep/data/LeptonEffScaleFactors/egammaEffi_ptAbove20.txt_EGM2D_UL_2017_RMS.root")
   eleSF_RECO_low=TFile.Open("/afs/cern.ch/work/a/asculac/REGULAR_CRZLL_TEST/CMSSW_10_6_26/src/ZZAnalysis/AnalysisStep/data/LeptonEffScaleFactors/egammaEffi_ptBelow20.txt_EGM2D_UL_2017_RMS.root")
elif (year=="2018"):
   lumi = 59.74
   #muSF=TFile.Open("/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/RunIILegacy/200205_CutBased/MuonSF/final_HZZ_SF_2018_rereco_mupogsysts_3010.root")
   # muSF=TFile.Open("$CMSSW_BASE/src/ZZAnalysis/AnalysisStep/data/LeptonEffScaleFactors/final_HZZ_muon_SF_2018RunA2D_ER_2702.root")   
   muSF=TFile.Open("/afs/cern.ch/user/t/tjavaid/public/UL_SFrms/final_HZZ_2017UL_SF_mupogSysts_newLoose_abseta3_rmsSystLooseSIPIso.root") #latest muons with "RMS"

# folder = "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/Run2UL_22_apr23/MC2017/ggH_Chunks/" # Define input folder name
# folder = "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/Run2UL_22_apr23/MC2018/" # Define input folder name
folder = "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/Run2UL_22_apr23/MC2016/" # Define input folder name

file_name = '/ZZ4lAnalysis.root' # Define input file name


correlated_leptons = True
uncorrelated_leptons = True
corr_factor = 0

# List of samples to run on
List = [
'ggH125',
#'ZZTo4l',
# 'ggTo4l'
]

def sigma_event(rho, SF1, SF2, SF3, SF4, sigma1, sigma2, sigma3, sigma4):
   rez = (sigma1/SF1)**2 + (sigma2/SF2)**2 + (sigma3/SF3)**2 + (sigma4/SF4)**2 + 2*rho*( sigma1*sigma2/SF1/SF2 + sigma1*sigma3/SF1/SF3 + sigma1*sigma4/SF1/SF4 + sigma2*sigma3/SF2/SF3 + sigma2*sigma4/SF2/SF4 + sigma3*sigma4/SF3/SF4)
   if rez < 0.000001:
      return 0
   else:
      return rez

for Type in List:
   # Open the root file and get tree
   mc=TFile.Open(folder+Type+file_name)
   tree = mc.Get("ZZTree/candTree")
   counters = mc.Get("ZZTree/Counters")
   NGen = counters.GetBinContent(40)
   
   MuonSFHisto=muSF.Get("FINAL")
   MuonUncHisto=muSF.Get("ERROR")

   # EleIDSFHisto=eleSF_ID_nogap.Get("EGamma_SF2D")
   # EleIDUncHisto=eleSF_ID_nogap.Get("EGamma_SF2D")

   # EleIDSFHisto_gap=eleSF_ID_gap.Get("EGamma_SF2D")
   # EleIDUncHisto_gap=eleSF_ID_gap.Get("EGamma_SF2D")

   # EleRECO_low_SFHisto=eleSF_RECO_low.Get("EGamma_SF2D")
   # EleRECO_low_UncHisto=eleSF_RECO_low.Get("EGamma_SF2D")
   
   # EleRECO_high_SFHisto=eleSF_RECO_high.Get("EGamma_SF2D")
   # EleRECO_high_UncHisto=eleSF_RECO_high.Get("EGamma_SF2D")
   # Set mass range and resolution
   lo=105
   hi=140
   nbins=(hi-lo)*1 # 1 GeV resolution


   # Calculate nominal value of yield
   h_nom = TH1F("h_nom", "h_nom", nbins, lo, hi)
   
   tree.Draw("ZZMass >> h_nom" , "(abs(LepLepId[0]) == 11 && abs(LepLepId[3]) == 11)*overallEventWeight*1000*xsec*"+str(lumi)+"/"+str(NGen))
   nom_yield_4e = h_nom.Integral()
   
   tree.Draw("ZZMass >> h_nom" , "(abs(LepLepId[0]) == 13 && abs(LepLepId[3]) == 13)*overallEventWeight*1000*xsec*"+str(lumi)+"/"+str(NGen))
   nom_yield_4mu = h_nom.Integral()
   
   tree.Draw("ZZMass >> h_nom" , "(abs(abs(LepLepId[0]) - (LepLepId[3])) == 2)      *overallEventWeight*1000*xsec*"+str(lumi)+"/"+str(NGen))
   nom_yield_2e2mu = h_nom.Integral()

   # print "Processing {} ...".format(Type)
   br_data = 0
   br_fs4e = 0
   br_fs4mu = 0
   br_fs2e2mu = 0
   
   # Set yields for up and dn variations to zero
   yield_4e_up = []
   yield_4mu_up = []
   yield_2e2mu_e_up = []
   yield_2e2mu_mu_up = []
   yield_4e_dn = []
   yield_4mu_dn = []
   yield_2e2mu_e_dn = []
   yield_2e2mu_mu_dn = []
   
   uncor_4e_up = []
   uncor_4e_dn = []
   uncor_4mu_up = []
   uncor_4mu_dn = []
   uncor_2e2mu_e_up = []
   uncor_2e2mu_e_dn = []
   uncor_2e2mu_mu_up = []
   uncor_2e2mu_mu_dn = []
   
   # Uncrrelate effect of trigger/reco/selection
   variation_name = []
   variation_name.append("TRIGGER")
   variation_name.append("RECO_SEL")
   for k in range (0,2):
      yield_4e_up.append(0.)
      yield_4mu_up.append(0.)
      yield_2e2mu_e_up.append(0.)
      yield_2e2mu_mu_up.append(0.)
      yield_4e_dn.append(0.)
      yield_4mu_dn.append(0.)
      yield_2e2mu_e_dn.append(0.)
      yield_2e2mu_mu_dn.append(0.)
   
      uncor_4e_up.append(0.)
      uncor_4e_dn.append(0.)
      uncor_4mu_up.append(0.)
      uncor_4mu_dn.append(0.)
      uncor_2e2mu_e_up.append(0.)
      uncor_2e2mu_e_dn.append(0.)
      uncor_2e2mu_mu_up.append(0.)
      uncor_2e2mu_mu_dn.append(0.)


   for event in tree:# Loop over all events in tree
      br_data+=1
      
      mass4l = tree.ZZMass
      if( mass4l < 105. or mass4l > 140.): continue # Skip events that are not in the mass window
      
      idL1 = abs(tree.LepLepId[0])
      idL3 = abs(tree.LepLepId[3])

      if (idL1==11 and idL3==11):
         br_fs4e +=1
      elif (idL1==13 and idL3==13):
         br_fs4mu += 1
      elif (abs(idL1-idL3)==2):
         br_fs2e2mu += 1

      # if(br_data % int((tree.GetEntries()/10)) == 0):
      #    print "{} %".format(str(100*br_data/tree.GetEntries() + 1))
      

    
      GENmassZZ = tree.GenHMass
      
      # Calculate nominal weigh using central value of SF
      weight_nom = event.overallEventWeight * 1000 * lumi * event.xsec * event.L1prefiringWeight / NGen
      # Nominal value of total SF, product of 4 lepton nominal SF
      SF_tot_nom = event.dataMCWeight
      
      SF_lep_trig = []
      err_lep_trig_up = []
      err_lep_trig_dn = []

      SF_lep = []
      err_lep_up = []
      err_lep_dn = []
      
      for i in range (0,4):
         # print(err_lep_trig_dn)
         SF_lep_trig.append(0.)
         # err_lep_trig_up.append(0.)
         # err_lep_trig_dn.append(0.)
     
         SF_lep.append(0.)
         err_lep_up.append(0.)
         err_lep_dn.append(0.)
         
      for i in range (0,4):
         # Hard-code trigger SF and unc
         err_lep_trig_dn = [0.01,0,0,0]
         SF_lep_trig[i] = 1.
 
         if(idL1==11 and idL3==11 and event.LepPt[3] < 12):
               # err_lep_trig_up[i] = 0.02
               err_lep_trig_up = [0.02,0,0,0]
               # err_lep_trig_dn[i] = 0.11
               err_lep_trig_dn = [0.11,0,0,0]
               # print(err_lep_trig_dn)
         if(idL1==11 and idL3==11 and event.LepPt[3] >= 12):
               # err_lep_trig_up[i] = 0.005
               err_lep_trig_up = [0.005,0,0,0]
               # err_lep_trig_dn[i] = 0.01
               err_lep_trig_dn = [0.01,0,0,0]
               # print(err_lep_trig_dn)
         if(idL1==13 and idL3==13 and event.LepPt[3] < 7):
               # err_lep_trig_up[i] = 0.001
               err_lep_trig_up = [0.001,0,0,0]
               # err_lep_trig_dn[i] = 0.032
               err_lep_trig_dn = [0.032,0,0,0]
               # print(err_lep_trig_dn)
         if(idL1==13 and idL3==13 and event.LepPt[3] >= 7 and event.LepPt[3] < 12):
               # err_lep_trig_up[i] = 0.001
               err_lep_trig_up = [0.001,0,0,0]
               # err_lep_trig_dn[i] = 0.015
               err_lep_trig_dn = [0.015,0,0,0]
               # print(err_lep_trig_dn)
         if(idL1==13 and idL3==13 and event.LepPt[3] >= 12):
               # err_lep_trig_up[i] = 0.001
               err_lep_trig_up = [0.001,0,0,0]
               # err_lep_trig_dn[i] = 0.015
               err_lep_trig_dn = [0.015,0,0,0]
               # print(err_lep_trig_dn)
                  
         if(abs(idL1-idL3)==2 and event.LepPt[3] < 7):
               # err_lep_trig_up[i] = 0.005
               err_lep_trig_up = [0.005,0,0,0]
               # err_lep_trig_dn[i] = 0.08
               err_lep_trig_dn = [0.08,0,0,0]
               # print(err_lep_trig_dn)
         if(abs(idL1-idL3)==2 and event.LepPt[3] >= 7 and event.LepPt[3] < 12):
               # err_lep_trig_up[i] = 0.005
               err_lep_trig_up = [0.005,0,0,0]
               # err_lep_trig_dn[i] = 0.032
               err_lep_trig_dn = [0.032,0,0,0]
               # print(err_lep_trig_dn)
         if(abs(idL1-idL3)==2 and event.LepPt[3] >= 12):
               # err_lep_trig_up[i] = 0.001
               err_lep_trig_up = [0.001,0,0,0]
               # err_lep_trig_dn[i] = 0.010
               err_lep_trig_dn = [0.01,0,0,0]
       

         # Load reco and selection SF and unc
         #for electrons, reconstrucion SF unc and HZZ selection SF unc are combined and provided in the branch stored in our NTuple;
         # print(err_lep_trig_up)
         if (abs(tree.LepLepId[i]) == 11):

            SF_lep[i] = event.LepSF[i]
            err_lep_up[i] = event.LepSF_Unc[i]
            err_lep_dn[i] = event.LepSF_Unc[i]

            # for electrons, RMS RECO AND ID SF directly retrieved from SF file
            # isGap=False
            # isGap = (abs(tree.LepEta[i])>1.443 and abs(tree.LepEta[i])<1.567) 
            # SF_lep[i] = event.LepSF[i]

            # if(tree.LepPt[i]>=20 and not(isGap)):
            #    SF_ID =  EleIDSFHisto.GetBinContent( EleIDSFHisto.GetXaxis().FindBin(tree.LepEta[i]), EleIDSFHisto.GetYaxis().FindBin(min(tree.LepPt[i],199.99))) 
            #    SF_RECO = EleRECO_high_SFHisto.GetBinContent( EleRECO_high_SFHisto.GetXaxis().FindBin(tree.LepEta[i]), EleRECO_high_SFHisto.GetYaxis().FindBin(min(tree.LepPt[i],199.99))) 
               
            #    SF_ID_UNC= EleIDUncHisto.GetBinError(EleIDUncHisto.GetXaxis().FindBin(tree.LepEta[i]),EleIDUncHisto.GetYaxis().FindBin(min(tree.LepPt[i],199.99)))
            #    SF_RECO_UNC = EleRECO_high_UncHisto.GetBinError(EleRECO_high_UncHisto.GetXaxis().FindBin(tree.LepEta[i]),EleRECO_high_UncHisto.GetYaxis().FindBin(min(tree.LepPt[i],199.99)))

            #    err_lep_up[i] = math.sqrt( SF_RECO_UNC*SF_RECO_UNC/(SF_RECO*SF_RECO) + SF_ID_UNC*SF_ID_UNC/(SF_ID*SF_ID))
            #    err_lep_dn[i] = math.sqrt( SF_RECO_UNC*SF_RECO_UNC/(SF_RECO*SF_RECO) + SF_ID_UNC*SF_ID_UNC/(SF_ID*SF_ID))
            #    # print(tree.LepPt[i],tree.LepEta[i])
            #    # print(SF_ID,SF_RECO)
            #    # print(SF_ID_UNC,SF_RECO_UNC)
            #    # print(err_lep_up[i])
            # elif(tree.LepPt[i]<20 and not(isGap)):
            #    SF_ID =  EleIDSFHisto.GetBinContent( EleIDSFHisto.GetXaxis().FindBin(tree.LepEta[i]), EleIDSFHisto.GetYaxis().FindBin(min(tree.LepPt[i],199.99))) 
            #    SF_RECO = EleRECO_low_SFHisto.GetBinContent( EleRECO_low_SFHisto.GetXaxis().FindBin(tree.LepEta[i]), EleRECO_low_SFHisto.GetYaxis().FindBin(15.)) 
            #    # print(tree.LepPt[i],tree.LepEta[i])
            #    # print(SF_ID,SF_RECO)
            #    SF_ID_UNC= EleIDUncHisto.GetBinError(EleIDUncHisto.GetXaxis().FindBin(tree.LepEta[i]),EleIDUncHisto.GetYaxis().FindBin(min(tree.LepPt[i],199.99)))
            #    SF_RECO_UNC = EleRECO_low_UncHisto.GetBinError(EleRECO_low_UncHisto.GetXaxis().FindBin(tree.LepEta[i]),EleRECO_low_UncHisto.GetYaxis().FindBin(15.))
            #    err_lep_up[i] = math.sqrt( SF_RECO_UNC*SF_RECO_UNC/(SF_RECO*SF_RECO) + SF_ID_UNC*SF_ID_UNC/(SF_ID*SF_ID))
            #    err_lep_dn[i] = math.sqrt( SF_RECO_UNC*SF_RECO_UNC/(SF_RECO*SF_RECO) + SF_ID_UNC*SF_ID_UNC/(SF_ID*SF_ID))
            #    # print(err_lep_up[i])
            # elif(tree.LepPt[i]<20 and isGap):
            #    SF_ID =  EleIDSFHisto_gap.GetBinContent( EleIDSFHisto_gap.GetXaxis().FindBin(tree.LepEta[i]), EleIDSFHisto_gap.GetYaxis().FindBin(min(tree.LepPt[i],199.99))) 
            #    SF_RECO = EleRECO_low_SFHisto.GetBinContent( EleRECO_low_SFHisto.GetXaxis().FindBin(tree.LepEta[i]), EleRECO_low_SFHisto.GetYaxis().FindBin(15.)) 
               
            #    SF_ID_UNC= EleIDUncHisto_gap.GetBinError(EleIDUncHisto_gap.GetXaxis().FindBin(tree.LepEta[i]),EleIDUncHisto_gap.GetYaxis().FindBin(min(tree.LepPt[i],199.99)))
            #    SF_RECO_UNC = EleRECO_low_UncHisto.GetBinError(EleRECO_low_UncHisto.GetXaxis().FindBin(tree.LepEta[i]),EleRECO_low_UncHisto.GetYaxis().FindBin(15.))
            #    err_lep_up[i] = math.sqrt( SF_RECO_UNC*SF_RECO_UNC/(SF_RECO*SF_RECO) + SF_ID_UNC*SF_ID_UNC/(SF_ID*SF_ID))
            #    err_lep_dn[i] = math.sqrt( SF_RECO_UNC*SF_RECO_UNC/(SF_RECO*SF_RECO) + SF_ID_UNC*SF_ID_UNC/(SF_ID*SF_ID))
            #    # print(err_lep_up[i])
            # elif(tree.LepPt[i]>=20 and isGap):
            #    SF_ID =  EleIDSFHisto_gap.GetBinContent( EleIDSFHisto_gap.GetXaxis().FindBin(tree.LepEta[i]), EleIDSFHisto_gap.GetYaxis().FindBin(min(tree.LepPt[i],199.99))) 
            #    SF_RECO = EleRECO_high_SFHisto.GetBinContent( EleRECO_high_SFHisto.GetXaxis().FindBin(tree.LepEta[i]), EleRECO_high_SFHisto.GetYaxis().FindBin(min(tree.LepPt[i],199.99))) 
               
            #    SF_ID_UNC= EleIDUncHisto_gap.GetBinError(EleIDUncHisto_gap.GetXaxis().FindBin(tree.LepEta[i]),EleIDUncHisto_gap.GetYaxis().FindBin(min(tree.LepPt[i],199.99)))
            #    SF_RECO_UNC = EleRECO_high_UncHisto.GetBinError(EleRECO_high_UncHisto.GetXaxis().FindBin(tree.LepEta[i]),EleRECO_high_UncHisto.GetYaxis().FindBin(min(tree.LepPt[i],199.99)))
            #    err_lep_up[i] = math.sqrt( SF_RECO_UNC*SF_RECO_UNC/(SF_RECO*SF_RECO) + SF_ID_UNC*SF_ID_UNC/(SF_ID*SF_ID))
            #    err_lep_dn[i] = math.sqrt( SF_RECO_UNC*SF_RECO_UNC/(SF_RECO*SF_RECO) + SF_ID_UNC*SF_ID_UNC/(SF_ID*SF_ID))
            #    # print(err_lep_up[i])
            
            #print "ELE SF: ", SF_lep, " ########## event.LepPt[i]", event.LepPt[i], " ### event.LepEta[i]", event.LepEta[i]
            #print "   ELE: err_lep_up", err_lep_up
            #print "   ELE: err_lep_dn", err_lep_dn
         elif (abs(tree.LepLepId[i]) == 13):
            # SF_lep[i] = event.LepSF[i]
            # err_lep_up[i] = event.LepSF_Unc[i]
            # err_lep_dn[i] = event.LepSF_Unc[i]

            #print "MUON SF: ", SF_lep, " ########## event.LepPt[i]"
            #print "   MUON: err_lep_up", err_lep_up
            #print "   MUON: err_lep_dn", err_lep_dn
            
            # for muons, directly retrieved from muon SF file
            SF_lep[i] = MuonSFHisto.GetBinContent( MuonSFHisto.GetXaxis().FindBin(tree.LepEta[i]), MuonSFHisto.GetYaxis().FindBin(min(tree.LepPt[i],199.99)))
            err_lep_up[i] = MuonUncHisto.GetBinContent(MuonUncHisto.GetXaxis().FindBin(tree.LepEta[i]),MuonUncHisto.GetYaxis().FindBin(min(tree.LepPt[i],199.99)))
            err_lep_dn[i] = MuonUncHisto.GetBinContent(MuonUncHisto.GetXaxis().FindBin(tree.LepEta[i]),MuonUncHisto.GetYaxis().FindBin(min(tree.LepPt[i],199.99)))

            #print "MUON SF test: ", SF_lep_test, " ########## event.Lep[i]"
            #print "   MUON test: err_lep_up", err_lep_up_test
            #print "   MUON test: err_lep_dn", err_lep_dn_test
         
      for k in range (0,2):
         # Calculate each variation independently
         if(k==0):
            TRIG = 1
            RECO_SEL = 0
         elif(k==1):
            TRIG = 0
            RECO_SEL = 1


         SF_var_up = 1.
         SF_var_dn = 1.
         SF_var_e_up = 1.
         SF_var_e_dn = 1.
         SF_var_mu_up = 1.
         SF_var_mu_dn = 1.
         
         if (idL1==11 and idL3==11):
            # Vary SF of each lepton up and down
            for i in range (0,4):
               if ( correlated_leptons ):
                  SF_var_up *= (SF_lep_trig[i] + TRIG*err_lep_trig_up[i]) * (SF_lep[i] + RECO_SEL*err_lep_up[i])
                  SF_var_dn *= (SF_lep_trig[i] - TRIG*err_lep_trig_dn[i]) * (SF_lep[i] - RECO_SEL*err_lep_dn[i])
           
            if (uncorrelated_leptons):
               uncor_4e_up[k] += TRIG*(sigma_event(corr_factor,SF_lep_trig[0],SF_lep_trig[1],SF_lep_trig[2],SF_lep_trig[3],err_lep_trig_up[0],err_lep_trig_up[1],err_lep_trig_up[2],err_lep_trig_up[3])) + RECO_SEL*(sigma_event(corr_factor,SF_lep[0],SF_lep[1],SF_lep[2],SF_lep[3],err_lep_up[0],err_lep_up[1],err_lep_up[2],err_lep_up[3]))
               uncor_4e_dn[k] += TRIG*(sigma_event(corr_factor,SF_lep_trig[0],SF_lep_trig[1],SF_lep_trig[2],SF_lep_trig[3],err_lep_trig_dn[0],err_lep_trig_dn[1],err_lep_trig_dn[2],err_lep_trig_dn[3])) + RECO_SEL*(sigma_event(corr_factor,SF_lep[0],SF_lep[1],SF_lep[2],SF_lep[3],err_lep_dn[0],err_lep_dn[1],err_lep_dn[2],err_lep_dn[3])) 
               
            yield_4e_up[k] += weight_nom/SF_tot_nom * SF_var_up
            yield_4e_dn[k] += weight_nom/SF_tot_nom * SF_var_dn
      

         elif (idL1==13 and idL3==13):
            for i in range (0,4):
               if ( correlated_leptons ):
                  SF_var_up *= (SF_lep_trig[i] + TRIG*err_lep_trig_up[i]) * (SF_lep[i] + RECO_SEL*err_lep_up[i])
                  SF_var_dn *= (SF_lep_trig[i] - TRIG*err_lep_trig_dn[i]) * (SF_lep[i] - RECO_SEL*err_lep_dn[i])
            if (uncorrelated_leptons):
               uncor_4mu_up[k] += TRIG*(sigma_event(corr_factor,SF_lep_trig[0],SF_lep_trig[1],SF_lep_trig[2],SF_lep_trig[3],err_lep_trig_up[0],err_lep_trig_up[1],err_lep_trig_up[2],err_lep_trig_up[3])) + RECO_SEL*(sigma_event(corr_factor,SF_lep[0],SF_lep[1],SF_lep[2],SF_lep[3],err_lep_up[0],err_lep_up[1],err_lep_up[2],err_lep_up[3]))
               uncor_4mu_dn[k] += TRIG*(sigma_event(corr_factor,SF_lep_trig[0],SF_lep_trig[1],SF_lep_trig[2],SF_lep_trig[3],err_lep_trig_dn[0],err_lep_trig_dn[1],err_lep_trig_dn[2],err_lep_trig_dn[3])) + RECO_SEL*(sigma_event(corr_factor,SF_lep[0],SF_lep[1],SF_lep[2],SF_lep[3],err_lep_dn[0],err_lep_dn[1],err_lep_dn[2],err_lep_dn[3])) 

            yield_4mu_up[k] += weight_nom/SF_tot_nom * SF_var_up
            yield_4mu_dn[k] += weight_nom/SF_tot_nom * SF_var_dn

         elif (abs(idL1-idL3)==2):
            # Vary electron SF while fixing muon and vice-versa
            if ( idL1 == 11):
               if ( correlated_leptons ):
                  SF_var_e_up  = (SF_lep_trig[0] + TRIG*err_lep_trig_up[0]) * (SF_lep[0] + RECO_SEL*err_lep_up[0]) * (SF_lep_trig[1] + TRIG*err_lep_trig_up[1]) * (SF_lep[1] + RECO_SEL*err_lep_up[1]) * SF_lep_trig[2] * SF_lep[2] * SF_lep_trig[3] * SF_lep[3]
                  SF_var_mu_up = (SF_lep_trig[2] + TRIG*err_lep_trig_up[2]) * (SF_lep[2] + RECO_SEL*err_lep_up[2]) * (SF_lep_trig[3] + TRIG*err_lep_trig_up[3]) * (SF_lep[3] + RECO_SEL*err_lep_up[3]) * SF_lep_trig[0] * SF_lep[0] * SF_lep_trig[1] * SF_lep[1]

                  SF_var_e_dn  = (SF_lep_trig[0] + TRIG*err_lep_trig_dn[0]) * (SF_lep[0] + RECO_SEL*err_lep_dn[0]) * (SF_lep_trig[1] + TRIG*err_lep_trig_dn[1]) * (SF_lep[1] + RECO_SEL*err_lep_dn[1]) * SF_lep_trig[2] * SF_lep[2] * SF_lep_trig[3] * SF_lep[3]
                  SF_var_mu_dn = (SF_lep_trig[2] + TRIG*err_lep_trig_dn[2]) * (SF_lep[2] + RECO_SEL*err_lep_dn[2]) * (SF_lep_trig[3] + TRIG*err_lep_trig_dn[3]) * (SF_lep[3] + RECO_SEL*err_lep_dn[3]) * SF_lep_trig[0] * SF_lep[0] * SF_lep_trig[1] * SF_lep[1]
               
               if (uncorrelated_leptons):
                  uncor_2e2mu_e_up[k] += TRIG*(sigma_event(corr_factor,SF_lep_trig[0],SF_lep_trig[1],SF_lep_trig[2],SF_lep_trig[3],err_lep_trig_up[0],err_lep_trig_up[1],0,0)) + RECO_SEL*(sigma_event(corr_factor,SF_lep[0],SF_lep[1],SF_lep[2],SF_lep[3],err_lep_up[0],err_lep_up[1],0,0))
                  uncor_2e2mu_e_dn[k] += TRIG*(sigma_event(corr_factor,SF_lep_trig[0],SF_lep_trig[1],SF_lep_trig[2],SF_lep_trig[3],err_lep_trig_dn[0],err_lep_trig_dn[1],0,0)) + RECO_SEL*(sigma_event(corr_factor,SF_lep[0],SF_lep[1],SF_lep[2],SF_lep[3],err_lep_dn[0],err_lep_dn[1],0,0))
                  uncor_2e2mu_mu_up[k] += TRIG*(sigma_event(corr_factor,SF_lep_trig[0],SF_lep_trig[1],SF_lep_trig[2],SF_lep_trig[3],0,0,err_lep_trig_up[2],err_lep_trig_up[3])) + RECO_SEL*(sigma_event(corr_factor,SF_lep[0],SF_lep[1],SF_lep[2],SF_lep[3],0,0,err_lep_up[2],err_lep_dn[3]))
                  uncor_2e2mu_mu_dn[k] += TRIG*(sigma_event(corr_factor,SF_lep_trig[0],SF_lep_trig[1],SF_lep_trig[2],SF_lep_trig[3],0,0,err_lep_trig_dn[2],err_lep_trig_dn[3])) + RECO_SEL*(sigma_event(corr_factor,SF_lep[0],SF_lep[1],SF_lep[2],SF_lep[3],0,0,err_lep_dn[2],err_lep_dn[3]))
            
            elif ( idL1 == 13):
               if ( correlated_leptons ):
                  SF_var_mu_up  = (SF_lep_trig[0] + TRIG*err_lep_trig_up[0]) * (SF_lep[0] + RECO_SEL*err_lep_up[0]) * (SF_lep_trig[1] + TRIG*err_lep_trig_up[1]) * (SF_lep[1] + RECO_SEL*err_lep_up[1]) * SF_lep_trig[2] * SF_lep[2] * SF_lep_trig[3] * SF_lep[3]
                  SF_var_e_up = (SF_lep_trig[2] + TRIG*err_lep_trig_up[2]) * (SF_lep[2] + RECO_SEL*err_lep_up[2]) * (SF_lep_trig[3] + TRIG*err_lep_trig_up[3]) * (SF_lep[3] + RECO_SEL*err_lep_up[3]) * SF_lep_trig[0] * SF_lep[0] * SF_lep_trig[1] * SF_lep[1]

                  SF_var_mu_dn  = (SF_lep_trig[0] + TRIG*err_lep_trig_dn[0]) * (SF_lep[0] + RECO_SEL*err_lep_dn[0]) * (SF_lep_trig[1] + TRIG*err_lep_trig_dn[1]) * (SF_lep[1] + RECO_SEL*err_lep_dn[1]) * SF_lep_trig[2] * SF_lep[2] * SF_lep_trig[3] * SF_lep[3]
                  SF_var_e_dn = (SF_lep_trig[2] + TRIG*err_lep_trig_dn[2]) * (SF_lep[2] + RECO_SEL*err_lep_dn[2]) * (SF_lep_trig[3] + TRIG*err_lep_trig_dn[3]) * (SF_lep[3] + RECO_SEL*err_lep_dn[3]) * SF_lep_trig[0] * SF_lep[0] * SF_lep_trig[1] * SF_lep[1]
                     
               if (uncorrelated_leptons):
                  uncor_2e2mu_mu_up[k] += TRIG*(sigma_event(corr_factor,SF_lep_trig[0],SF_lep_trig[1],SF_lep_trig[2],SF_lep_trig[3],err_lep_trig_up[0],err_lep_trig_up[1],0,0)) + RECO_SEL*(sigma_event(corr_factor,SF_lep[0],SF_lep[1],SF_lep[2],SF_lep[3],err_lep_up[0],err_lep_up[1],0,0))
                  uncor_2e2mu_mu_dn[k] += TRIG*(sigma_event(corr_factor,SF_lep_trig[0],SF_lep_trig[1],SF_lep_trig[2],SF_lep_trig[3],err_lep_trig_dn[0],err_lep_trig_dn[1],0,0)) + RECO_SEL*(sigma_event(corr_factor,SF_lep[0],SF_lep[1],SF_lep[2],SF_lep[3],err_lep_dn[0],err_lep_dn[1],0,0))
                  uncor_2e2mu_e_up[k] += TRIG*(sigma_event(corr_factor,SF_lep_trig[0],SF_lep_trig[1],SF_lep_trig[2],SF_lep_trig[3],0,0,err_lep_trig_up[2],err_lep_trig_up[3])) + RECO_SEL*(sigma_event(corr_factor,SF_lep[0],SF_lep[1],SF_lep[2],SF_lep[3],0,0,err_lep_up[2],err_lep_dn[3]))
                  uncor_2e2mu_e_dn[k] += TRIG*(sigma_event(corr_factor,SF_lep_trig[0],SF_lep_trig[1],SF_lep_trig[2],SF_lep_trig[3],0,0,err_lep_trig_dn[2],err_lep_trig_dn[3])) + RECO_SEL*(sigma_event(corr_factor,SF_lep[0],SF_lep[1],SF_lep[2],SF_lep[3],0,0,err_lep_dn[2],err_lep_dn[3]))

            yield_2e2mu_e_up[k]  += weight_nom/SF_tot_nom * SF_var_e_up
            yield_2e2mu_mu_up[k] += weight_nom/SF_tot_nom * SF_var_mu_up
               
            yield_2e2mu_e_dn[k]  += weight_nom/SF_tot_nom * SF_var_e_dn
            yield_2e2mu_mu_dn[k] += weight_nom/SF_tot_nom * SF_var_mu_dn

   comb_4e_up = 0.
   comb_4mu_up = 0.
   comb_2e2mu_e_up = 0.
   comb_2e2mu_mu_up = 0.

   comb_4e_dn = 0.
   comb_4mu_dn = 0.
   comb_2e2mu_e_dn = 0.
   comb_2e2mu_mu_dn = 0.
   
   comb_uncor_4e_up = 0.
   comb_uncor_4mu_up = 0.
   comb_uncor_2e2mu_e_up = 0.
   comb_uncor_2e2mu_mu_up = 0.
   
   comb_uncor_4e_dn = 0.
   comb_uncor_4mu_dn = 0.
   comb_uncor_2e2mu_e_dn = 0.
   comb_uncor_2e2mu_mu_dn = 0.

   # Sum relative uncertainties in quadrature
   for k in range (0,2):
      if ( correlated_leptons ):
         comb_4e_up += ((yield_4e_up[k]-nom_yield_4e)/nom_yield_4e*100)**2
         comb_4mu_up += ((yield_4mu_up[k]-nom_yield_4mu)/nom_yield_4mu*100)**2
         comb_2e2mu_e_up += ((yield_2e2mu_e_up[k]-nom_yield_2e2mu)/nom_yield_2e2mu*100)**2
         comb_2e2mu_mu_up += ((yield_2e2mu_mu_up[k]-nom_yield_2e2mu)/nom_yield_2e2mu*100)**2
         
         comb_4e_dn += ((-yield_4e_dn[k]+nom_yield_4e)/nom_yield_4e*100)**2
         comb_4mu_dn += ((-yield_4mu_dn[k]+nom_yield_4mu)/nom_yield_4mu*100)**2
         comb_2e2mu_e_dn += ((-yield_2e2mu_e_dn[k]+nom_yield_2e2mu)/nom_yield_2e2mu*100)**2
         comb_2e2mu_mu_dn += ((-yield_2e2mu_mu_dn[k]+nom_yield_2e2mu)/nom_yield_2e2mu*100)**2
         
      if (uncorrelated_leptons):
         comb_uncor_4e_up += (uncor_4e_up[k]/br_fs4e)
         comb_uncor_4mu_up += (uncor_4mu_up[k]/br_fs4mu)
         comb_uncor_2e2mu_e_up += (uncor_2e2mu_e_up[k]/br_fs2e2mu)
         comb_uncor_2e2mu_mu_up += (uncor_2e2mu_mu_up[k]/br_fs2e2mu)
         
         comb_uncor_4e_dn += (uncor_4e_dn[k]/br_fs4e)
         comb_uncor_4mu_dn += (uncor_4mu_dn[k]/br_fs4mu)
         comb_uncor_2e2mu_e_dn += (uncor_2e2mu_e_dn[k]/br_fs2e2mu)
         comb_uncor_2e2mu_mu_dn += (uncor_2e2mu_mu_dn[k]/br_fs2e2mu)
      
      # Print the results
      if ( correlated_leptons ):
         print "=================================================================="
         print "| CORRELATED LEPTON UNCERTAINTIES FOR {} IN {} SAMPLE |".format(variation_name[k],Type)
         print "=================================================================="
         print "4e = +{:.1f} % -{:.1f}%  4mu = +{:.1f}% -{:.1f} %".format(((yield_4e_up[k]-nom_yield_4e)/nom_yield_4e*100),((-yield_4e_dn[k]+nom_yield_4e)/nom_yield_4e*100),((yield_4mu_up[k]-nom_yield_4mu)/nom_yield_4mu*100),((-yield_4mu_dn[k]+nom_yield_4mu)/nom_yield_4mu*100))
         print "2e2mu/ele = +{:.1f} % -{:.1f} % 2e2mu/mu = +{:.1f}% -{:.1f}%".format(((yield_2e2mu_e_up[k]-nom_yield_2e2mu)/nom_yield_2e2mu*100),((-yield_2e2mu_e_dn[k]+nom_yield_2e2mu)/nom_yield_2e2mu*100),((yield_2e2mu_mu_up[k]-nom_yield_2e2mu)/nom_yield_2e2mu*100),((-yield_2e2mu_mu_dn[k]+nom_yield_2e2mu)/nom_yield_2e2mu*100))
   
      if (uncorrelated_leptons):
         print "=================================================================="
         print "| UNCORRELATED LEPTON UNCERTAINTIES FOR {} IN {} SAMPLE |".format(variation_name[k],Type)
         print "=================================================================="
         print "4e = +{:.1f}% -{:.1f}%  4mu = +{:.1f}% -{:.1f}%".format(100*(uncor_4e_up[k]/br_fs4e)**0.5,100*(uncor_4e_dn[k]/br_fs4e)**0.5,100*(uncor_4mu_up[k]/br_fs4mu)**0.5,100*(uncor_4mu_dn[k]/br_fs4mu)**0.5)
         print "2e2mu/ele = +{:.1f}% -{:.1f}% 2e2mu/mu = +{:.1f}% -{:.1f}%".format(100*(uncor_2e2mu_e_up[k]/br_fs2e2mu)**0.5,100*(uncor_2e2mu_e_dn[k]/br_fs2e2mu)**0.5,100*(uncor_2e2mu_mu_up[k]/br_fs2e2mu)**0.5,100*(uncor_2e2mu_mu_dn[k]/br_fs2e2mu)**0.5)

   if ( correlated_leptons ):
      print "==============================="
      print "| CORRELATED COMBINATION |"
      print "==============================="
      print "4e = +{:.1f}% -{:.1f}%  4mu = +{:.1f}% -{:.1f}%".format(comb_4e_up**0.5,comb_4e_dn**0.5,comb_4mu_up**0.5,comb_4mu_dn**0.5)
      print "2e2mu/ele = +{:.1f}% -{:.1f}% 2e2mu/mu = +{:.1f}% -{:.1f}%".format(comb_2e2mu_e_up**0.5,comb_2e2mu_e_dn**0.5,comb_2e2mu_mu_up**0.5,comb_2e2mu_mu_dn**0.5)
      print "=============================================================================="

   if (uncorrelated_leptons):
      print "================================"
      print "| UNCORRELATED UNCORRELATED |"
      print "================================"
      print "4e = +{:.1f}% -{:.1f}%  4mu = +{:.1f}% -{:.1f}%".format(100*(comb_uncor_4e_up**0.5),100*(comb_uncor_4e_dn**0.5),100*(comb_uncor_4mu_up**0.5),100*(comb_uncor_4mu_dn**0.5))
      print "2e2mu/ele = +{:.1f}% -{:.1f}% 2e2mu/mu = +{:.1f}% -{:.1f}%".format(100*(comb_uncor_2e2mu_e_up**0.5),100*(comb_uncor_2e2mu_e_dn**0.5),100*(comb_uncor_2e2mu_mu_up**0.5),100*(comb_uncor_2e2mu_mu_dn**0.5))
      print "=============================================================================="



