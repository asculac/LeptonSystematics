import ROOT
ROOT.gROOT.SetBatch(True)
from ROOT import *

gRandom.SetSeed(101)


lumi = 59.7 # Set lumi to be used for MC scaling

folder = '/eos/user/t/tsculac/BigStuff/Run2/2018/' # Define input folder name
file_name = '/ZZ4lAnalysis.root' # Define input dile name

correlated_leptons = True
uncorrelated_leptons = True
corr_factor = 0

# List of samples to run on
List = [
'ggH125',
#'ZZTo4l',
#'ggTo4l'
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
   
   muSF=TFile.Open("RunBCDEF_SF_ID_syst.root")
   MuonUncHisto=muSF.Get("NUM_LooseID_DEN_genTracks_pt_abseta_syst")

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

   print "Processing {} ...".format(Type)
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
   variation_name.append("RECO")
   variation_name.append("SELECTION")
   for k in range (0,3):
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

      if(br_data % int((tree.GetEntries()/10)) == 0):
         print "{} %".format(str(100*br_data/tree.GetEntries() + 1))
      

    
      GENmassZZ = tree.GenHMass
      
      # Calculate nominal weigh using central value of SF
      weight_nom = event.overallEventWeight*1000*lumi*event.xsec/NGen
      # Nominal value of total SF, product of 4 lepton nominal SF
      SF_tot_nom = event.dataMCWeight
      
      SF_lep_trig = []
      err_lep_trig_up = []
      err_lep_trig_dn = []
      SF_lep_reco = []
      err_lep_reco_up = []
      err_lep_reco_dn = []
      SF_lep_sel = []
      err_lep_sel_up = []
      err_lep_sel_dn = []
      
      for i in range (0,4):
         SF_lep_trig.append(0.)
         err_lep_trig_up.append(0.)
         err_lep_trig_dn.append(0.)
         SF_lep_reco.append(0.)
         err_lep_reco_up.append(0.)
         err_lep_reco_dn.append(0.)
         SF_lep_sel.append(0.)
         err_lep_sel_up.append(0.)
         err_lep_sel_dn.append(0.)
         
      for i in range (0,4):
         # Hard-code trigger SF and unc
         SF_lep_trig[i] = 1.
         if(i == 0):
            if(idL1==11 and idL3==11 and event.LepPt[3] < 12):
               err_lep_trig_up[i] = 0.02
               err_lep_trig_dn[i] = 0.11
            if(idL1==11 and idL3==11 and event.LepPt[3] >= 12):
               err_lep_trig_up[i] = 0.005
               err_lep_trig_dn[i] = 0.01

            if(idL1==13 and idL3==13 and event.LepPt[3] < 7):
               err_lep_trig_up[i] = 0.001
               err_lep_trig_dn[i] = 0.032
            if(idL1==13 and idL3==13 and event.LepPt[3] >= 7 and event.LepPt[3] < 12):
               err_lep_trig_up[i] = 0.001
               err_lep_trig_dn[i] = 0.015
            if(idL1==13 and idL3==13 and event.LepPt[3] >= 12):
               err_lep_trig_up[i] = 0.001
               err_lep_trig_dn[i] = 0.015
                  
            if(abs(idL1-idL3)==2 and event.LepPt[3] < 7):
               err_lep_trig_up[i] = 0.005
               err_lep_trig_dn[i] = 0.08
            if(abs(idL1-idL3)==2 and event.LepPt[3] >= 7 and event.LepPt[3] < 12):
               err_lep_trig_up[i] = 0.005
               err_lep_trig_dn[i] = 0.032
            if(abs(idL1-idL3)==2 and event.LepPt[3] >= 12):
               err_lep_trig_up[i] = 0.001
               err_lep_trig_dn[i] = 0.010
         else:
            err_lep_trig_up[i] = 0.
            err_lep_trig_dn[i] = 0.
         # Load reco and selection SF and unc
         SF_lep_reco[i] = event.LepRecoSF[i]
         err_lep_reco_up[i] = event.LepRecoSF_Unc[i]
         err_lep_reco_dn[i] = event.LepRecoSF_Unc[i]
         if (abs(tree.LepLepId[i]) == 11):
            SF_lep_sel[i] = event.LepSelSF[i]
            err_lep_sel_up[i] = event.LepSelSF_Unc[i]
            err_lep_sel_dn[i] = event.LepSelSF_Unc[i]
         elif (abs(tree.LepLepId[i]) == 13):
            SF_lep_sel[i] = event.LepSelSF[i]
            err_lep_sel_up[i] = (event.LepSelSF_Unc[i]**2+MuonUncHisto.GetBinError(MuonUncHisto.GetXaxis().FindBin(min(tree.LepPt[i],199.99)),MuonUncHisto.GetYaxis().FindBin(abs(tree.LepEta[i])))**2)**0.5
            err_lep_sel_dn[i] = (event.LepSelSF_Unc[i]**2+MuonUncHisto.GetBinError(MuonUncHisto.GetXaxis().FindBin(min(tree.LepPt[i],199.99)),MuonUncHisto.GetYaxis().FindBin(abs(tree.LepEta[i])))**2)**0.5
               #print err_lep_sel_up[i]


      for k in range (0,3):
         # Calculate each variation independently
         if(k==0):
            TRIG = 1
            RECO = 0
            SEL  = 0
         elif(k==1):
            TRIG = 0
            RECO = 1
            SEL  = 0
         elif(k==2):
            TRIG = 0
            RECO = 0
            SEL  = 1

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
                  SF_var_up *= (SF_lep_trig[i] + TRIG*err_lep_trig_up[i]) * (SF_lep_reco[i] + RECO*err_lep_reco_up[i]) * (SF_lep_sel[i] + SEL*err_lep_sel_up[i])
                  SF_var_dn *= (SF_lep_trig[i] - TRIG*err_lep_trig_dn[i]) * (SF_lep_reco[i] - RECO*err_lep_reco_dn[i]) * (SF_lep_sel[i] - SEL*err_lep_sel_dn[i])
            if (uncorrelated_leptons):
               uncor_4e_up[k] +=  TRIG*(sigma_event(corr_factor,SF_lep_trig[0],SF_lep_trig[1],SF_lep_trig[2],SF_lep_trig[3],err_lep_trig_up[0],err_lep_trig_up[1],err_lep_trig_up[2],err_lep_trig_up[3]))+RECO*(sigma_event(corr_factor,SF_lep_reco[0],SF_lep_reco[1],SF_lep_reco[2],SF_lep_reco[3],err_lep_reco_up[0],err_lep_reco_up[1],err_lep_reco_up[2],err_lep_reco_up[3]))+SEL*(sigma_event(corr_factor,SF_lep_sel[0],SF_lep_sel[1],SF_lep_sel[2],SF_lep_sel[3],err_lep_sel_up[0],err_lep_sel_up[1],err_lep_sel_up[2],err_lep_sel_up[3]))
               uncor_4e_dn[k] += TRIG*(sigma_event(corr_factor,SF_lep_trig[0],SF_lep_trig[1],SF_lep_trig[2],SF_lep_trig[3],err_lep_trig_dn[0],err_lep_trig_dn[1],err_lep_trig_dn[2],err_lep_trig_dn[3]))+RECO*(sigma_event(corr_factor,SF_lep_reco[0],SF_lep_reco[1],SF_lep_reco[2],SF_lep_reco[3],err_lep_reco_dn[0],err_lep_reco_dn[1],err_lep_reco_dn[2],err_lep_reco_dn[3]))+SEL*(sigma_event(corr_factor,SF_lep_sel[0],SF_lep_sel[1],SF_lep_sel[2],SF_lep_sel[3],err_lep_sel_dn[0],err_lep_sel_dn[1],err_lep_sel_dn[2],err_lep_sel_dn[3]))
               
            yield_4e_up[k] += weight_nom/SF_tot_nom * SF_var_up
            yield_4e_dn[k] += weight_nom/SF_tot_nom * SF_var_dn
      

         elif (idL1==13 and idL3==13):
            for i in range (0,4):
               if ( correlated_leptons ):
                  SF_var_up *= (SF_lep_trig[i] + TRIG*err_lep_trig_up[i]) * (SF_lep_reco[i] + RECO*err_lep_reco_up[i]) * (SF_lep_sel[i] + SEL*err_lep_sel_up[i])
                  SF_var_dn *= (SF_lep_trig[i] - TRIG*err_lep_trig_dn[i]) * (SF_lep_reco[i] - RECO*err_lep_reco_dn[i]) * (SF_lep_sel[i] - SEL*err_lep_sel_dn[i])
            if (uncorrelated_leptons):
               uncor_4mu_up[k] += TRIG*(sigma_event(corr_factor,SF_lep_trig[0],SF_lep_trig[1],SF_lep_trig[2],SF_lep_trig[3],err_lep_trig_up[0],err_lep_trig_up[1],err_lep_trig_up[2],err_lep_trig_up[3]))+RECO*(sigma_event(corr_factor,SF_lep_reco[0],SF_lep_reco[1],SF_lep_reco[2],SF_lep_reco[3],err_lep_reco_up[0],err_lep_reco_up[1],err_lep_reco_up[2],err_lep_reco_up[3]))+SEL*(sigma_event(corr_factor,SF_lep_sel[0],SF_lep_sel[1],SF_lep_sel[2],SF_lep_sel[3],err_lep_sel_up[0],err_lep_sel_up[1],err_lep_sel_up[2],err_lep_sel_up[3]))
               uncor_4mu_dn[k] += TRIG*(sigma_event(corr_factor,SF_lep_trig[0],SF_lep_trig[1],SF_lep_trig[2],SF_lep_trig[3],err_lep_trig_dn[0],err_lep_trig_dn[1],err_lep_trig_dn[2],err_lep_trig_dn[3]))+RECO*(sigma_event(corr_factor,SF_lep_reco[0],SF_lep_reco[1],SF_lep_reco[2],SF_lep_reco[3],err_lep_reco_dn[0],err_lep_reco_dn[1],err_lep_reco_dn[2],err_lep_reco_dn[3]))+SEL*(sigma_event(corr_factor,SF_lep_sel[0],SF_lep_sel[1],SF_lep_sel[2],SF_lep_sel[3],err_lep_sel_dn[0],err_lep_sel_dn[1],err_lep_sel_dn[2],err_lep_sel_dn[3]))

            yield_4mu_up[k] += weight_nom/SF_tot_nom * SF_var_up
            yield_4mu_dn[k] += weight_nom/SF_tot_nom * SF_var_dn

         elif (abs(idL1-idL3)==2):
            # Vary electron SF while fixing muon and vice-versa
            if ( idL1 == 11):
               if ( correlated_leptons ):
                  SF_var_e_up  = (SF_lep_trig[0] + TRIG*err_lep_trig_up[0]) * (SF_lep_reco[0] + RECO*err_lep_reco_up[0]) * (SF_lep_sel[0] + SEL*err_lep_sel_up[0]) * (SF_lep_trig[1] + TRIG*err_lep_trig_up[1]) * (SF_lep_reco[1] + RECO*err_lep_reco_up[1]) * (SF_lep_sel[1] + SEL*err_lep_sel_up[1]) * SF_lep_trig[2] * SF_lep_reco[2] * SF_lep_sel[2] * SF_lep_trig[3] * SF_lep_reco[3] * SF_lep_sel[3]
                  SF_var_mu_up = (SF_lep_trig[2] + TRIG*err_lep_trig_up[2]) * (SF_lep_reco[2] + RECO*err_lep_reco_up[2]) * (SF_lep_sel[2] + SEL*err_lep_sel_up[2]) * (SF_lep_trig[3] + TRIG*err_lep_trig_up[3]) * (SF_lep_reco[3] + RECO*err_lep_reco_up[3]) * (SF_lep_sel[3] + SEL*err_lep_sel_up[3]) * SF_lep_trig[0] * SF_lep_reco[0] * SF_lep_sel[0] * SF_lep_trig[1] * SF_lep_reco[1] * SF_lep_sel[1]
               
                  SF_var_e_dn  = (SF_lep_trig[0] - TRIG*err_lep_trig_dn[0]) * (SF_lep_reco[0] - RECO*err_lep_reco_dn[0]) * (SF_lep_sel[0] - SEL*err_lep_sel_dn[0]) * (SF_lep_trig[1] - TRIG*err_lep_trig_dn[1]) * (SF_lep_reco[1] - RECO*err_lep_reco_dn[1]) * (SF_lep_sel[1] - SEL*err_lep_sel_dn[1]) * SF_lep_trig[2] * SF_lep_reco[2] * SF_lep_sel[2] * SF_lep_trig[3] * SF_lep_reco[3] * SF_lep_sel[3]
                  SF_var_mu_dn = (SF_lep_trig[2] - TRIG*err_lep_trig_dn[2]) * (SF_lep_reco[2] - RECO*err_lep_reco_dn[2]) * (SF_lep_sel[2] - SEL*err_lep_sel_dn[2]) * (SF_lep_trig[3] - TRIG*err_lep_trig_dn[3]) * (SF_lep_reco[3] - RECO*err_lep_reco_dn[3]) * (SF_lep_sel[3] - SEL*err_lep_sel_dn[3]) * SF_lep_trig[0] * SF_lep_reco[0] * SF_lep_sel[0] * SF_lep_trig[1] * SF_lep_reco[1] * SF_lep_sel[1]
               
               if (uncorrelated_leptons):
                  uncor_2e2mu_e_up[k] += TRIG*(sigma_event(corr_factor,SF_lep_trig[0],SF_lep_trig[1],SF_lep_trig[2],SF_lep_trig[3],err_lep_trig_up[0],err_lep_trig_up[1],0,0))+RECO*(sigma_event(corr_factor,SF_lep_reco[0],SF_lep_reco[1],SF_lep_reco[2],SF_lep_reco[3],err_lep_reco_up[0],err_lep_reco_up[1],0,0))+SEL*(sigma_event(corr_factor,SF_lep_sel[0],SF_lep_sel[1],SF_lep_sel[2],SF_lep_sel[3],err_lep_sel_up[0],err_lep_sel_up[1],0,0))
                  uncor_2e2mu_e_dn[k] += TRIG*(sigma_event(corr_factor,SF_lep_trig[0],SF_lep_trig[1],SF_lep_trig[2],SF_lep_trig[3],err_lep_trig_dn[0],err_lep_trig_dn[1],0,0))+RECO*(sigma_event(corr_factor,SF_lep_reco[0],SF_lep_reco[1],SF_lep_reco[2],SF_lep_reco[3],err_lep_reco_dn[0],err_lep_reco_dn[1],0,0))+SEL*(sigma_event(corr_factor,SF_lep_sel[0],SF_lep_sel[1],SF_lep_sel[2],SF_lep_sel[3],err_lep_sel_dn[0],err_lep_sel_dn[1],0,0))
                  uncor_2e2mu_mu_up[k] += TRIG*(sigma_event(corr_factor,SF_lep_trig[0],SF_lep_trig[1],SF_lep_trig[2],SF_lep_trig[3],0,0,err_lep_trig_up[2],err_lep_trig_up[3]))+RECO*(sigma_event(corr_factor,SF_lep_reco[0],SF_lep_reco[1],SF_lep_reco[2],SF_lep_reco[3],0,0,err_lep_reco_up[2],err_lep_reco_dn[3]))+SEL*(sigma_event(corr_factor,SF_lep_sel[0],SF_lep_sel[1],SF_lep_sel[2],SF_lep_sel[3],0,0,err_lep_sel_up[2],err_lep_sel_up[3]))
                  uncor_2e2mu_mu_dn[k] += TRIG*(sigma_event(corr_factor,SF_lep_trig[0],SF_lep_trig[1],SF_lep_trig[2],SF_lep_trig[3],0,0,err_lep_trig_dn[2],err_lep_trig_dn[3]))+RECO*(sigma_event(corr_factor,SF_lep_reco[0],SF_lep_reco[1],SF_lep_reco[2],SF_lep_reco[3],0,0,err_lep_reco_dn[2],err_lep_reco_dn[3]))+SEL*(sigma_event(corr_factor,SF_lep_sel[0],SF_lep_sel[1],SF_lep_sel[2],SF_lep_sel[3],0,0,err_lep_sel_dn[2],err_lep_sel_dn[3]))
            
            elif ( idL1 == 13):
               if ( correlated_leptons ):
                  SF_var_mu_up = (SF_lep_trig[0] + TRIG*err_lep_trig_up[0]) * (SF_lep_reco[0] + RECO*err_lep_reco_up[0]) * (SF_lep_sel[0] + SEL*err_lep_sel_up[0]) * (SF_lep_trig[1] + TRIG*err_lep_trig_up[1]) * (SF_lep_reco[1] + RECO*err_lep_reco_up[1]) * (SF_lep_sel[1] + SEL*err_lep_sel_up[1]) * SF_lep_trig[2] * SF_lep_reco[2] * SF_lep_sel[2] * SF_lep_trig[3] * SF_lep_reco[3] * SF_lep_sel[3]
                  SF_var_e_up  = (SF_lep_trig[2] + TRIG*err_lep_trig_up[2]) * (SF_lep_reco[2] + RECO*err_lep_reco_up[2]) * (SF_lep_sel[2] + SEL*err_lep_sel_up[2]) * (SF_lep_trig[3] + TRIG*err_lep_trig_up[3]) * (SF_lep_reco[3] + RECO*err_lep_reco_up[3]) * (SF_lep_sel[3] + SEL*err_lep_sel_up[3]) * SF_lep_trig[0] * SF_lep_reco[0] * SF_lep_sel[0] * SF_lep_trig[1] * SF_lep_reco[1] * SF_lep_sel[1]
                  
                  SF_var_mu_dn = (SF_lep_trig[0] - TRIG*err_lep_trig_dn[0]) * (SF_lep_reco[0] - RECO*err_lep_reco_dn[0]) * (SF_lep_sel[0] - SEL*err_lep_sel_dn[0]) * (SF_lep_trig[1] - TRIG*err_lep_trig_dn[1]) * (SF_lep_reco[1] - RECO*err_lep_reco_dn[1]) * (SF_lep_sel[1] - SEL*err_lep_sel_dn[1]) * SF_lep_trig[2] * SF_lep_reco[2] * SF_lep_sel[2] * SF_lep_trig[3] * SF_lep_reco[3] * SF_lep_sel[3]
                  SF_var_e_dn  = (SF_lep_trig[2] - TRIG*err_lep_trig_dn[2]) * (SF_lep_reco[2] - RECO*err_lep_reco_dn[2]) * (SF_lep_sel[2] - SEL*err_lep_sel_dn[2]) * (SF_lep_trig[3] - TRIG*err_lep_trig_dn[3]) * (SF_lep_reco[3] - RECO*err_lep_reco_dn[3]) * (SF_lep_sel[3] - SEL*err_lep_sel_dn[3]) * SF_lep_trig[0] * SF_lep_reco[0] * SF_lep_sel[0] * SF_lep_trig[1] * SF_lep_reco[1] * SF_lep_sel[1]
                     
               if (uncorrelated_leptons):
                  uncor_2e2mu_mu_up[k] += TRIG*(sigma_event(corr_factor,SF_lep_trig[0],SF_lep_trig[1],SF_lep_trig[2],SF_lep_trig[3],err_lep_trig_up[0],err_lep_trig_up[1],0,0))+RECO*(sigma_event(corr_factor,SF_lep_reco[0],SF_lep_reco[1],SF_lep_reco[2],SF_lep_reco[3],err_lep_reco_up[0],err_lep_reco_up[1],0,0))+SEL*(sigma_event(corr_factor,SF_lep_sel[0],SF_lep_sel[1],SF_lep_sel[2],SF_lep_sel[3],err_lep_sel_up[0],err_lep_sel_up[1],0,0))
                  uncor_2e2mu_mu_dn[k] += TRIG*(sigma_event(corr_factor,SF_lep_trig[0],SF_lep_trig[1],SF_lep_trig[2],SF_lep_trig[3],err_lep_trig_dn[0],err_lep_trig_dn[1],0,0))+RECO*(sigma_event(corr_factor,SF_lep_reco[0],SF_lep_reco[1],SF_lep_reco[2],SF_lep_reco[3],err_lep_reco_dn[0],err_lep_reco_dn[1],0,0))+SEL*(sigma_event(corr_factor,SF_lep_sel[0],SF_lep_sel[1],SF_lep_sel[2],SF_lep_sel[3],err_lep_sel_dn[0],err_lep_sel_dn[1],0,0))
                  uncor_2e2mu_e_up[k] += TRIG*(sigma_event(corr_factor,SF_lep_trig[0],SF_lep_trig[1],SF_lep_trig[2],SF_lep_trig[3],0,0,err_lep_trig_up[2],err_lep_trig_up[3]))+RECO*(sigma_event(corr_factor,SF_lep_reco[0],SF_lep_reco[1],SF_lep_reco[2],SF_lep_reco[3],0,0,err_lep_reco_up[2],err_lep_reco_dn[3]))+SEL*(sigma_event(corr_factor,SF_lep_sel[0],SF_lep_sel[1],SF_lep_sel[2],SF_lep_sel[3],0,0,err_lep_sel_up[2],err_lep_sel_up[3]))
                  uncor_2e2mu_e_dn[k] += TRIG*(sigma_event(corr_factor,SF_lep_trig[0],SF_lep_trig[1],SF_lep_trig[2],SF_lep_trig[3],0,0,err_lep_trig_dn[2],err_lep_trig_dn[3]))+RECO*(sigma_event(corr_factor,SF_lep_reco[0],SF_lep_reco[1],SF_lep_reco[2],SF_lep_reco[3],0,0,err_lep_reco_dn[2],err_lep_reco_dn[3]))+SEL*(sigma_event(corr_factor,SF_lep_sel[0],SF_lep_sel[1],SF_lep_sel[2],SF_lep_sel[3],0,0,err_lep_sel_dn[2],err_lep_sel_dn[3]))

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
   for k in range (0,3):
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





