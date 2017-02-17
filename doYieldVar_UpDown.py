import ROOT
ROOT.gROOT.SetBatch(True)
from ROOT import *

gRandom.SetSeed(101)


lumi = 35.876 # Set lumi to be used for MC scaling

folder = './Moriond_2017_v2/' # Define input folder name
file_name = '/ZZ4lAnalysis.root' # Define input dile name

# List of samples to run on
List = [
'ggH125',
#'ZZTo4l',
#'ggTo4l'
]

for Type in List:
   # Open the root file and get tree
   mc=TFile.Open(folder+Type+file_name)
   tree = mc.Get("ZZTree/candTree")
   counters = mc.Get("ZZTree/Counters")
   NGen = counters.GetBinContent(40)

   # Set mass range and resolution
   lo=105
   hi=140
   nbins=(hi-lo)*1 # 1 GeV resolution


   # Calculate central value of yield so we know how to center variation histograms
   h_nom = TH1F("h_nom", "h_nom", nbins, lo, hi)
   
   tree.Draw("ZZMass >> h_nom" , "(abs(LepLepId[0]) == 11 && abs(LepLepId[3]) == 11)*overallEventWeight*1000*xsec*"+str(lumi)+"/"+str(NGen))
   nom_yield_4e = h_nom.Integral()
   
   tree.Draw("ZZMass >> h_nom" , "(abs(LepLepId[0]) == 13 && abs(LepLepId[3]) == 13)*overallEventWeight*1000*xsec*"+str(lumi)+"/"+str(NGen))
   nom_yield_4mu = h_nom.Integral()
   
   tree.Draw("ZZMass >> h_nom" , "(abs(abs(LepLepId[0]) - (LepLepId[3])) == 2)      *overallEventWeight*1000*xsec*"+str(lumi)+"/"+str(NGen))
   nom_yield_2e2mu = h_nom.Integral()
   

   print "Processing {} ...".format(Type)
   br_data = 0
   
   # Set yields for toys to zero
   yield_4e_up = 0.
   yield_4mu_up = 0.
   yield_2e2mu_e_up = 0.
   yield_2e2mu_mu_up = 0.
   
   yield_4e_dn = 0.
   yield_4mu_dn = 0.
   yield_2e2mu_e_dn = 0.
   yield_2e2mu_mu_dn = 0.
   

   for event in tree:# Loop over all events in tree
      br_data+=1
      if(br_data % int((tree.GetEntries()/10)) == 0):
         print "{} %".format(str(100*br_data/tree.GetEntries() + 1))
      
      mass4l = tree.ZZMass
      if( mass4l < 105. or mass4l > 140.): continue # Skip events that are not in the mass window
      idL1 = abs(tree.LepLepId[0])
      idL3 = abs(tree.LepLepId[3])
    
      GENmassZZ = tree.GenHMass
      
      # Calculate nominal weigh using central value of SF
      weight_nom = event.overallEventWeight*1000*lumi*event.xsec/NGen
      # Nominal value of total SF, product of 4 lepton nominal SF
      SF_tot_nom = event.dataMCWeight
      
      SF_lep_trig = []
      err_lep_trig = []
      SF_lep_reco = []
      err_lep_reco = []
      SF_lep_sel = []
      err_lep_sel = []
      
      for i in range (0,4):
         SF_lep_trig.append(0.)
         err_lep_trig.append(0.)
         SF_lep_reco.append(0.)
         err_lep_reco.append(0.)
         SF_lep_sel.append(0.)
         err_lep_sel.append(0.)
         
      for i in range (0,4):
         # Read lepton SF and unc from tree
         SF_lep_trig[i] = 1. #hard-code value for trigger SF at the moment
         err_lep_trig[i] = 0. #hard-code value for trigger unc at the moment
         SF_lep_reco[i] = event.LepRecoSF[i]
         err_lep_reco[i] = event.LepRecoSF_Unc[i]
         SF_lep_sel[i] = event.LepSelSF[i]
         err_lep_sel[i] = event.LepSelSF_Unc[i]
      
      SF_var_up = 1.
      SF_var_dn = 1.
      SF_var_e_up = 1.
      SF_var_e_dn = 1.
      SF_var_mu_up = 1.
      SF_var_mu_dn = 1.
      
      if (idL1==11 and idL3==11):
         # Vary SF of each lepton up and down
         for i in range (0,4):
            SF_var_up *= (SF_lep_trig[i] + err_lep_trig[i]) * (SF_lep_reco[i] + err_lep_reco[i]) * (SF_lep_sel[i] + err_lep_sel[i])
            SF_var_dn *= (SF_lep_trig[i] - err_lep_trig[i]) * (SF_lep_reco[i] - err_lep_reco[i]) * (SF_lep_sel[i] - err_lep_sel[i])
         yield_4e_up += weight_nom/SF_tot_nom * SF_var_up
         yield_4e_dn += weight_nom/SF_tot_nom * SF_var_dn

      elif (idL1==13 and idL3==13):
         for i in range (0,4):
            SF_var_up *= (SF_lep_trig[i] + err_lep_trig[i]) * (SF_lep_reco[i] + err_lep_reco[i]) * (SF_lep_sel[i] + err_lep_sel[i])
            SF_var_dn *= (SF_lep_trig[i] - err_lep_trig[i]) * (SF_lep_reco[i] - err_lep_reco[i]) * (SF_lep_sel[i] - err_lep_sel[i])
         yield_4mu_up += weight_nom/SF_tot_nom * SF_var_up
         yield_4mu_dn += weight_nom/SF_tot_nom * SF_var_dn

      elif (abs(idL1-idL3)==2):
         # Vary electron SF while fixing muon and vice-versa
         if ( idL1 == 11):
            SF_var_e_up  = (SF_lep_trig[0] + err_lep_trig[0]) * (SF_lep_reco[0] + err_lep_reco[0]) * (SF_lep_sel[0] + err_lep_sel[0]) * (SF_lep_trig[1] + err_lep_trig[1]) * (SF_lep_reco[1] + err_lep_reco[1]) * (SF_lep_sel[1] + err_lep_sel[1]) * SF_lep_trig[2] * SF_lep_reco[2] * SF_lep_sel[2] * SF_lep_trig[3] * SF_lep_reco[3] * SF_lep_sel[3]
            SF_var_mu_up = (SF_lep_trig[2] + err_lep_trig[2]) * (SF_lep_reco[2] + err_lep_reco[2]) * (SF_lep_sel[2] + err_lep_sel[2]) * (SF_lep_trig[3] + err_lep_trig[3]) * (SF_lep_reco[3] + err_lep_reco[3]) * (SF_lep_sel[3] + err_lep_sel[3]) * SF_lep_trig[0] * SF_lep_reco[0] * SF_lep_sel[0] * SF_lep_trig[1] * SF_lep_reco[1] * SF_lep_sel[1]
         
            SF_var_e_dn  = (SF_lep_trig[0] - err_lep_trig[0]) * (SF_lep_reco[0] - err_lep_reco[0]) * (SF_lep_sel[0] - err_lep_sel[0]) * (SF_lep_trig[1] - err_lep_trig[1]) * (SF_lep_reco[1] - err_lep_reco[1]) * (SF_lep_sel[1] - err_lep_sel[1]) * SF_lep_trig[2] * SF_lep_reco[2] * SF_lep_sel[2] * SF_lep_trig[3] * SF_lep_reco[3] * SF_lep_sel[3]
            SF_var_mu_dn = (SF_lep_trig[2] - err_lep_trig[2]) * (SF_lep_reco[2] - err_lep_reco[2]) * (SF_lep_sel[2] - err_lep_sel[2]) * (SF_lep_trig[3] - err_lep_trig[3]) * (SF_lep_reco[3] - err_lep_reco[3]) * (SF_lep_sel[3] - err_lep_sel[3]) * SF_lep_trig[0] * SF_lep_reco[0] * SF_lep_sel[0] * SF_lep_trig[1] * SF_lep_reco[1] * SF_lep_sel[1]
         
         elif ( idL1 == 13):
            SF_var_mu_up  = (SF_lep_trig[0] + err_lep_trig[0]) * (SF_lep_reco[0] + err_lep_reco[0]) * (SF_lep_sel[0] + err_lep_sel[0]) * (SF_lep_trig[1] + err_lep_trig[1]) * (SF_lep_reco[1] + err_lep_reco[1]) * (SF_lep_sel[1] + err_lep_sel[1]) * SF_lep_trig[2] * SF_lep_reco[2] * SF_lep_sel[2] * SF_lep_trig[3] * SF_lep_reco[3] * SF_lep_sel[3]
            SF_var_e_up = (SF_lep_trig[2] + err_lep_trig[2]) * (SF_lep_reco[2] + err_lep_reco[2]) * (SF_lep_sel[2] + err_lep_sel[2]) * (SF_lep_trig[3] + err_lep_trig[3]) * (SF_lep_reco[3] + err_lep_reco[3]) * (SF_lep_sel[3] + err_lep_sel[3]) * SF_lep_trig[0] * SF_lep_reco[0] * SF_lep_sel[0] * SF_lep_trig[1] * SF_lep_reco[1] * SF_lep_sel[1]
            
            SF_var_mu_dn  = (SF_lep_trig[0] - err_lep_trig[0]) * (SF_lep_reco[0] - err_lep_reco[0]) * (SF_lep_sel[0] - err_lep_sel[0]) * (SF_lep_trig[1] - err_lep_trig[1]) * (SF_lep_reco[1] - err_lep_reco[1]) * (SF_lep_sel[1] - err_lep_sel[1]) * SF_lep_trig[2] * SF_lep_reco[2] * SF_lep_sel[2] * SF_lep_trig[3] * SF_lep_reco[3] * SF_lep_sel[3]
            SF_var_e_dn = (SF_lep_trig[2] - err_lep_trig[2]) * (SF_lep_reco[2] - err_lep_reco[2]) * (SF_lep_sel[2] - err_lep_sel[2]) * (SF_lep_trig[3] - err_lep_trig[3]) * (SF_lep_reco[3] - err_lep_reco[3]) * (SF_lep_sel[3] - err_lep_sel[3]) * SF_lep_trig[0] * SF_lep_reco[0] * SF_lep_sel[0] * SF_lep_trig[1] * SF_lep_reco[1] * SF_lep_sel[1]
         
         yield_2e2mu_e_up  += weight_nom/SF_tot_nom * SF_var_e_up
         yield_2e2mu_mu_up += weight_nom/SF_tot_nom * SF_var_mu_up
            
         yield_2e2mu_e_dn  += weight_nom/SF_tot_nom * SF_var_e_dn
         yield_2e2mu_mu_dn += weight_nom/SF_tot_nom * SF_var_mu_dn

#   print "4e"
#   print "Nominal yield = {} Yield up = {} Yields down = {}".format(nom_yield_4e,yield_4e_up,yield_4e_dn)
#   print "\n"
#   print "4mu"
#   print "Nominal yield = {} Yield up = {} Yields down = {}".format(nom_yield_4mu,yield_4mu_up,yield_4mu_dn)
#   print "\n"
#   print "2e2mu"
#   print "Muons"
#   print "Nominal yield = {} Yield up = {} Yields down = {}".format(nom_yield_2e2mu,yield_2e2mu_mu_up,yield_2e2mu_mu_dn)
#   print "Electrons"
#   print "Nominal yield = {} Yield up = {} Yields down = {}".format(nom_yield_2e2mu,yield_2e2mu_e_up,yield_2e2mu_e_dn)
#   print "\n"

   print "LEPTON UNCERTAINTIES"
   print "4e = {:.1f} % 4mu = {:.1f}%".format(((yield_4e_up-nom_yield_4e)/nom_yield_4e*100),((yield_4mu_up-nom_yield_4mu)/nom_yield_4mu*100))
   print "2e2mu/ele = {:.1f} % 2e2mu/mu = {:.1f} %".format(((yield_2e2mu_e_up-nom_yield_2e2mu)/nom_yield_2e2mu*100),((yield_2e2mu_mu_up-nom_yield_2e2mu)/nom_yield_2e2mu*100))



