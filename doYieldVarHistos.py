import ROOT
ROOT.gROOT.SetBatch(True)
from ROOT import *

gRandom.SetSeed(101)


lumi = 35.876 # Set lumi to be used for MC scaling
n_toys = 200 # Number of toys to be generated for each lepton

folder = './Moriond_2017_v2/' # Define input folder name
file_name = '/ZZ4lAnalysis.root' # Define input dile name
output_file_name = "Yield_distributions.root"

# List of samples to run on
List = [
'ggH125',
#'ZZTo4l',
#'ggTo4l'
]

output=TFile.Open(output_file_name, "RECREATE") # root file to save histograms

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
   
   #Declare histograms where we store varied yields
   nbins_yields = 100
   histo_width = 0.15
   
   yield_4mu = TH1F("yield_4mu_" + Type, "yield_4mu_" + Type, nbins_yields, ( 1 - histo_width ) * nom_yield_4mu , ( 1 + histo_width ) * nom_yield_4mu )

   yield_4e = TH1F("yield_4e_" + Type, "yield_4e_" + Type, nbins_yields   , ( 1 - histo_width ) * nom_yield_4e,   ( 1 + histo_width ) * nom_yield_4e )
   
   yield_2e2mu_e = TH1F("yield_2e2mu_e_" + Type, "yield_2e2mu_e_" + Type, nbins_yields,    ( 1 - histo_width ) * nom_yield_2e2mu, ( 1 + histo_width ) * nom_yield_2e2mu )
   
   yield_2e2mu_mu = TH1F("yield_2e2mu_mu_" + Type, "yield_2e2mu_mu_" + Type, nbins_yields, ( 1 - histo_width ) * nom_yield_2e2mu, ( 1 + histo_width ) * nom_yield_2e2mu )
   
   #control_histo_nom = TH1F("Lep_sf_nom_" + Type, "Lep_sf_nom_" + Type, 100, ( 1 - histo_width )  , ( 1 + histo_width )  )
   #control_histo_var = TH1F("Lep_sf_var_" + Type, "Lep_sf_var_" + Type, 100, ( 1 - histo_width )  , ( 1 + histo_width )  )

   print "Processing {} ...".format(Type)
   br_data = 0
   
   # Set yields for toys to zero
   yield_4e_sum = []
   yield_4mu_sum = []
   yield_2e2mu_e_sum = []
   yield_2e2mu_mu_sum = []
   var_trig = []
   var_reco = []
   var_sel = []
   
   for i_toy in range(0,n_toys):
      yield_4e_sum.append(0.)
      yield_4mu_sum.append(0.)
      yield_2e2mu_e_sum.append(0.)
      yield_2e2mu_mu_sum.append(0.)
      var_trig.append(gRandom.Gaus(0.,1.));
      var_reco.append(gRandom.Gaus(0.,1.));
      var_sel.append(gRandom.Gaus(0.,1.));

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

      # Produce toy MC and add to integral yields
      for i_toy in range(0,n_toys):
         SF_var = 1.
         SF_var_e = 1.
         SF_var_mu = 1.
         
         for i in range (0,4):
            # Read lepton SF and unc from tree
            SF_lep_trig[i] = 1. #hard-code value for trigger SF at the moment
            err_lep_trig[i] = 0. #hard-code value for trigger unc at the moment
            SF_lep_reco[i] = event.LepRecoSF[i]
            err_lep_reco[i] = event.LepRecoSF_Unc[i]
            SF_lep_sel[i] = event.LepSelSF[i]
            err_lep_sel[i] = event.LepSelSF_Unc[i]
         
         
         if (idL1==11 and idL3==11):
            # Vary SF of each lepton using gaus distribution with mean as SF nominal value and sigma as SF unc
            #control_histo_nom.Fill(SF_lep_sel[0])
            #control_histo_var.Fill(gRandom.Gaus(SF_lep_sel[0], err_lep_sel[0]))
            for i in range (0,4):
               SF_var *= (SF_lep_trig[i] + var_trig[i_toy]*err_lep_trig[i]) * (SF_lep_reco[i] +var_reco[i_toy]*err_lep_reco[i]) * (SF_lep_sel[i] + var_sel[i_toy]*err_lep_sel[i])
            yield_4e_sum[i_toy] += weight_nom/SF_tot_nom * SF_var

         elif (idL1==13 and idL3==13):
            for i in range (0,4):
               SF_var *= (SF_lep_trig[i] + var_trig[i_toy]*err_lep_trig[i]) * (SF_lep_reco[i] +var_reco[i_toy]*err_lep_reco[i]) * (SF_lep_sel[i] + var_sel[i_toy]*err_lep_sel[i])
            yield_4mu_sum[i_toy] += weight_nom/SF_tot_nom * SF_var

         elif (abs(idL1-idL3)==2):
            # Vary electron SF while fixing muon and vice-versa
            if ( idL1 == 11):
               SF_var_e  = (SF_lep_trig[0] + var_trig[i_toy]*err_lep_trig[0]) * (SF_lep_reco[0] + var_reco[i_toy]*err_lep_reco[0]) * (SF_lep_sel[0] + var_sel[i_toy]*err_lep_sel[0]) * (SF_lep_trig[1] + var_trig[i_toy]*err_lep_trig[1]) * (SF_lep_reco[1] + var_reco[i_toy]*err_lep_reco[1]) * (SF_lep_sel[1] + var_sel[i_toy]*err_lep_sel[1]) * SF_lep_trig[2] * SF_lep_reco[2] * SF_lep_sel[2] * SF_lep_trig[3] * SF_lep_reco[3] * SF_lep_sel[3]
               SF_var_mu = (SF_lep_trig[2] + var_trig[i_toy]*err_lep_trig[2]) * (SF_lep_reco[2] + var_reco[i_toy]*err_lep_reco[2]) * (SF_lep_sel[2] + var_sel[i_toy]*err_lep_sel[2]) * (SF_lep_trig[3] + var_trig[i_toy]*err_lep_trig[3]) * (SF_lep_reco[3] + var_reco[i_toy]*err_lep_reco[3]) * (SF_lep_sel[3] + var_sel[i_toy]*err_lep_sel[3]) * SF_lep_trig[0] * SF_lep_reco[0] * SF_lep_sel[0] * SF_lep_trig[1] * SF_lep_reco[1] * SF_lep_sel[1]
            
            elif ( idL1 == 13):
               SF_var_mu = (SF_lep_trig[0] + var_trig[i_toy]*err_lep_trig[0]) * (SF_lep_reco[0] + var_reco[i_toy]*err_lep_reco[0]) * (SF_lep_sel[0] + var_sel[i_toy]*err_lep_sel[0]) * (SF_lep_trig[1] + var_trig[i_toy]*err_lep_trig[1]) * (SF_lep_reco[1] + var_reco[i_toy]*err_lep_reco[1]) * (SF_lep_sel[1] + var_sel[i_toy]*err_lep_sel[1]) * SF_lep_trig[2] * SF_lep_reco[2] * SF_lep_sel[2] * SF_lep_trig[3] * SF_lep_reco[3] * SF_lep_sel[3]
               SF_var_e  = (SF_lep_trig[2] + var_trig[i_toy]*err_lep_trig[2]) * (SF_lep_reco[2] + var_reco[i_toy]*err_lep_reco[2]) * (SF_lep_sel[2] + var_sel[i_toy]*err_lep_sel[2]) * (SF_lep_trig[3] + var_trig[i_toy]*err_lep_trig[3]) * (SF_lep_reco[3] + var_reco[i_toy]*err_lep_reco[3]) * (SF_lep_sel[3] + var_sel[i_toy]*err_lep_sel[3]) * SF_lep_trig[0] * SF_lep_reco[0] * SF_lep_sel[0] * SF_lep_trig[1] * SF_lep_reco[1] * SF_lep_sel[1]
            
            yield_2e2mu_e_sum[i_toy] += weight_nom/SF_tot_nom * SF_var_e
            yield_2e2mu_mu_sum[i_toy] += weight_nom/SF_tot_nom * SF_var_mu
   
   #Save all toy yields in histograms
   for i_toy in range(0,n_toys):
      yield_4e.Fill(yield_4e_sum[i_toy])
      yield_4mu.Fill(yield_4mu_sum[i_toy])
      yield_2e2mu_e.Fill(yield_2e2mu_e_sum[i_toy])
      yield_2e2mu_mu.Fill(yield_2e2mu_mu_sum[i_toy])


   output.cd()
   yield_4e.Write()
   yield_4mu.Write()
   yield_2e2mu_e.Write()
   yield_2e2mu_mu.Write()
   #control_histo_nom.Write()
   #control_histo_var.Write()

   print "Processing of {} finished.".format(Type)
   print "\n"

output.Close()
print "Done! Everything saved in {}.".format(output_file_name)
