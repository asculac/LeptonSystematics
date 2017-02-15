import ROOT
ROOT.gROOT.SetBatch(True)
from ROOT import *

gRandom.SetSeed(101)


lumi = 35.8
n_events = 10000 #number of events in tree to run on
n_toys = 10 #number of toys to be generated for each lepton

folder = './Moriond_2017/'
file_name = '/ZZ4lAnalysis.root'

List = [
'ggH125', #signal sample
#'ZZTo4l', #backgroung sample
#'ggTo4l'
]

for Type in List:
   
   mc=TFile.Open(folder+Type+file_name)
   tree = mc.Get("ZZTree/candTree")
   counters = mc.Get("ZZTree/Counters")
   NGen = counters.GetBinContent(40)

   nbins=35 #number of bins in histogram
   lo=105   #Mass range for yields
   hi=140

   #Declare histogram to calculate nominal value of yields
   h_nom = TH1F("h_nom", "h_nom", nbins, lo, hi)
   
   #Declare histograms where we store yields
   nbins_yields = 100
   tree.Draw("ZZMass >> h_nom" , "(abs(LepLepId[0]) == 11 && abs(LepLepId[3]) == 11)*overallEventWeight*1000*xsec*"+str(lumi)+"/"+str(NGen))
   nom_yield_4e = h_nom.Integral()
   
   tree.Draw("ZZMass >> h_nom" , "(abs(LepLepId[0]) == 13 && abs(LepLepId[3]) == 13)*overallEventWeight*1000*xsec*"+str(lumi)+"/"+str(NGen))
   nom_yield_4mu = h_nom.Integral()
   
   tree.Draw("ZZMass >> h_nom" , "(abs(abs(LepLepId[0]) - (LepLepId[3])) == 2)      *overallEventWeight*1000*xsec*"+str(lumi)+"/"+str(NGen))
   nom_yield_2e2mu = h_nom.Integral()
   
   yield_4mu = TH1F("yield_4mu", "yield_4mu", nbins_yields, nom_yield_4mu - 0.5, nom_yield_4mu + 0.5)

   yield_4e = TH1F("yield_4e", "yield_4e", nbins_yields, nom_yield_4e - 0.5, nom_yield_4e + 0.5)
   
   yield_2e2mu_e = TH1F("yield_2e2mu_e", "yield_2e2mu_e", nbins_yields, nom_yield_2e2mu - 0.5, nom_yield_2e2mu + 0.5)
   
   yield_2e2mu_mu = TH1F("yield_2e2mu_mu", "yield_2e2mu_mu", nbins_yields, nom_yield_2e2mu - 0.5, nom_yield_2e2mu + 0.5)

   print "Processing {} ...".format(Type)
   br_data = 0
   
   #Set yields for toys to zero
   yield_4e_sum = []
   yield_4mu_sum = []
   yield_2e2mu_e_sum = []
   yield_2e2mu_mu_sum = []
   
   for i_toy in range(0,n_toys):
      yield_4e_sum.append(0.)
      yield_4mu_sum.append(0.)
      yield_2e2mu_e_sum.append(0.)
      yield_2e2mu_mu_sum.append(0.)

   for event in tree:#loop over all events in tree
      br_data+=1
      if(br_data % int((tree.GetEntries()/10)) == 0):
         print "{} %".format(str(100*br_data/tree.GetEntries() + 1))
      
      mass4l = tree.ZZMass
      if( mass4l < 105. or mass4l > 140.): continue #Skip events that are not in the mass window
      idL1 = abs(tree.LepLepId[0])
      idL3 = abs(tree.LepLepId[3])
    
      GENmassZZ = tree.GenHMass
      
      #Calculate nominal weigh using central value of SF
      weight_nom = event.overallEventWeight*1000*lumi*event.xsec/NGen
      
      # Produce toy MC and add to integral yields
      for i_toy in range(0,n_toys):
         if (idL1==11 and idL3==11):
            yield_4e_sum[i_toy] += weight_nom/event.dataMCWeight * gRandom.Gaus(event.dataMCWeight, 0.1)
         elif (idL1==13 and idL3==13):
            yield_4mu_sum[i_toy] += weight_nom/event.dataMCWeight * gRandom.Gaus(event.dataMCWeight, 0.1)
         elif (abs(idL1-idL3)==2):
            yield_2e2mu_e_sum[i_toy] += weight_nom/event.dataMCWeight * gRandom.Gaus(event.dataMCWeight, 0.1)
            yield_2e2mu_mu_sum[i_toy] += weight_nom/event.dataMCWeight * gRandom.Gaus(event.dataMCWeight, 0.1)
               
   #Save all toy yields in histograms
   for i_toy in range(0,n_toys):
      yield_4e.Fill(yield_4e_sum[i_toy])
      yield_4mu.Fill(yield_4mu_sum[i_toy])
      yield_2e2mu_e.Fill(yield_2e2mu_e_sum[i_toy])
      yield_2e2mu_mu.Fill(yield_2e2mu_mu_sum[i_toy])

#      recoerre=1.0
#      recoerrmu=1.0
#
#    # Reconstruction efficiencies errors
#    for l in range(0,4):
#      if abs(tree.lep_id[tree.lep_Hindex[l]])==11:
#        # uncertainty on ID/Iso/SIP, stored in tree per lepton
#        recoerre *= (1.0+gRandom.Gaus(0.0,abs(tree.lep_dataMCErr[tree.lep_Hindex[l]])))
#        # uncertainty on GSF track eff. for crack electrons, hard coded
#        if (abs(tree.lep_eta[tree.lep_Hindex[l]])>1.4442 and abs(tree.lep_eta[tree.lep_Hindex[l]])<1.566): 
#          recoerre *= (1.0+gRandom.Gaus(0.0,0.08))
#        # uncertainty on GSF track eff. for low pt non-crack electrons, hard coded
#        elif (tree.lep_pt[tree.lep_Hindex[l]]<20): 
#          recoerre *= (1.0+gRandom.Gaus(0.0,0.08))
#        # uncertainty on GSF track eff. for high pt non-crack electrons, hard coded
#        else:
#          recoerre *= (1.0+gRandom.Gaus(0.0,0.04))
#      else:
#        # total uncertainty on low pt muons, hard coded
#        if (tree.lep_pt[tree.lep_Hindex[l]]<10.0 or abs(tree.lep_eta[tree.lep_Hindex[l]])>2.0): recoerrmu *= (1.0+gRandom.Gaus(0.0,0.02))
#        # total uncertainty on high pt muons, hard coded
#        else: recoerrmu *= (1.0+gRandom.Gaus(0.0,0.005))    
#
#    if (recoerre<1.0): 
#      recoerrupe=1.0/recoerre
#      recoerrdne=recoerre
#    else:
#      recoerrupe=recoerre
#      recoerrdne=(1.0/recoerre)
#
#    if (recoerrmu<1.0): 
#      recoerrupmu=1.0/recoerrmu
#      recoerrdnmu=recoerrmu
#    else:
#      recoerrupmu=recoerrmu
#      recoerrdnmu=(1.0/recoerrmu)
#
#    # Trigger uncertainties, hard coded
#    if(idL1==11 and idL3==11):
#      recoerrupe = 1.0+sqrt((1.0-recoerrupe)**2+0.01**2)
#      recoerrdne = 1.0-sqrt((1.0-recoerrdne)**2+0.04**2)
#    if(idL1==13 and idL3==13):
#      recoerrupmu = 1.0+sqrt((1.0-recoerrupmu)**2+0.00**2)
#      recoerrdnmu = 1.0-sqrt((1.0-recoerrdnmu)**2+0.01**2)
#    if(abs(idL1-idL3)==2):
#      recoerrupe = 1.0+sqrt((1.0-recoerrupe)**2+0.01**2)
#      recoerrdne = 1.0-sqrt((1.0-recoerrdne)**2+0.01**2)
#      recoerrupmu = 1.0+sqrt((1.0-recoerrupmu)**2+0.01**2)
#      recoerrdnmu = 1.0-sqrt((1.0-recoerrdnmu)**2+0.01**2)
#
#    if (idL1==11 and idL3==11): 
#      hrecoup4e.Fill(mass4l,recoerrupe*dataMCWeight*genWeight*crossSection*(2763./sumweights))
#      hrecodn4e.Fill(mass4l,recoerrdne*dataMCWeight*genWeight*crossSection*(2763./sumweights))
#    elif (idL1==13 and idL3==13): 
#      hrecoup4mu.Fill(mass4l,recoerrupmu*dataMCWeight*genWeight*crossSection*(2763./sumweights))
#      hrecodn4mu.Fill(mass4l,recoerrdnmu*dataMCWeight*genWeight*crossSection*(2763./sumweights))
#    elif (abs(idL1-idL3)==2): 
#      hrecoup2e2mu_e.Fill(mass4l,recoerrupe*dataMCWeight*genWeight*crossSection*(2763./sumweights))
#      hrecodn2e2mu_e.Fill(mass4l,recoerrdne*dataMCWeight*genWeight*crossSection*(2763./sumweights))
#      hrecoup2e2mu_mu.Fill(mass4l,recoerrupmu*dataMCWeight*genWeight*crossSection*(2763./sumweights))
#      hrecodn2e2mu_mu.Fill(mass4l,recoerrdnmu*dataMCWeight*genWeight*crossSection*(2763./sumweights))
#
#  #print Type,'Reco/Id'
#  #print '4e +',round((hrecoup4e.Integral()/hnom4e.Integral()-1.0),3),'-',round((1.0-hrecodn4e.Integral()/hnom4e.Integral()),3)
#  #print '4mu +',round((hrecoup4mu.Integral()/hnom4mu.Integral()-1.0),3),'-',round((1.0-hrecodn4mu.Integral()/hnom4mu.Integral()),3)
#  #print '2e2mu +',round((hrecoup2e2mu.Integral()/hnom2e2mu.Integral()-1.0),3),'-',round((1.0-hrecodn2e2mu.Integral()/hnom2e2mu.Integral()),3)
#
#  print Type,'Reco/Id/Trig'
#  print '4e +',round((hrecoup4e.Integral()/hnom4e.Integral()-1.0),3),'-',round((1.0-hrecodn4e.Integral()/hnom4e.Integral()),3)
#  print '4mu +',round((hrecoup4mu.Integral()/hnom4mu.Integral()-1.0),3),'-',round((1.0-hrecodn4mu.Integral()/hnom4mu.Integral()),3)
#  print '2e2mu el +',round((hrecoup2e2mu_e.Integral()/hnom2e2mu.Integral()-1.0),3),'-',round((1.0-hrecodn2e2mu_e.Integral()/hnom2e2mu.Integral()),3)
#  print '2e2mu mu +',round((hrecoup2e2mu_mu.Integral()/hnom2e2mu.Integral()-1.0),3),'-',round((1.0-hrecodn2e2mu_mu.Integral()/hnom2e2mu.Integral()),3)
   output=TFile.Open(Type+"_Yield_var.root", "RECREATE")
   
   yield_4e.Write()
   yield_4mu.Write()
   yield_2e2mu_e.Write()
   yield_2e2mu_mu.Write()
   output.Close()
