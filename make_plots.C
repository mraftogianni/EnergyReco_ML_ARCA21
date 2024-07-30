#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TF1.h"
#include "TMath.h"   
#if not defined(__CINT__) || defined(__MAKECINT__)
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#endif

double median1(TH1D* h1){
  double x[200]={0}; double y[200]={0}; double median=0;
 //compute the median for 1-d histogram h1
 int nbins2 = h1->GetXaxis()->GetNbins();
 //cout<<"nbins2: "<<nbins2<<endl;
 for(int i=0; i<nbins2; i++){
    x[i] = h1->GetBinCenter(i+1);
    y[i] = h1->GetBinContent(i+1);
    //cout<<"x: "<<x[i]<<" y: "<<y[i]<<endl;
  }
  median = TMath::Median(nbins2, x, y);
  return median;
}
/*  
void median2(TH2D* h2){
   //compute and print the median for each slice along X of h
   int nbins = h2->GetXaxis()->GetNbins();
   for(int i=0; i<nbins; i++){
      h1=h2->ProjectionY("",i,i+1);
      median = median1(h1);
      mean = h1->GetMean();
      cout<<"i: "<<i<<" median: "<<median<<" mean: "<<mean<<endl;
      delete h1;
    }
}
*/
  
void  make_plots()
{
  //input filename
  TFile* file = new TFile("NuMuCC_ML_ARCA21_v8.1_ml_13522_193_preproc.root"); 

  //input tree
  TTree *regTree = (TTree*)file->Get("ProcessedEvents");
  Float_t logEdepos, Log_No_PMTs_with_dist_weights_MAXhits, No_OMs_length_potential_length_chosen, PMTs_hit_to_nohit, MaxLenPos_OvrMaxDist, TrLenIT_3_OvrMaxDist, ToT_trig;

  regTree->SetBranchAddress("logEdepos", &logEdepos);
  regTree->SetBranchAddress("Log_No_PMTs_with_dist_weights_MAXhits", &Log_No_PMTs_with_dist_weights_MAXhits);
  regTree->SetBranchAddress("No_OMs_length_potential_length_chosen", &No_OMs_length_potential_length_chosen);
  regTree->SetBranchAddress("PMTs_hit_to_nohit", &PMTs_hit_to_nohit);
  regTree->SetBranchAddress("MaxLenPos_OvrMaxDist", &MaxLenPos_OvrMaxDist);
  regTree->SetBranchAddress("TrLenIT_3_OvrMaxDist", &TrLenIT_3_OvrMaxDist);
  regTree->SetBranchAddress("ToT_trig", &ToT_trig);

  //define plots
  TH2F *Log_emuMC_No_PMTs_with_dist_weights_MAXhits = new TH2F("plotLog_emuMC_No_PMTs_with_dist_weights_MAXhits_chosen","MC muon energy vs Number of PMTs with Pulses per Track (Adding weights from OM Distance);Log(MC deposited muon energy (GeV));Log(Number of PMTs with Pulses);",70, 1., 8., 40,0.,4.); 
  TH2F *emuNo_OMs_length_potential_length_chosen = new TH2F("plotemuNo_OMs_length_potential_length_chosen","MC muon energy vs Number of PMTs with Pulses per Track (Adding weights from OM Distance);Log(MC deposited muon energy (GeV));Log(Number of PMTs with Pulses);",70, 1., 8., 130,0.,1300.); 
    TH2F *emuPMTs_hit_to_nohit = new TH2F("plotemuPMTs_hit_to_nohit","MC muon energy vs Number of PMTs with Pulses per Track (Adding weights from OM Distance);Log(MC deposited muon energy (GeV));Log(Number of PMTs with Pulses);",70, 1., 8., 30,0.,3.);
    TH2F *emuMaxLenPos_OvrMaxDist = new TH2F("plotemuMaxLenPos_OvrMaxDist","MC muon energy vs Number of PMTs with Pulses per Track (Adding weights from OM Distance);Log(MC deposited muon energy (GeV));Log(Number of PMTs with Pulses);",70, 1., 8., 30,0.,0.0003);
    TH2F *emuTrLenIT_3_OvrMaxDist = new TH2F("plotemuTrLenIT_3_OvrMaxDist","MC muon energy vs Number of PMTs with Pulses per Track (Adding weights from OM Distance);Log(MC deposited muon energy (GeV));Log(Number of PMTs with Pulses);",70, 1., 8., 2,0.,0.2);
    TH2F *emuToT_trig = new TH2F("plotemuToT_trig","MC muon energy vs Number of PMTs with Pulses per Track (Adding weights from OM Distance);Log(MC deposited muon energy (GeV));Log(Number of PMTs with Pulses);",70, 1., 8., 2,0.,);

  //get entries form .root file
  for (Long64_t ievt=0; ievt<regTree->GetEntries();ievt++) {
     if (ievt%1000 == 0) {
        std::cout << "--- ... Processing event: "<<ievt<<" logEdepos= "<<logEdepos<<std::endl;
     }
     regTree->GetEntry(ievt);
       
     Log_emuMC_No_PMTs_with_dist_weights_MAXhits -> Fill(logEdepos, Log_No_PMTs_with_dist_weights_MAXhits);
     emuNo_OMs_length_potential_length_chosen -> Fill(logEdepos, No_OMs_length_potential_length_chosen);
     emuPMTs_hit_to_nohit -> Fill(logEdepos, PMTs_hit_to_nohit);
     emuMaxLenPos_OvrMaxDist -> Fill(logEdepos, MaxLenPos_OvrMaxDist);
     emuTrLenIT_3_OvrMaxDist -> Fill(logEdepos, TrLenIT_3_OvrMaxDist);
     emuToT_trig -> Fill(logEdepos,ToT_trig);
    
  }//end of entries
  //new TCanvas();


  //save plots in output file:
  TFile *outFile = new TFile("histogramsForEnergy.root", "RECREATE");

  // Write the histogram to the file
  Log_emuMC_No_PMTs_with_dist_weights_MAXhits -> Write(); 
  emuNo_OMs_length_potential_length_chosen -> Write();
  emuPMTs_hit_to_nohit -> Write();
  emuMaxLenPos_OvrMaxDist -> Write();
  emuTrLenIT_3_OvrMaxDist -> Write();
  emuToT_trig -> Write();

  outFile -> Close();  

}
