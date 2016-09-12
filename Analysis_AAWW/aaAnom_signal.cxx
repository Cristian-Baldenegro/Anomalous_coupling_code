#include "aaAnom_tree.h"
#define aaAnom_tree_cxx
#include "aaAnom_analysis.h"
#include "TFileCollection.h"
#include <iostream>


int main(){
  double L = 300000; //Expected Integrated Luminosity in pb-1
  int N = 100000; //Events Generated 
  //-------SIGNAL-------//
  cout << "Signal analysis " << endl; 

//Buffer is the input root-file name
//Bufferout is the output root-file name

  TFile * f_sgn = new TFile("AAZZ_anom.root","READ");


  TH1F * h1 = (TH1F*)f_sgn->Get("h999");
  double sigma_sgn = h1->GetBinContent(1);  //Gets the cross-section in pb
  f_sgn->Close();
    
  TFileCollection * fc_sgn = new TFileCollection();
  fc_sgn->Add("AAZZ_anom.root");
    
  TChain * c_sgn = new TChain("h777");
  c_sgn->AddFileInfoList((TCollection*)fc_sgn->GetList());
    
  aaAnom_analysis * a_sgn = new aaAnom_analysis(c_sgn); //Feed in a TTree/TChain
  TFile * out_sgn = new TFile("AAZZ_anom_pu.root","recreate");

//  void initialize(normalization factor, channel, pileup flag, signal_flag flag, high flag)
  a_sgn->initialize(sigma_sgn*L/N,0,1,1,0);
  a_sgn->execute(out_sgn);

  out_sgn->Close();
}
