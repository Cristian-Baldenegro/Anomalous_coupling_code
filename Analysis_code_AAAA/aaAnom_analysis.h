#include "aaAnom_tree.h" // Here is where the TTree is fed in.
#include "TStopwatch.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <math.h>
#include <cmath>
#include <string>
#include "TROOT.h"
#include "TTree.h"
#include "TH2.h"
#include "TH1.h"
#include "TString.h"
#include "TLegend.h"
#include "TSystem.h"
#include "TAxis.h"
#include "TLatex.h"
#include "TLorentzVector.h"
#include "TMinuit.h"
#include "TFile.h"            
#include "TCanvas.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TF1.h"
#include "TStyle.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TSlider.h"
#include <iomanip>
#include "TGaxis.h"

void aaAnom_tree::Loop(){};

using namespace std;

class aaAnom_analysis : public aaAnom_tree 
{
   public : 

   int part_id;
   int pu_flag;
   int event_count;
   int dimGen;
   int high_flag;

   double norm_f;
   double mass_window;

   TStopwatch m_timer;
   double rtime;
   double ctime;


  // dummy constructor //
  aaAnom_analysis(TTree * tree): aaAnom_tree(tree){}
  //--------Calculation of physical quantities----//
  //--------Rapidities, transverse momentum, azimuthal angle
  double CalcE(double e)
  {
    double eout = e*gRandom->Gaus(1,0.01);
    return eout;
  } 
 
  double CalcEta(double pz, double e)
  {
    double eta = 0.5*log((e + pz)/(e - pz));
    eta = gRandom->Gaus(eta,0.001);
    return eta;
  }
    
  double CalcPhi(double px, double py)
  {
     double phi = 0;
     if(py>0)  phi = atan2(py,px); 
     else phi = atan2(py,px) + 2*3.14159;
     phi = gRandom->Gaus(phi,0.001);
     return phi;
  }
 
   double CalcPt(double e, double eta)
  {
    double pt = e/cosh(eta);
    return pt;
  }

 
  //-------analysis methods---------//
  
  //initialize the constant attribute
  //chan=0 photons, chan=1 electrons, chan=2 jets
  
  void initialize(double nf,int chan, int pu, int signal_flag,int hf)
  {//Starting the watch
    cout << "*****************************************" << endl;
  
    m_timer.Start(); //Start the clock
  
    //Defining class attributes
    if(signal_flag!=1)  norm_f=nf/*/40*/;   //VERTEX CUT /40 at mu=50
    else                norm_f=nf;

    event_count=1;
    dimGen=50000;
    mass_window=0.03;
    high_flag=hf; 
  
    if(chan==0) part_id=22; //photon regular case
  
    else if(chan==1)
    { 
      part_id=11;
      norm_f*=0.01*0.01; // fake electron->photon  2ELECTRONS!
    }
  
    else if(chan==2)
    {
      part_id=999;   //  jets
    }
    
    if(pu==1)  pu_flag=1;    //pile-up background
    else       pu_flag=0;
  }
  
  //if not specified, double tag is assumed
  void initialize(double nf, int chan)
  {
    this->initialize(nf,chan,0,0,0);
  }
  
  
  //if not specified, photons are spotted and double tag is assumed
  void initialize(double nf)
  {
    this->initialize(nf,0,0,0,0);
  }
  
  //if not specified, signal_flag==0 is assumed
  void initialize(double nf,int chan, int pu)
  {
    this->initialize(nf,chan,pu,0,0);
  }
  
//------------------------------------------------------------------------------//  
  //Loop on the events and on the particles
  void execute(TFile * out) //This will be the output root file
  {
  
    if (fChain == 0) return;
  
     Long64_t nentries = fChain->GetEntriesFast();
     Long64_t nbytes = 0, nb = 0;
  
 
     //1-Dim Histograms definitions
  
     TH1F * h_eta = new TH1F("h_eta","",1000,-5,5);
     TH1F * h_delta_eta = new TH1F("h_delta_eta","",1000,0,10);

     TH1F * h_phi = new TH1F("h_phi","",1000,-0.5,6.5);
     TH1F * h_delta_phi = new TH1F("h_delta_phi","",6500,-3.15,3.15);

     TH1F * h_pt1 = new TH1F("h_pt1","",1000,0,1000);
     TH1F * h_pt2 = new TH1F("h_pt2","",1000,0,1000);
     TH1F * h_pt1_no = new TH1F("h_pt1_no","",200,0,2000);
     TH1F * h_pt2_no = new TH1F("h_pt2_no","",200,0,2000);
     TH1F * h_pt_ratio = new TH1F("h_pt_ratio","",1000,0,1);

     TH1F * h_w = new TH1F("h_w","",2200,0,2200);

     TH1F * h_na = new TH1F("h_na","",11,-0.5,10.5);
    
     TH1F * h_xi_no = new TH1F("h_xi_no","",1000,0,1);
     TH1F * h_xi = new TH1F("h_xi","",1000,0,1);
     TH1F * h_w_diff_no = new TH1F("h_w_diff_no","",1000,0,5000);
     TH1F * h_w_diff = new TH1F("h_w_diff","",5000,0,5000);
  
     TH1F * h_w_ratio = new TH1F("h_w_ratio","",10000,0,10);
  
     TH1F * h_rapid = new TH1F ("h_rapid","",40000,-2,2); 

     TH1F * h_mcweight = new TH1F("h_mcweight","",10000,0,1);

     TH1F * h_Npu = new TH1F("h_Npu","",500,0,500);
     TH1F * h_Np_left = new TH1F("h_Np_left","",500,0,500);
     TH1F * h_Np_right = new TH1F("h_Np_right","",500,0,500);
     
     //Loop on the events
  
     for (Long64_t jentry=0; jentry<nentries;jentry++) 
     {
  
  
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
        // if (Cut(ientry) < 0) continue;
  
  
        //declaration and initialization of local variables
   
        double eta1=0;
        double eta2=0;
        double delta_eta=0;

        double phi1=0;
        double phi2=0;
        double delta_phi=0;

        double pt1=0;
        double pt2=0;
        double pt_ratio=0;

        double w=0; //Invariant mass

        double xi1=0;
        double xi2=0;

        double w_diff=0;

        int a_index=0;

        double norm_f_conv=norm_f; 
        double norm_f_eff=1;
  
        double e1=0;
        double e2=0;

        double e1_raw=0;
        double e2_raw=0;

        double pz1=0;
        double pz2=0;

        double rapid=0; 
        vector<double> xi1_vect;
        vector<double> xi2_vect;

        //Loop on the particles
  
        if(ngen>dimGen) cout << "WARNING: EVENT CUT" << endl;
  
        for(int i=0; i<ngen; i++) //This loop is to calculate basic kinematical variables
        {
  
          if((fabs(id[i])==part_id || id[i]==94 /*Gamma case*/ ) && a_index==0)
          {
            eta1 = CalcEta(pz[i],e[i]);
            phi1 = CalcPhi(px[i],py[i]);
            e1 = CalcE(e[i]);
            pt1 = CalcPt(e1,eta1);      

            e1_raw=e[i]; 
            pz1 = pz[i];

            a_index++;
          }
  
          else if((fabs(id[i])==part_id || id[i]==94) && a_index>0)
          {
  
            if(pt1<CalcPt(e[i],CalcEta(pz[i],e[i])))
            {
              eta2 = eta1;
              phi2 = phi1;
              e2 = e1;
              pt2 = pt1;     
  
              e2_raw = e1_raw;
              pz2=pz1;
             

              eta1 = CalcEta(pz[i],e[i]);
              phi1 = CalcPhi(px[i],py[i]);         
              e1 = CalcE(e[i]);

              pt1 = CalcPt(e1,eta1);    

              e1_raw = e[i];   
              pz1 = pz[i];
            }
  
            else if(pt2<CalcPt(e[i],CalcEta(pz[i],e[i])))
            {
              eta2 = CalcEta(pz[i],e[i]);
              phi2 = CalcPhi(px[i],py[i]);
              e2 = CalcE(e[i]);

              pt2 = CalcPt(e2,eta2);  

              e2_raw = e[i];     
              pz2 = pz[i];
            }
  
            a_index++;
  
          } //endif
        }  // end loop on particles
  
///////////////////////////////////////////////////////////////////////////////              
        //diffractive mass calculation - no pile-up case
        if(pu_flag!=1)
        {
          xi1 = (6500-e[0])/6500; 
          xi1 *= gRandom->Gaus(1,0.02); 
          xi2 = (6500-e[1])/6500; 
          xi2 *= gRandom->Gaus(1,0.02); 
          w_diff = sqrt(xi1*xi2)*13000;  //GeV

          h_xi_no->Fill(xi1,norm_f);
          h_xi_no->Fill(xi2,norm_f);
          h_w_diff_no->Fill(w_diff,norm_f);
        }
        
        //-------PHOTON OR ELECTRON CASE-------//
        if(part_id!=999)
        {
 
          //Warning and rejection if less than 2 central object 
          //(photons or electrons)
          if(a_index<2)
          {
            cout << "WARNING: EVENT " << event_count << 
            " WITH LESS THAN 2 PHOTONS/ELECTRONS" << endl;
            pt1=-999;
            pt2=-999;
          }

          //rejecting bad events from low mass DY prod 
          if(high_flag==1 && pu_flag==1 && part_id==11 && w>300)
          {
            pt1=-999;
            pt2=-999;
          }

          //rejecting bad events low pt ND diphoton production
          if(high_flag==1 && pu_flag==1 && part_id==22 && pt1>200)
          {
            pt1=-999;
            pt2=-999;
          } 
 
          //diphoton mass
          w = sqrt(2*pt1*pt2*(cosh(eta1-eta2)-cos(phi1-phi2)));

          //pT before any selection
          h_pt1_no->Fill(pt1,norm_f);
          h_pt2_no->Fill(pt2,norm_f);
  
          //pt ratio 
          pt_ratio = pt2/pt1;
 
          //delta_eta and delta_phi
          delta_eta = fabs(eta1-eta2);    
          if(fabs(phi1-phi2)<3.14159265359) delta_phi = fabs(phi1-phi2);
          else delta_phi = 2*3.14159265359 - fabs(phi1-phi2);
          //Rapidity of the system
          rapid = CalcEta(pz1+pz2,e1_raw+e2_raw); //eta1+eta2;          

          //for pile-up, selection of the xi1/xi2 which are the closest to
          //the diphoton mass/rap
          if(pu_flag==1)
          {

            double w_ratio_test = 999;
            double rapid_test = 999;

            if(xi1_vect.size()>0 && xi2_vect.size()>0)
            {

              for(int l=0; l<xi1_vect.size(); l++)
              {   
                for(int m=0; m<xi2_vect.size(); m++)
                {
                  double test1=fabs(1-sqrt(xi1_vect[l]*xi2_vect[m])*14000/w);
                  double test2=fabs(rapid-0.5*log(xi1_vect[l]/xi2_vect[m]));

                  //selecting the good protons on mass/rap criteria
                  if(pow(test1,2)+pow(test2,2) < 
                     pow(w_ratio_test,2)+pow(rapid_test,2))
                  {
                    xi1=xi1_vect[l];
                    xi2=xi2_vect[m];
                    w_diff = sqrt(xi1*xi2)*14000;
                    w_ratio_test = 1-w_diff/w;
                    rapid_test = rapid-0.5*log(xi1/xi2);
                  }

                }
              }
            }

            else
            {
              xi1=-999;
              xi2=-999;
            }
          }


          //selection: pT, eta, xi_1, xi_2, deltaphi, ptratio
          if(pt1>200 && pt2>100 && w>600 && fabs(eta1)<2.37 && 
             fabs(eta2)<2.37 && max(xi1,xi2)<0.15 && 
             min(xi1,xi2)>0.015 && fabs(3.14159265359-delta_phi)<0.01 && 
             pt_ratio>0.95 /*&& w_diff/w<1+mass_window && 
             w_diff/w>1-mass_window && fabs(rapid-0.5*log(xi1/xi2))<0.03*/)
          {        
    
             //converted photon factor as an eta function
             //we request 1 converted photon
             //for electrons, already converted ...
             if(part_id==22)
             {
               if(fabs(eta1)<1.4 && fabs(eta2)<1.4)
               { 
                 norm_f_conv *= (1-(0.85*0.85)); 
               }
               else if(fabs(eta1)>=1.4 && fabs(eta2)>=1.4)
               {
                 norm_f_conv *= (1-(0.7*0.7)); 
               }
               else norm_f_conv *= (1-(0.7*0.85));
             }
   
             //efficiency factor for real photons-electrons
             norm_f_eff=norm_f_conv*(0.76-1.98*exp(-pt1/16.1))*
                                    (0.76-1.98*exp(-pt2/16.1));
    
 
              //Filling histos
             h_eta->Fill(eta1,norm_f_eff);
             h_eta->Fill(eta2,norm_f_eff);
             h_delta_eta->Fill(delta_eta,norm_f_eff);
   
             h_phi->Fill(phi1,norm_f_eff);
             h_phi->Fill(phi2,norm_f_eff);
             h_delta_phi->Fill(3.14159265359-delta_phi,norm_f_eff);
   
             h_pt1->Fill(pt1,norm_f_eff);
             h_pt2->Fill(pt2,norm_f_eff);
             h_pt_ratio->Fill(1-pt_ratio,norm_f_eff);
   
             h_w->Fill(w,norm_f_eff);
             h_na->Fill(a_index,norm_f_eff);
   
             h_xi->Fill(xi1,norm_f_eff);
             h_xi->Fill(xi2,norm_f_eff);
             h_w_diff->Fill(w_diff,norm_f_eff);
   
             h_w_ratio->Fill(w_diff/w,norm_f_eff);
  
             h_rapid->Fill(rapid-0.5*log(xi1/xi2),norm_f_eff);
            
             h_mcweight->Fill(norm_f);
 
          } //endif 2nd part of selection: pT, eta, xi_1, xi_2, deltaphi, ptratio
        } //endif parti_id!=999
       event_count++;
      }//Ends event loop
  
       //Compute MC error
  
       cout << "NORM MC TO 300 FB-1 (Wrong for jets): " << norm_f << endl;
  
       out->cd();
       out->Write(); //Here's where the sweet magic happens.
  
       m_timer.Stop();
       rtime = m_timer.RealTime();
       ctime = m_timer.CpuTime();
    
       cout << "********************************" << endl;
       cout << "INFORMATION :: Performance : " << endl;
       cout << "RealTime= " << rtime << "seconds, CpuTime= " << ctime << 
                  "seconds" << endl;
    }
  
     private:
  
};
