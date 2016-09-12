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
#include "fastjet/ClusterSequence.hh"
void aaAnom_tree::Loop(){};

using namespace std;
using namespace fastjet;

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
  aaAnom_analysis(TTree * tree): aaAnom_tree(tree)
  {}


  //-------trivial calculation methods-----//
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
   double CalcWsq_X(double pt1, double pt2, double eta1, double eta2, double phi1, double phi2){
	double Wsq_X=2*pt1*pt2*(cosh(eta1-eta2)-cos(phi1-phi2));
	return Wsq_X;
  }

	double LorenzProd(double e1, double e2, double px1, double px2, double py1, double py2, double pz1, double pz2){
	double LorenzProd=e1*e2-px1*px2-py1*py2-pz1*pz2;
	return LorenzProd;
  }

 
  //-------analysis methods---------//
  
  //initialize the constant attribute
  //chan=0 photons, chan=1 electrons, chan=2 jets
  
  void initialize(double nf,int chan, int pu, int signal_flag,int hf)
  {
  
    //Starting the watch
    cout << "*****************************************" << endl;
  
    m_timer.Start();
  
    //Defining class attributes
    if(signal_flag!=1)  norm_f=nf/*/40*/;   //VERTEX CUT /40 at mu=50
    else                norm_f=nf;

    event_count=1;
    dimGen=50000;
    mass_window=0.1;
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
  
  
  
  
  //loop on the events and on the particles
  void execute(TFile * out)
  {
  
    if (fChain == 0) return;
  
     Long64_t nentries = fChain->GetEntriesFast();
     Long64_t nbytes = 0, nb = 0;
  
 
     //local histo
  
     TH1F * h_eta = new TH1F("h_eta","",100,-10,10);
     TH1F * h_delta_eta = new TH1F("h_delta_eta","",500,0,10);

     TH1F * h_phi = new TH1F("h_phi","",500,-0.5,6.5);
     TH1F * h_delta_phi = new TH1F("h_delta_phi","",300,-3.15,3.15);

     TH1F * h_pt1 = new TH1F("h_pt1","",100,0,1000);
     TH1F * h_pt2 = new TH1F("h_pt2","",100,0,1000);
     TH1F * h_pt1_no = new TH1F("h_pt1_no","",100,0,2000);
     TH1F * h_pt2_no = new TH1F("h_pt2_no","",100,0,2000);
     TH1F * h_pt_ratio = new TH1F("h_pt_ratio","",1000,0,1);

     TH1F * h_w = new TH1F("h_w","",100,0,1800);

     TH1F * h_na = new TH1F("h_na","",11,-0.5,10.5);
    
     TH1F * h_xi_no = new TH1F("h_xi_no","",1000,0,1);
     TH1F * h_xi = new TH1F("h_xi","",1000,0,1);
     TH1F * h_w_diff_no = new TH1F("h_w_diff_no","",1000,0,5000);
     TH1F * h_w_diff = new TH1F("h_w_diff","",5000,0,5000);
  
     TH1F * h_w_ratio = new TH1F("h_w_ratio","",10000,0,10);
  
     TH1F * h_rapid = new TH1F ("h_rapid","",50,-2,2); 

     TH1F * h_mcweight = new TH1F("h_mcweight","",10000,0,1);

     TH1F * h_Npu = new TH1F("h_Npu","",500,0,100);
     TH1F * h_Np_left = new TH1F("h_Np_left","",300,0,0.20);
     TH1F * h_Np_right = new TH1F("h_Np_right","",300,0,0.20);
     
     //loop on the events
  
     for (Long64_t jentry=0; jentry<nentries;jentry++) 
     {
  
  
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
        // if (Cut(ientry) < 0) continue;
  
  
        //declaration and initialization of local variables
   
        double eta1=0;
        double eta2=0;
		double eta3=0;
		double eta4=0;
        double delta_eta=0;

        double phi1=0;
        double phi2=0;
		double phi3=0;
		double phi4=0;
        double delta_phi=0;

        double e1=0;
        double e2=0;
		double e3=0;
		double e4=0;

        double e1_raw=0;
        double e2_raw=0;
		double e3_raw=0;
		double e4_raw=0;

        double px1=0;
        double px2=0;
		double px3=0;
		double px4=0;

        double py1=0;
        double py2=0;
		double py3=0;
		double py4=0;

        double pz1=0;
        double pz2=0;
		double pz3=0;
		double pz4=0;


        double pt1=0;
        double pt2=0;
		double pt3=0;
		double pt4=0;

        double pt_ratio=0;

        double w=0;

        double xi1=0;//Proton's fractional energy loss
        double xi2=0;

        double w_diff=0;

        int a_index=0;//Counts particle of interest
		int e_index=0; //Counts e+e-
		int mu_index=0; //Counts mu+mu-

        double norm_f_conv=norm_f; 
        double norm_f_eff=1;
  

		double px_system1=0;
		double px_system2=0;
		double py_system1=0;
		double py_system2=0;
		double pz_system1=0;
		double pz_system2=0;
		double e_system1=0;
		double e_system2=0;

		double pt_system1=0;
		double pt_system2=0;
		double eta_system1=0;
		double eta_system2=0;
		double phi_system1=0;
		double phi_system2=0;

        double rapid=0; //Rapidity of system 
 
        vector<PseudoJet> particles;
        vector<double> xi1_vect;
        vector<double> xi2_vect;

        //loop on the particles
  
        if(ngen>dimGen) cout << "WARNING: EVENT CUT" << endl;
  
        for(int i=0; i<ngen;i++)
        {
		  if (fabs(id[i])==16)
		  {// Id = 16 tau neutrino, we reject any tau production from Z decay for this analysis.
            cout << "WARNING: EVENT " << event_count << 
            " WITH NU_TAU and TAU" << endl;
            pt1=-999;
            pt2=-999;
			pt3=-999;
			pt4=-999;
          i=ngen;
		  e_index=5;
		  mu_index=5;
		  event_count++;
		  }
  
          if(fabs(id[i])== 11 && e_index==0)
          {
			if(!(e[i] > 4000) && CalcPt(e[i],CalcEta(pz[i],e[i]))>10){
            eta1 = CalcEta(pz[i],e[i]);
            phi1 = CalcPhi(px[i],py[i]);
            e1 = CalcE(e[i]);
			px1=px[i];
			py1=py[i];
            pt1 = CalcPt(e1,eta1);      


            e1_raw=e[i]; 
            pz1 = pz[i];

            e_index++;
			cout << "EVENT:" << event_count << "WITH 1 electron" << endl;
			}
          }

          if(fabs(id[i])== 11 && e_index==1)
          {
			if(!(e[i] > 4000) && CalcPt(e[i],CalcEta(pz[i],e[i]))>10){
            eta2 = CalcEta(pz[i],e[i]);
            phi2 = CalcPhi(px[i],py[i]);
            e2 = CalcE(e[i]);
			px2=px[i];
			py2=py[i];
            pt2 = CalcPt(e2,eta2);      

            e2_raw=e[i]; 
            pz2 = pz[i];
			cout << "EVENT:" << event_count << "WITH 2 electrons" << endl;
            e_index++;}
          }

          if(fabs(id[i])==13 && mu_index==0)
          {
			if(!(e[i] > 4000) && CalcPt(e[i],CalcEta(pz[i],e[i]))>10){
            eta3 = CalcEta(pz[i],e[i]);
            phi3 = CalcPhi(px[i],py[i]);
            e3 = CalcE(e[i]);
			px3=px[i];
			py3=py[i];
            pt3 = CalcPt(e1,eta1);      

            e3_raw=e[i]; 
            pz3 = pz[i];
			cout << "EVENT:" << event_count << "WITH 1 muon" << endl;
            mu_index++;

			}
          }

          if(fabs(id[i])==13 && mu_index==1)
          {
			if(!(e[i] > 4000) && CalcPt(e[i],CalcEta(pz[i],e[i]))>10){
            eta4 = CalcEta(pz[i],e[i]);
            phi4 = CalcPhi(px[i],py[i]);
            e4 = CalcE(e[i]);
			px4=px[i];
			py4=py[i];
            pt4 = CalcPt(e1,eta1);      
			cout << "EVENT:" << event_count << "WITH 2 muons" << endl;
            e4_raw=e[i]; 
            pz4 = pz[i];

            mu_index++;
			}
          }


          //Jet case - Excluding neutrinos and muons from the cluster 
          //         - and high eta
          else if(part_id==999 && fabs(id[i]) != 12 && fabs(id[i]) != 13 && 
                  fabs(id[i]) != 14 && fabs(id[i]) != 16 && 
                  fabs(id[i]) != 18 && CalcEta(pz[i],e[i])<2.8) //Central jet?
          {
            particles.push_back(PseudoJet(px[i],py[i],pz[i],e[i])); //These are then added to the Jet object
          } 
  
        }  // end loop on particles
  
              
        //Diffractive mass calculation - no pile-up case
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
        else
        {//pu_flag=1 pile-up case.
          //Prob to get s/d tag in AFP for 1 MB event
          TFile * input = new TFile("outPileup.root","READ");
          TH1F * h1 = (TH1F*)input->Get("hp");
          double p_s = h1->GetBinContent(2); // probability to get a single tag.
          double p_d = h1->GetBinContent(4); //probability to get a double tag from one MB event

          //pile-up draw
          double puEvents = gRandom->Poisson(50);
          h_Npu->Fill(puEvents);

          double ximin=0.015; //Standard TOTEM/AFS acceptance
          double ximax=0.15;
  
          for(int k=0;k<puEvents;k++)
          {//Draw for n/s/t tag
            double tag = gRandom->Uniform(1); //Random variable between [0,1]

            if(tag <= p_d)
            { //Inverse transform sampling for f(xi)~1/xi
              xi1_vect.push_back(ximin *
                                 pow( ximax/ximin, gRandom->Uniform(1))*
                                 gRandom->Gaus(1,0.02));
              xi2_vect.push_back(ximin * 
                                 pow( ximax/ximin, gRandom->Uniform(1))*
                                 gRandom->Gaus(1,0.02));
            }
            else if(tag>p_d && tag<=p_d+p_s)
            {
              xi1_vect.push_back(ximin *
                                 pow( ximax/ximin, gRandom->Uniform(1))* 
                                 gRandom->Gaus(1,0.02));
            }
            else if(tag>p_d+p_s && tag<=p_d+2*p_s)
            {
              xi2_vect.push_back(ximin *
                                 pow( ximax/ximin, gRandom->Uniform(1))*
                                 gRandom->Gaus(1,0.02));
            } 

          }//End for, Pileup case.
 
          h_Np_left->Fill(xi1_vect.size()); //Count number of pile-up protons with xi1
          h_Np_right->Fill(xi2_vect.size()); //Count number of pile-up protons with xi2
          input->Close();
        }  // end pile-up simulation
 
  
        //-------E+E-,MU+MU- case-------//
        if(part_id!=999)
        {//Warning and rejection if less than 4 leptons
          if(e_index+mu_index<4)
          {
            cout << "WARNING: EVENT " << event_count << 
            " WITH LESS THAN 4 CHARGED LEPTONS" << endl;
            pt1=-999;
            pt2=-999;
			pt3=-999;
			pt4=-999;
          }

          //rejecting bad events from low mass DY prod 
          if(high_flag==1 && pu_flag==1 && part_id==11 && w>400)
          {
            pt1=-999;
            pt2=-999;
			pt3=-999;
			pt4=-999;
          }

          //rejecting bad events low pt ND diphoton production
          if(high_flag==1 && pu_flag==1 && pt1>200 && pt2>200 && pt3>200 && pt4>200)
          {
            pt1=-999;
            pt2=-999;
			pt3=-999;
			pt4=-999;
          }

//        System 1 (Two electrons)
          eta_system1=CalcEta(pz1+pz2, e1_raw+e2_raw);
		  pt_system1=CalcPt(e1_raw+e2_raw, eta1+eta2);
		  phi_system1=CalcPhi(px1+px2,py1+py2);
		  e_system1=e1_raw+e2_raw;
		  px_system1=px1+px2;
		  py_system1=py1+py2;
		  pz_system1=pz1+pz2;

//        System 2 (Two muons)
          eta_system2=CalcEta(pz3+pz4, e3_raw+e4_raw);
		  pt_system2=CalcPt(e3_raw+e4_raw, eta3+eta4);
		  phi_system2=CalcPhi(px3+px4,py3+py4);
		  e_system2=e3_raw+e4_raw;
		  px_system2=px3+px4;
		  py_system2=py3+py4;
		  pz_system2=pz3+pz4;


          //Mass of the four-body system, invariant mass \approx 0 assumed for ee, mumu.
/*          w = CalcWsq_X(pt1, pt2, eta1, eta2, phi1, phi2)+CalcWsq_X(pt3, pt4, eta3, eta4, phi3, phi4)-2*LorenzProd(e_system1,e_system2,px_system1,px_system2,py_system1,py_system1,pz_system1,pz_system2);
*/

		  w=CalcWsq_X(pt_system1, pt_system2, eta_system1, eta_system2, phi_system1, phi_system2);
		  w=sqrt(w);
		  if (e_index+mu_index==4){
		  cout << w << endl;
		  cout <<pt_system1 <<endl;
		  }
          //pT before any selection
//          h_pt1_no->Fill(pt_system1,norm_f);
//          h_pt2_no->Fill(pt_system2,norm_f);
/*  
          //pt ratio 

		  if (pt_ratio!=0){
		  cout << pt_ratio << endl;
		  }

*/
          pt_ratio = pt_system2/pt_system1;
          //Delta_eta and Delta_phi

          delta_eta = fabs(eta_system1-eta_system2);
//          if(fabs(phi1-phi2)<3.14159265359) delta_phi = fabs(phi1-phi2);
//          else delta_phi = 2*3.14159265359 - fabs(phi1-phi2);

          //rapidity of the system
          rapid = CalcEta(pz1+pz2+pz3+pz4,e1_raw+e2_raw+e3_raw+e4_raw);          

          //for pile-up, selection of the xi1/xi2 which are the closest to
          //the diphoton mass/rap
          if(pu_flag==1)
          { double w_ratio_test = 999;
            double rapid_test = 999;
            if(xi1_vect.size()>0 && xi2_vect.size()>0)
            {
              for(int l=0;l<xi1_vect.size();l++)
              {   
                for(int m=0;m<xi2_vect.size();m++)
                {
                  double test1=fabs(1-sqrt(xi1_vect[l]*xi2_vect[m])*13000/w);
                  double test2=fabs(rapid-0.5*log(xi1_vect[l]/xi2_vect[m]));
                  //Selecting the good protons on mass/rapidity criteria
                  if(pow(test1,2)+pow(test2,2) < pow(w_ratio_test,2)+pow(rapid_test,2))
                  {
                    xi1=xi1_vect[l];
                    xi2=xi2_vect[m];
                    w_diff = sqrt(xi1*xi2)*13000;
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
          if(pt_system1>10 && pt_system2>10 /* && w>100 && fabs(eta1)<2.37 && 
             fabs(eta2)<2.37 && max(xi1,xi2)<0.15 && 
             min(xi1,xi2)>0.015  && 
             pt_ratio>0.85 && w_diff/w<1+mass_window && 
             w_diff/w>1-mass_window && fabs(rapid-0.5*log(xi1/xi2))<0.10/*&& fabs(3.14159265359-delta_phi)<0.01 */)
          {            
             //converted photon factor as an eta function
             //we request 1 converted photon
             //for electrons, already converted ...
             /*if(part_id==22)
             {*/
               if(fabs(eta1)<1.4 && fabs(eta2)<1.4)
               { 
                 norm_f_conv *= 0.85*0.85; 
               }
               else if(fabs(eta1)>=1.4 && fabs(eta2)>=1.4)
               {
                 norm_f_conv *= 0.7*0.7; 
               }
               else norm_f_conv *= 0.7*0.85;
             //}
   
             //efficiency factor for real photons-electrons
             //norm_f_eff=norm_f_conv*(0.76-1.98*exp(-pt1/16.1))*
             //                       (0.76-1.98*exp(-pt2/16.1));
    
			   norm_f_eff=norm_f_conv; 
              //Filling histos
             h_eta->Fill(eta_system1,norm_f_eff);
             h_eta->Fill(eta_system2,norm_f_eff);
             h_delta_eta->Fill(delta_eta,norm_f_eff);
   
//             h_phi->Fill(phi1,norm_f_eff);
//             h_phi->Fill(phi2,norm_f_eff);
//             h_delta_phi->Fill(3.14159265359-delta_phi,norm_f_eff);
   
             h_pt1->Fill(pt_system1,norm_f_eff);
             h_pt2->Fill(pt_system2,norm_f_eff);
             h_pt_ratio->Fill(1-pt_ratio,norm_f_eff);
   
             h_w->Fill(w,norm_f_eff);
//             h_na->Fill(a_index,norm_f_eff);
   
             h_xi->Fill(xi1,norm_f_eff);
             h_xi->Fill(xi2,norm_f_eff);
             h_w_diff->Fill(w_diff,norm_f_eff);
   
             h_w_ratio->Fill(w_diff/w,norm_f_eff);
  
             h_rapid->Fill(rapid-0.5*log(xi1/xi2),norm_f_eff);
            
             h_mcweight->Fill(norm_f);
 
           } //endif 2nd part of selection
         } //endif parti_id!=999
      
  
         //-------JET CASE-------//
         else //parti_id=999
         {
           // choose a jet definition
           double R = 0.4;
           JetDefinition jet_def(antikt_algorithm, R); //Check FastJet for more algorithms.
   
           // run the clustering, extract the jets > pTmin = 50 GeV 
           // at the moment
           ClusterSequence cs(particles, jet_def);
           vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets(50));
   
           a_index = jets.size();
   
           //select the two leading jets
           if(a_index>1)
           {
   
             e1 = jets[0].e();
             eta1 = jets[0].eta();
             phi1 = jets[0].phi();
   
             e2 = jets[1].e();       
             eta2 = jets[1].eta();
             phi2 = jets[1].phi();
        
             //retrieving fake gamma E
             //produce random numbers according to a normal distrib 
             e1 = e1*gRandom->Gaus(0.75,0.13);  
             e2 = e2*gRandom->Gaus(0.75,0.13); 
 
             eta1 = gRandom->Gaus(eta1,0.001); 
             eta2 = gRandom->Gaus(eta2,0.001); 
             phi1 = gRandom->Gaus(phi1,0.001); 
             phi2 = gRandom->Gaus(phi2,0.001); 
 
             pt1 = CalcPt(e1,eta1);
             pt2 = CalcPt(e2,eta2);

//Probably need to add some smearing to pt1 and pt2?... Or already included in eta, e1?...
 
             //Jet pTs
             h_pt1_no->Fill(pt1,norm_f);//No pile-up considered
             h_pt2_no->Fill(pt2,norm_f);//No pile-up considered
 
             pz1 = jets[0].pz();
             pz2 = jets[1].pz();
 
             e1_raw = jets[0].e();
             e2_raw = jets[1].e();

             w = sqrt(2*pt1*pt2*(cosh(eta1-eta2)-cos(phi1-phi2)));
 
             //pt_ratio
             pt_ratio = pt2/pt1;
 
             //delta_eta and delta_phi
             delta_eta = fabs(eta1-eta2);
     
             if(fabs(phi1-phi2)<3.14159) delta_phi = fabs(phi1-phi2);
             else delta_phi = 2*3.14159 - fabs(phi1-phi2);
 
             //rapidity of the jet system
             rapid = CalcEta(pz1+pz2,e1_raw+e2_raw);
    
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
                     double test1=fabs(1-sqrt(xi1_vect[l]*xi2_vect[m])*
                                       13000/w);
                     double test2=fabs(rapid-0.5*log(xi1_vect[l]/xi2_vect[m]));
   
                     //selecting the good protons on mass/rap criteria
                     if(pow(test1,2)+pow(test2,2) < 
                        pow(w_ratio_test,2)+pow(rapid_test,2))
                     {
                       xi1=xi1_vect[l];
                       xi2=xi2_vect[m];
                       w_diff = sqrt(xi1*xi2)*13000;
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
   

             //Online Selection: pT, eta, sqrt(xi_1xi_2s) 
             if(pt1>200 && pt2>100 && w>400 && fabs(eta1)<2.37 && 
                fabs(eta2)<2.37 && max(xi1,xi2)<0.15 && min(xi1,xi2)>0.015 && 
                fabs(3.14159265359-delta_phi)<0.01 && pt_ratio>0.95 &&
                 fabs(rapid-0.5*log(xi1/xi2))<0.03 /*&& w_diff/w>1-mass_window &&
                w_diff/w<1+mass_window */)
             {
 
  
               //converted photon factor as an eta function, 
               //we request 1 converted photon
              /* if(fabs(eta1)<1.4 && fabs(eta2)<1.4)
               {
                 norm_f_conv *= (1-(0.85*0.85)); 
               }
               else if(fabs(eta1)>=1.4 && fabs(eta2)>=1.4)
               {
                 norm_f_conv *= (1-(0.7*0.7)); 
               }
               else norm_f_conv *= (1-(0.7*0.85));
				*/
 
   
               //fake photons reconstruction efficiency until 200GeV
               //norm_f_eff=norm_f_conv*(0.0093*exp(-min(pt1,200.0)/27.5)*
               //                        0.0093*exp(-min(pt2,200.0)/27.5));
   
               //drawing
               h_eta->Fill(eta1,norm_f_eff);
               h_eta->Fill(eta2,norm_f_eff);
               h_delta_eta->Fill(delta_eta,norm_f_eff);
   
               h_phi->Fill(phi1,norm_f_eff);
               h_phi->Fill(phi2,norm_f_eff);
//               h_delta_phi->Fill(3.14159265359-delta_phi,norm_f_eff);
   
               h_pt1->Fill(pt1,norm_f_eff);
               h_mcweight->Fill(0.0093*exp(-min(pt1,200.0)/27.5)*
                                0.0093*exp(-min(pt2,200.0)/27.5)*norm_f);
               h_pt2->Fill(pt2,norm_f_eff);
               h_pt_ratio->Fill(1-pt_ratio,norm_f_eff);
   
               h_w->Fill(w,norm_f_eff);
               h_na->Fill(a_index,norm_f_eff);
    
               h_xi->Fill(xi1,norm_f_eff);
               h_xi->Fill(xi2,norm_f_eff);
               h_w_diff->Fill(w_diff,norm_f_eff);
               h_w_ratio->Fill(w_diff/w,norm_f_eff);
   
               h_rapid->Fill(rapid-0.5*log(xi1/xi2),norm_f_eff);
  
             } // End selection
     
     
           } // End if > 2 jets 
     
         } // End jet case
    
         event_count++;
       }
  
       //Compute MC error
  
       cout << "NORM MC TO 300 FB-1 (Wrong for jets): " << norm_f << endl;
	     
       out->cd();
       out->Write();
  
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
