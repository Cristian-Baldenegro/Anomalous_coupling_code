//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Oct  2 16:04:16 2013 by ROOT version 5.34/02
// from TTree h777/ntuple
// found on file: aa_anom_a1a_ff_1.root
//////////////////////////////////////////////////////////

#ifndef aaAnom_tree_h
#define aaAnom_tree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class aaAnom_tree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           ngen;
   Float_t         px[5000];   //[ngen]
   Float_t         py[5000];   //[ngen]
   Float_t         pz[5000];   //[ngen]
   Float_t         e[5000];   //[ngen]
   Float_t         rm[5000];   //[ngen]
   Int_t           id[5000];   //[ngen]

   // List of branches
   TBranch        *b_ngen;   //!
   TBranch        *b_px;   //!
   TBranch        *b_py;   //!
   TBranch        *b_pz;   //!
   TBranch        *b_e;   //!
   TBranch        *b_rm;   //!
   TBranch        *b_id;   //!

   aaAnom_tree(TTree *tree=0);
   virtual ~aaAnom_tree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef aaAnom_tree_cxx
aaAnom_tree::aaAnom_tree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("aa_anom_a1a_ff_1.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("aa_anom_a1a_ff_1.root");
      }
      f->GetObject("h777",tree);

   }
   Init(tree);
}

aaAnom_tree::~aaAnom_tree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t aaAnom_tree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t aaAnom_tree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void aaAnom_tree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("ngen", &ngen, &b_ngen);
   fChain->SetBranchAddress("px", px, &b_px);
   fChain->SetBranchAddress("py", py, &b_py);
   fChain->SetBranchAddress("pz", pz, &b_pz);
   fChain->SetBranchAddress("e", e, &b_e);
   fChain->SetBranchAddress("rm", rm, &b_rm);
   fChain->SetBranchAddress("id", id, &b_id);
   Notify();
}

Bool_t aaAnom_tree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void aaAnom_tree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t aaAnom_tree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef aaAnom_tree_cxx
