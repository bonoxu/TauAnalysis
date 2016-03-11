/**********************************************************************************
 * Project   : TMVA - a Root-integrated toolkit for multivariate data analysis    *
 * Package   : TMVA                                                               *
 * Exectuable: TMVAMulticlassApplication                                          *
 *                                                                                *
 * This macro provides a simple example on how to use the trained multiclass      *
 * classifiers within an analysis module                                          *
 **********************************************************************************/

#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TH1F.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/Config.h"

using namespace TMVA;
//------------------------------------------------------------------------------------------------------------------------------------------

class Parameters
{
public:
    Parameters(const TString &analysisName, const TString &energyName, const TString &tfileName, const TString ttreeName, const TString &trainedAnalysisName, 
        const TString &trainedEnergyName);

    TString     m_analysisName;
    TString     m_trainedEnergyName;
    TString     m_trainedAnalysisName;
    TString     m_energyName;
    TString     m_tfileName;
    TString     m_ttreeName;
    TChain     *m_pTChain;
};

//------------------------------------------------------------------------------------------------------------------------------------------

Parameters::Parameters(const TString &analysisName, const TString &energyName, const TString &tfileName, const TString ttreeName, const TString &trainedAnalysisName, 
        const TString &trainedEnergyName):
m_analysisName(analysisName),
m_energyName(energyName),
m_tfileName(tfileName),
m_ttreeName(ttreeName),
m_trainedEnergyName(trainedEnergyName),
m_trainedAnalysisName(trainedAnalysisName)
{
    m_pTChain = new TChain(m_ttreeName.Data());
    m_pTChain->Add(m_tfileName.Data());
}

//------------------------------------------------------------------------------------------------------------------------------------------

int main(int argc, char** argv )
{
   Parameters parameters("TauTauAnalysisBono", "200GeV", "/r06/lc/xu/TautauAnalysis/PandoraBono/200GeV/rootCustom/tauAnalysisTemplate_10.root", "sel", "TauTauAnalysisBono", "100GeV");
   TMVA::Tools::Instance();
   
   //---------------------------------------------------------------
   // default MVA methods to be trained + tested
   std::map<std::string,int> Use;
   Use["MLP"]             = 1;
   Use["BDTG"]            = 1;
   Use["FDA_GA"]          = 0;
   Use["PDEFoam"]         = 0;
   //---------------------------------------------------------------
  
   std::cout << std::endl;
   std::cout << "==> Start TMVAMulticlassApplication" << std::endl; 

   if (argc>1) {
      for (std::map<std::string,int>::iterator it = Use.begin();
           it != Use.end(); it++) {
         it->second = 0;
      }
   }
   for (int i=1; i<argc; i++) {
      std::string regMethod(argv[i]);
      if (Use.find(regMethod) == Use.end()) {
         std::cout << "Method " << regMethod << " not known in TMVA under this name. Please try one of:" << std::endl;
         for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) std::cout << it->first << " ";
         std::cout << std::endl;
         return 1;
      }
      Use[regMethod] = kTRUE;
   }
   
   // book output histograms
   UInt_t nbin = 100;
   TH1F *histMLP_signal(0), *histBDTG_signal(0), *histFDAGA_signal(0), *histPDEFoam_signal(0);
   if (Use["MLP"])    
      histMLP_signal    = new TH1F( "MVA_MLP_signal",    "MVA_MLP_signal",    nbin, 0., 1.1 );
   if (Use["BDTG"])
      histBDTG_signal  = new TH1F( "MVA_BDTG_signal",   "MVA_BDTG_signal",   nbin, 0., 1.1 );
   if (Use["FDA_GA"])
      histFDAGA_signal = new TH1F( "MVA_FDA_GA_signal", "MVA_FDA_GA_signal", nbin, 0., 1.1 );
   if (Use["PDEFoam"])
      histPDEFoam_signal = new TH1F( "MVA_PDEFoam_signal", "MVA_PDEFoam_signal", nbin, 0., 1.1 );

   // create the Reader object
   TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );    

   // create a set of variables and declare them to the reader
   // - the variable names must corresponds in name and type to 
   // those given in the weight file(s) that you use
   Float_t float1, float2, float3, float4, float5, float6, float7, float8, float9, float10, float11, float12, float13, float14, float15, float16
   , float17, float18, float19;
   Float_t int1, int2, int3, int4, int5, int6, int7, int8, int9, int10, int11, int12;

    reader->AddVariable( "thrustPrinciple", &float1);
    reader->AddVariable( "EEHCalRatio > 1. ? 0. : EEHCalRatio", &float2);
    reader->AddVariable( "nPfo", &int1);
    reader->AddVariable( "mVis", &float3);
    reader->AddVariable( "nCharge", &int4);
    reader->AddVariable( "nNeutral", &int5);
    reader->AddVariable( "nMuon", &int6);
    reader->AddVariable( "eMuon", &float4);
    reader->AddVariable( "nElectron", &int7);
    reader->AddVariable( "eElectron", &float5);
    reader->AddVariable( "nPhoton", &int8);
    reader->AddVariable( "mPhoton", &float6);
    reader->AddVariable( "ePhoton", &float7);
    reader->AddVariable( "nPionCharge", &int9);
    reader->AddVariable( "mPionCharge", &float8);
    reader->AddVariable( "ePionCharge", &float9);
    reader->AddVariable( "logChi2RhoFit < -1000 ? -TMath::Log(-logChi2RhoFit) : logChi2RhoFit", &float10);
    reader->AddVariable( "mPionRhoFit", &float11);
    reader->AddVariable( "mRhoRhoFit", &float12);
    reader->AddVariable( "cosStarRhoFit", &float13);
    reader->AddVariable( "logChi2A1Fit < -1000 ? -TMath::Log(-logChi2A1Fit) : logChi2A1Fit", &float14);
    reader->AddVariable( "mPionLhsA1Fit", &float15);
    reader->AddVariable( "mPionRhsA1Fit", &float16);
    reader->AddVariable( "mA1A1Fit", &float17);
    reader->AddVariable( "cosStarLhsA1Fit", &float18);
    reader->AddVariable( "cosStarRhsA1Fit", &float19);

    reader->AddSpectator( "nEvt", &int9 );
    reader->AddSpectator( "eventType", &int10 );
    reader->AddSpectator( "photonEC", &int11);
    reader->AddSpectator( "mcCloseToZ", &int12);
   // book the MVA methods
   const TString dir    = "weights/";
   const TString inputfileNamePrefix(parameters.m_trainedAnalysisName + "_" + parameters.m_trainedEnergyName + "_TMVAMulticlass");
   for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) {
      if (it->second) {
         TString methodName = it->first + " method";
         TString weightfile = dir + inputfileNamePrefix + "_" + TString(it->first) + ".weights.xml";
         reader->BookMVA( methodName, weightfile ); 
      }
   }

   std::cout << "--- TMVAMulticlassApp : Using input file: " << parameters.m_tfileName << std::endl;
   
   // prepare the tree
   // - here the variable names have to corresponds to your tree
   // - you can use the same variables as above which is slightly faster,
   //   but of course you can use different ones and copy the values inside the event loop
  
   std::cout << "--- Select signal sample" << std::endl;

    parameters.m_pTChain->SetBranchAddress( "thrustPrinciple", &float1);
    parameters.m_pTChain->SetBranchAddress( "EEHCalRatio", &float2);
    parameters.m_pTChain->SetBranchAddress( "nPfo", &int1);
    parameters.m_pTChain->SetBranchAddress( "mVis", &float3);
    parameters.m_pTChain->SetBranchAddress( "nCharge", &int4);
    parameters.m_pTChain->SetBranchAddress( "nNeutral", &int5);
    parameters.m_pTChain->SetBranchAddress( "nMuon", &int6);
    parameters.m_pTChain->SetBranchAddress( "eMuon", &float4);
    parameters.m_pTChain->SetBranchAddress( "nElectron", &int7);
    parameters.m_pTChain->SetBranchAddress( "eElectron", &float5);
    parameters.m_pTChain->SetBranchAddress( "nPhoton", &int8);
    parameters.m_pTChain->SetBranchAddress( "mPhoton", &float6);
    parameters.m_pTChain->SetBranchAddress( "ePhoton", &float7);
    parameters.m_pTChain->SetBranchAddress( "nPionCharge", &int9);
    parameters.m_pTChain->SetBranchAddress( "mPionCharge", &float8);
    parameters.m_pTChain->SetBranchAddress( "ePionCharge", &float9);
    parameters.m_pTChain->SetBranchAddress( "logChi2RhoFit", &float10);
    parameters.m_pTChain->SetBranchAddress( "mPionRhoFit", &float11);
    parameters.m_pTChain->SetBranchAddress( "mRhoRhoFit", &float12);
    parameters.m_pTChain->SetBranchAddress( "cosStarRhoFit", &float13);
    parameters.m_pTChain->SetBranchAddress( "logChi2A1Fit", &float14);
    parameters.m_pTChain->SetBranchAddress( "mPionLhsA1Fit", &float15);
    parameters.m_pTChain->SetBranchAddress( "mPionRhsA1Fit", &float16);
    parameters.m_pTChain->SetBranchAddress( "mA1A1Fit", &float17);
    parameters.m_pTChain->SetBranchAddress( "cosStarLhsA1Fit", &float18);
    parameters.m_pTChain->SetBranchAddress( "cosStarRhsA1Fit", &float19);

    parameters.m_pTChain->SetBranchAddress( "nEvt", &int9 );
    parameters.m_pTChain->SetBranchAddress( "eventType", &int10 );
    parameters.m_pTChain->SetBranchAddress( "photonEC", &int11);
    parameters.m_pTChain->SetBranchAddress( "mcCloseToZ", &int12);
    
   std::cout << "--- Processing: " << parameters.m_pTChain->GetEntries() << " events" << std::endl;
   TStopwatch sw;
   sw.Start();

   for (Long64_t ievt=0; ievt<parameters.m_pTChain->GetEntries();ievt++) {
      if (ievt%1000 == 0){
         std::cout << "--- ... Processing event: " << ievt << std::endl;
      }
      
      parameters.m_pTChain->GetEntry(ievt);
      if (Use["MLP"])
         histMLP_signal->Fill((reader->EvaluateMulticlass( "MLP method" ))[0]);
      if (Use["BDTG"])
         histBDTG_signal->Fill((reader->EvaluateMulticlass( "BDTG method" ))[0]);
      if (Use["FDA_GA"])
         histFDAGA_signal->Fill((reader->EvaluateMulticlass( "FDA_GA method" ))[0]);
      if (Use["PDEFoam"])
         histPDEFoam_signal->Fill((reader->EvaluateMulticlass( "PDEFoam method" ))[0]);
      
   }
   
   // get elapsed time
   sw.Stop();
   std::cout << "--- End of event loop: "; sw.Print();
   const TString outputfileNamePrefix(parameters.m_analysisName + "_" + parameters.m_energyName + "_use_" + parameters.m_trainedAnalysisName + "_" + parameters.m_trainedEnergyName);
   TFile *target  = new TFile( outputfileNamePrefix + ".root","RECREATE" );
   if (Use["MLP"])
      histMLP_signal->Write();
   if (Use["BDTG"])
      histBDTG_signal->Write(); 
   if (Use["FDA_GA"])
      histFDAGA_signal->Write();
   if (Use["PDEFoam"])
      histPDEFoam_signal->Write();

   target->Close();
   std::cout << "--- Created root file: \"" << outputfileNamePrefix <<".root\" containing the MVA output histograms" << std::endl;

   delete reader;
   
   std::cout << "==> TMVAClassificationApplication is done!" << std::endl << std::endl;
}
