/**********************************************************************************
 * Project   : TMVA - a Root-integrated toolkit for multivariate data analysis    *
 * Package   : TMVA                                                               *
 * Executable: TMVAMulticlass                                                     *
 *                                                                                *
 * This macro provides a simple example for the training and testing of the TMVA  *
 * multiclass classification                                                      *
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

#include "TMVAMultiClassGui.C"

#ifndef __CINT__
#include "TMVA/Tools.h"
#include "TMVA/Factory.h"
#endif

using namespace TMVA;

typedef std::map<std::string, int> StringIntMap;

//------------------------------------------------------------------------------------------------------------------------------------------

class Parameters
{
public:
    Parameters(const TString &analysisName, const TString &energyName, const TString &tfileName, const TString ttreeName, const float nEvent);

    TString     m_analysisName;
    TString     m_energyName;
    TString     m_tfileName;
    TString     m_ttreeName;
    TChain     *m_pTChain;
    float       m_weight;
};

//------------------------------------------------------------------------------------------------------------------------------------------

Parameters::Parameters(const TString &analysisName, const TString &energyName, const TString &tfileName, const TString ttreeName, const float nEvent):
m_analysisName(analysisName),
m_energyName(energyName),
m_tfileName(tfileName),
m_ttreeName(ttreeName)
{
    m_pTChain = new TChain(m_ttreeName.Data());
    m_pTChain->Add(m_tfileName.Data());
    m_weight = nEvent / static_cast<float>(m_pTChain->GetEntries());
}

//------------------------------------------------------------------------------------------------------------------------------------------

StringIntMap InitialiseMVAMethods(int argc, char** argv)
{
    //---------------------------------------------------------------
    // default MVA methods to be trained + tested
    StringIntMap Use;
    Use["MLP"]             = 1;
    Use["BDTG"]            = 1;
    Use["FDA_GA"]          = 0;
    Use["PDEFoam"]         = 0;
    //---------------------------------------------------------------
    if (argc>1)
    {
        for (std::map<std::string,int>::iterator it = Use.begin();
                it != Use.end(); it++)
        {
            it->second = 0;
        }
    }
    for (int i=1; i<argc; i++)
    {
        std::string regMethod(argv[i]);
        if (Use.find(regMethod) == Use.end())
        {
            std::cout << "Method " << regMethod << " not known in TMVA under this name. Please try one of:" << std::endl;
            for (StringIntMap::const_iterator it = Use.begin(); it != Use.end(); std::cout << it++->first << " ");
            std::cout << std::endl;
            return StringIntMap();
        }
        else
        {
            Use[regMethod] = true;
        }
        
    }

    return Use;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void BookMethod(TMVA::Factory *&pFactory, const StringIntMap &mvaMethodVec)
{
    if (mvaMethodVec.at("BDTG")) // gradient boosted decision trees
        pFactory->BookMethod( TMVA::Types::kBDT, "BDTG", "!H:!V:NTrees=1000:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.50:nCuts=20:MaxDepth=2");
    if (mvaMethodVec.at("MLP")) // neural network
        pFactory->BookMethod( TMVA::Types::kMLP, "MLP", "!H:!V:NeuronType=tanh:NCycles=300:HiddenLayers=N+5,5:TestRate=5:EstimatorType=MSE");
    if (mvaMethodVec.at("FDA_GA")) // functional discriminant with GA minimizer
        pFactory->BookMethod( TMVA::Types::kFDA, "FDA_GA", "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:PopSize=300:Cycles=3:Steps=20:Trim=True:SaveBestGen=1" );
    if (mvaMethodVec.at("PDEFoam")) // PDE-Foam approach
        pFactory->BookMethod( TMVA::Types::kPDEFoam, "PDEFoam", "!H:!V:TailCut=0.001:VolFrac=0.0666:nActiveCells=500:nSampl=2000:nBin=5:Nmin=100:Kernel=None:Compress=T" );
}

//------------------------------------------------------------------------------------------------------------------------------------------

int main(int argc, char** argv)
{
    Parameters parameters("TauTauAnalysis", "100GeV", "/r06/lc/xu/TautauAnalysis/PandoraDefault/100GeV/rootCustom/tauAnalysisTemplate_10.root", "sel", 52869);
    // ATTN assume eventtype starts at 1 ends in 1 + nEvtType
    const int nEvtType(7);
    TMVA::Tools::Instance();

    const StringIntMap &mvaMethodVec(InitialiseMVAMethods(argc, argv));
    if (mvaMethodVec.empty())
        return 1;

    std::cout << std::endl;
    std::cout << "==> Start TMVAMulticlass" << std::endl;

    // Create a new root output file.

    const TString outfileNamePrefix(parameters.m_analysisName + "_" + parameters.m_energyName + "_TMVAMulticlass");
    
    const TString outfileName(outfileNamePrefix + ".root");
    TFile* pOutputFile = TFile::Open(outfileName, "RECREATE");
    TMVA::Factory *pFactory = new TMVA::Factory(outfileNamePrefix, pOutputFile,
        "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=multiclass" );

    //pFactory->AddVariable( "mVis", "visible mass", "GeV", 'F' );
    //pFactory->AddVariable( "nCharge", "number of charged PFOs", "", 'I' );
    pFactory->AddVariable( "thrustPrinciple", "thrust", "", 'F' );
    pFactory->AddVariable( "EEHCalRatioFixed:= EEHCalRatio > 1. ? 0. : EEHCalRatio", "E ECal to All ratio", "", 'F' );
    pFactory->AddVariable( "nPfo", "nPfo", "", 'I' );
    pFactory->AddVariable( "mVis", "mVis", "", 'I' );
    pFactory->AddVariable( "nCharge", "nCharge", "", 'I' );
    pFactory->AddVariable( "nNeutral", "nNeutral", "", 'I' );
    pFactory->AddVariable( "nMuon", "nMuon", "", 'I' );
    pFactory->AddVariable( "eMuon", "eMuon", "GeV", 'F' );
    pFactory->AddVariable( "nElectron", "nElectron", "", 'I' );
    pFactory->AddVariable( "eElectron", "eElectron", "GeV", 'F' );
    pFactory->AddVariable( "nPhoton", "nPhoton", "", 'I' );
    pFactory->AddVariable( "mPhoton", "mPhoton", "GeV", 'F' );
    pFactory->AddVariable( "ePhoton", "ePhoton", "GeV", 'F' );
    pFactory->AddVariable( "nPionCharge", "nPionCharge", "", 'I' );
    pFactory->AddVariable( "mPionCharge", "mPionCharge", "GeV", 'F' );
    pFactory->AddVariable( "ePionCharge", "ePionCharge", "GeV", 'F' );
    pFactory->AddVariable( "logChi2RhoFitFix:= logChi2RhoFit < -1000 ? -TMath::Log(-logChi2RhoFit) : logChi2RhoFit", "-logChi2 Rho Fit", "", 'F' );
    pFactory->AddVariable( "mPionRhoFit", "mPion Rho Fit", "GeV", 'F' );
    pFactory->AddVariable( "mRhoRhoFit", "mRho Fit", "GeV", 'F' );
    pFactory->AddVariable( "cosStarRhoFit", "cos* Rho Fit", "GeV", 'F' );
    pFactory->AddVariable( "logChi2A1FitFix:= logChi2A1Fit < -1000 ? -TMath::Log(-logChi2A1Fit) : logChi2A1Fit", "-logChi2 A1 Fit", "", 'F' );
    pFactory->AddVariable( "mPionLhsA1Fit", "mPion lhs A1 Fit", "GeV", 'F' );
    pFactory->AddVariable( "mPionRhsA1Fit", "mPion Rhs A1 Fit", "GeV", 'F' );
    pFactory->AddVariable( "mA1A1Fit", "mA1 Fit", "GeV", 'F' );
    pFactory->AddVariable( "cosStarLhsA1Fit", "cos* lhs A1 Fit", "GeV", 'F' );
    pFactory->AddVariable( "cosStarRhsA1Fit", "cos* rhs A1 Fit", "GeV", 'F' );

    pFactory->AddSpectator( "nEvt", 'I' );
    pFactory->AddSpectator( "eventType", 'I' );
    pFactory->AddSpectator( "photonEC", "Photon early conversion", "", 'I' );
    pFactory->AddSpectator( "mcCloseToZ", "MC close to beam pipe", "", 'I' );

    const TCut commonCut("photonEC < 1 && mcCloseToZ < 1");
    for (int iter = 1, iterEnd = 1 + nEvtType; iter != iterEnd; ++iter)
    {
        pFactory->AddTree(parameters.m_pTChain,  TString::Format("eventType%d", iter), parameters.m_weight,
            commonCut + TCut("eventType == " + TString::Format("%d", iter)));
    }

    pFactory->PrepareTrainingAndTestTree("", "SplitMode=Random:NormMode=NumEvents:!V" );

    BookMethod(pFactory, mvaMethodVec);
    
    // Train MVAs using the set of training events
    pFactory->TrainAllMethods();

    // ---- Evaluate all MVAs using the set of test events
    pFactory->TestAllMethods();

    // ----- Evaluate and compare performance of all configured MVAs
    pFactory->EvaluateAllMethods();

    // --------------------------------------------------------------

    // Save the output
    pOutputFile->Close();

    std::cout << "==> Wrote root file: " << pOutputFile->GetName() << std::endl;
    std::cout << "==> TMVAClassification is done!" << std::endl
    << std::endl
    << "==> To view the results, launch the GUI: \"root -l ./TMVAMultiClassGui.C\"" << std::endl
    << std::endl;
    delete pFactory;

}


