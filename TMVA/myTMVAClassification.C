// @(#)root/tmva $Id$
/**********************************************************************************
 * Project   : TMVA - a ROOT-integrated toolkit for multivariate data analysis    *
 * Package   : TMVA                                                               *
 * Root Macro: TMVAClassification                                                 *
 *                                                                                *
 * This macro provides examples for the training and testing of the               *
 * TMVA classifiers.                                                              *
 *                                                                                *
 * As input data is used a toy-MC sample consisting of four Gaussian-distributed  *
 * and linearly correlated input variables.                                       *
 *                                                                                *
 * The methods to be used can be switched on and off by means of booleans, or     *
 * via the prompt command, for example:                                           *
 *                                                                                *
 *    root -l ./TMVAClassification.C\(\"Fisher,Likelihood\"\)                     *
 *                                                                                *
 * (note that the backslashes are mandatory)                                      *
 * If no method given, a default set of classifiers is used.                      *
 *                                                                                *
 * The output file "TMVA.root" can be analysed with the use of dedicated          *
 * macros (simply say: root -l <macro.C>), which can be conveniently              *
 * invoked through a GUI that will appear at the end of the run of this macro.    *
 * Launch the GUI via the command:                                                *
 *                                                                                *
 *    root -l ./TMVAGui.C                                                         *
 *                                                                                *
 **********************************************************************************/

#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include "dirent.h"

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "TMVAGui.C"

#if not defined(__CINT__) || defined(__MAKECINT__)
// needs to be included when makecint runs (ACLIC)
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#endif

typedef std::vector<TString> TStringVec;
typedef std::map<TChain*, double> TChainDoubleMap;

//-----------------------------------------------------------------------------------------------------------------------------------------

double GetWeight(const TString &fileName, const TChain* pTChain)
{
    const double integratedLumi1400GeV(1500.f); //unit fb
    double effectiveCrossSection(1.f); // unit fb^-1
    const double lumiCorrectionEGamma(0.75f),lumiCorrectionGammaE(0.75f), lumiCorrectionGammaGamma(0.64f);
    if (fileName.Contains("prod6022"))
    {
        effectiveCrossSection = 0.149; // see https://twiki.cern.ch/twiki/bin/view/CLIC/MonteCarloSamplesForTheHiggsPaper
    }
    else if (fileName.Contains("prod5527"))
    {
        effectiveCrossSection = 23.2;
    }
    else if (fileName.Contains("prod5572"))
    {
        effectiveCrossSection = 62.1;
    }
    else if (fileName.Contains("prod5594"))
    {
        effectiveCrossSection = 110.4;
    }
    else if (fileName.Contains("prod4034"))
    {
        effectiveCrossSection = 1245.1;
    }
    else if (fileName.Contains("prod3366"))
    {
        effectiveCrossSection = 287.1 * lumiCorrectionEGamma;
    }
    else if (fileName.Contains("prod3369"))
    {
        effectiveCrossSection = 1160.7 * lumiCorrectionEGamma;
    }
    else if (fileName.Contains("prod3372"))
    {
        effectiveCrossSection = 286.9  * lumiCorrectionGammaE;
    }
    else if (fileName.Contains("prod3375"))
    {
        effectiveCrossSection = 1156.3 * lumiCorrectionGammaE;
    }
    else if (fileName.Contains("prod3378"))
    {
        effectiveCrossSection = 32.6 * lumiCorrectionEGamma;
    }
    else if (fileName.Contains("prod3381"))
    {
        effectiveCrossSection = 136.9 * lumiCorrectionEGamma;
    }
    else if (fileName.Contains("prod3384"))
    {
        effectiveCrossSection = 32.6 * lumiCorrectionGammaE;
    }
    else if (fileName.Contains("prod3387"))
    {
        effectiveCrossSection = 136.4 * lumiCorrectionGammaE;
    }
    else if (fileName.Contains("prod3414"))
    {
        effectiveCrossSection = 753.0 * lumiCorrectionGammaGamma;
    }
    else if (fileName.Contains("prod3417"))
    {
        effectiveCrossSection = 4034.8 * lumiCorrectionGammaGamma;
    }
    else if (fileName.Contains("prod3420"))
    {
        effectiveCrossSection = 4018.7 * lumiCorrectionGammaGamma;
    }
    else if (fileName.Contains("prod3423"))
    {
        effectiveCrossSection = 21406.2 * lumiCorrectionGammaGamma;
    }
    else if (fileName.Contains("prod3243"))
    {
        effectiveCrossSection = 787.7;
    }
    else if (fileName.Contains("prod3246"))
    {
        effectiveCrossSection = 2725.8;
    }
    else if (fileName.Contains("prod3249"))
    {
        effectiveCrossSection = 4309.7;
    }
    else
    {
        std::cout<<"GetWeight: do not recognise the file for weight "<<fileName<<std::endl;
        return 1.f;
    }
    return (2 * effectiveCrossSection * integratedLumi1400GeV / static_cast<double>(pTChain->GetEntries()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool fileFilter(TString &fileName, const TStringVec &filters )
{

    for (unsigned int i=0; i<filters.size(); ++i) {
        //std::cout<<"fileName "<<fileName<<" filters[i] " <<filters[i]<<std::endl;
        Bool_t found = fileName.Contains(filters[i]);
        if (!found) return false;
    }
    return true;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void GetRootFileList(const TString &inputFolderName, TStringVec &filters, TStringVec &filesNames)
{
    TString folder(inputFolderName);
    TStringVec folderVec = TMVA::gTools().SplitString( folder, '*' );
    
    TString filenameFirstPart(""), folderName("");
    filenameFirstPart = folderVec[0];
    if (folderVec.size()>1){
        for (unsigned int i = 1; i<folderVec.size(); ++i){
            filters.push_back(folderVec[i]);
        }
    }
    
    std::string folderStr(std::string(folder.Data()));
    const size_t idx(folderStr.find_last_of("\\/"));
    if (std::string::npos != idx) {
        filenameFirstPart = (folderStr.substr(idx + 1));
        folderName = (folderStr.substr(0,idx));
    }
    //std::cout<<folder.substr(0,idx)<<" "<<filename<<std::endl;

    DIR *dir;
    struct dirent *ent;
    filters.push_back(filenameFirstPart);
    //cout<<folderName<<endl;
    if ((dir = opendir (folderName.Data() )) != NULL) {
        /* print all the files and directories within directory */
        while ((ent = readdir (dir)) != NULL) {
            // printf ("%s\n", ent->d_name);
            TString fileName(ent->d_name);
            if ( fileFilter(fileName,filters) ) {
                filesNames.push_back(folderName + TString("/") + ent->d_name);
            }
        }
        closedir (dir);
    } else {
        /* could not open directory */
        //perror ("");
        std::cout<<"could not open directory "<<folderName <<std::endl;
        return;
    }
    return;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void GetFileList(const TString &mySginalList, TStringVec &filters, TStringVec &filesNames)
{
    TStringVec folderVec = TMVA::gTools().SplitString( mySginalList, ',' );
    for (std::vector<TString>::const_iterator iter = folderVec.begin(), iterEnd = folderVec.end(); iter != iterEnd; ++iter)
    {
        TStringVec newfilters = filters;
        (void) GetRootFileList(*iter, newfilters, filesNames);
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void myTMVAClassification(TString myMethodList = "")
{
    // The explicit loading of the shared libTMVA is done in TMVAlogon.C, defined in .rootrc
    // if you use your private .rootrc, or run from a different directory, please copy the
    // corresponding lines from .rootrc

    // methods to be processed can be given as an argument; use format:
    //
    // mylinux~> root -l TMVAClassification.C\(\"myMethod1,myMethod2,myMethod3\"\)
    //
    // if you like to use a method via the plugin mechanism, we recommend using
    //
    // mylinux~> root -l TMVAClassification.C\(\"P_myMethod\"\)
    // (an example is given for using the BDT as plugin (see below),
    // but of course the real application is when you write your own
    // method based)

    //---------------------------------------------------------------
    // This loads the library
    TMVA::Tools::Instance();

    // Default MVA methods to be trained + tested
    std::map<std::string,int> Use;

    // --- Cut optimisation
    Use["Cuts"]            = 1;
    Use["CutsD"]           = 1;
    Use["CutsPCA"]       = 0;
    Use["CutsGA"]          = 0;
    Use["CutsSA"]          = 0;
    //
    // --- 1-dimensional likelihood ("naive Bayes estimator")
    Use["Likelihood"]      = 1;
    Use["LikelihoodD"]     = 0; // the "D" extension indicates decorrelated input variables (see option strings)
    Use["LikelihoodPCA"]   = 1; // the "PCA" extension indicates PCA-transformed input variables (see option strings)
    Use["LikelihoodKDE"]   = 0;
    Use["LikelihoodMIX"]   = 0;
    //
    // --- Mutidimensional likelihood and Nearest-Neighbour methods
    Use["PDERS"]           = 1;
    Use["PDERSD"]          = 0;
    Use["PDERSPCA"]        = 0;
    Use["PDEFoam"]         = 1;
    Use["PDEFoamBoost"]    = 0; // uses generalised MVA method boosting
    Use["KNN"]             = 1; // k-nearest neighbour method
    //
    // --- Linear Discriminant Analysis
    Use["LD"]              = 1; // Linear Discriminant identical to Fisher
    Use["Fisher"]          = 0;
    Use["FisherG"]         = 0;
    Use["BoostedFisher"]   = 0; // uses generalised MVA method boosting
    Use["HMatrix"]         = 0;
    //
    // --- Function Discriminant analysis
    Use["FDA_GA"]          = 1; // minimisation of user-defined function using Genetics Algorithm
    Use["FDA_SA"]          = 0;
    Use["FDA_MC"]          = 0;
    Use["FDA_MT"]          = 0;
    Use["FDA_GAMT"]        = 0;
    Use["FDA_MCMT"]        = 0;
    //
    // --- Neural Networks (all are feed-forward Multilayer Perceptrons)
    Use["MLP"]             = 0; // Recommended ANN
    Use["MLPBFGS"]         = 0; // Recommended ANN with optional training method
    Use["MLPBNN"]          = 1; // Recommended ANN with BFGS training method and bayesian regulator
    Use["CFMlpANN"]        = 0; // Depreciated ANN from ALEPH
    Use["TMlpANN"]         = 0; // ROOT's own ANN
    //
    // --- Support Vector Machine
    Use["SVM"]             = 1;
    //
    // --- Boosted Decision Trees
    Use["BDT"]             = 1; // uses Adaptive Boost
    Use["BDTG"]            = 0; // uses Gradient Boost
    Use["BDTB"]            = 0; // uses Bagging
    Use["BDTD"]            = 0; // decorrelation + Adaptive Boost
    Use["BDTF"]            = 0; // allow usage of fisher discriminant for node splitting
    //
    // --- Friedman's RuleFit method, ie, an optimised series of cuts ("rules")
    Use["RuleFit"]         = 1;
    // ---------------------------------------------------------------

    std::cout << std::endl;
    std::cout << "==> Start TMVAClassification" << std::endl;

    // Select methods (don't look at this code - not of interest)
    if ("" != myMethodList) {
        for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;

        std::vector<TString> mlist = TMVA::gTools().SplitString( myMethodList, ',' );
        for (UInt_t i=0; i<mlist.size(); i++) {
            std::string regMethod(mlist[i]);

            if (Use.find(regMethod) == Use.end()) {
                std::cout << "Method \"" << regMethod << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
                for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) std::cout << it->first << " ";
                std::cout << std::endl;
                return;
            }
            Use[regMethod] = 1;
        }
    }

    // --------------------------------------------------------------------------------------------------

    // --- Here the preparation phase begins
    TString mySginalList("/r02/lc/xu/tau/PandoraBono/20160126am12/rootCustom/tauAnalysisTemplate_");
    TString myBackgroundList(
        "" 
        );
    // Create a ROOT output file where TMVA will store ntuples, histograms, etc.
    
    TString outfileName("TMVAElectron.root");
    TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

    // Create the factory object. Later you can choose the methods
    // whose performance you'd like to investigate. The factory is
    // the only TMVA object you have to interact with
    //
    // The first argument is the base of the name of all the
    // weightfiles in the directory weight/
    //
    // The second argument is the output file for the training results
    // All TMVA output can be suppressed by removing the "!" (not) in
    // front of the "Silent" argument in the option string
    TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification", outputFile,
            "!V:!Silent:Color:DrawProgressBar:Transformations=I:AnalysisType=Classification" );//;D;P;G,D

    // If you wish to modify default settings
    // (please check "src/Config.h" to see all available global options)
    //    (TMVA::gConfig().GetVariablePlotting()).fTimesRMS = 8.0;
    //    (TMVA::gConfig().GetIONames()).fWeightFileDir = "myWeightDirectory";

    // Define the input variables that shall be used for the MVA training
    // note that you may also use variable expressions, such as: "3*var1/var2*abs(var3)"
    // [all types of expressions that can also be parsed by TTree::Draw( "expression" )]
    //factory->AddVariable( "myvar1 := var1+var2", 'F' );
    //factory->AddVariable( "myvar2 := var1-var2", "Expression 2", "", 'F' );
    //factory->AddVariable( "var3",                "Variable 3", "units", 'F' );
    //factory->AddVariable( "var4",                "Variable 4", "units", 'F' );

    
    //factory->AddVariable( "nElectron","nElectron", "", 'I' );
    factory->AddVariable( "nPhoton","nPhoton", "", 'I' );
    //factory->AddVariable( "mVis","mVis", "", 'D' );
    factory->AddVariable( "mPhoton","mPhoton", "", 'D' );
    //factory->AddVariable( "EHCalRatio:= EEHCalRatio > 1 ? 0 : EEHCalRatio","EEHCalRatio", "", 'D' );
    
    //factory->AddVariable( "eBeamElectronPlus := nR1_0_6jet_btag2_Recoil_M * nR1_0_6jet_btag2_Recoil_M / 2 / (nR1_0_6jet_btag2_Recoil_E - nR1_0_6jet_btag2_Recoil_Pz)","energy of a potential electron in plus beam pipe", "", 'D' );
    //factory->AddVariable( "eBeamElectronMinus := nR1_0_6jet_btag2_Recoil_M * nR1_0_6jet_btag2_Recoil_M / 2 / (nR1_0_6jet_btag2_Recoil_E + nR1_0_6jet_btag2_Recoil_Pz)","energy of a potential electron in minus beam pipe", "", 'D' );
    
 
    // You can add so-called "Spectator variables", which are not used in the MVA training,
    // but will appear in the final "TestTree" produced by TMVA. This TestTree will contain the
    // input variables, the response values of all trained MVAs, and the spectator variables
    //factory->AddSpectator( "spec1 := var1*2",  "Spectator 1", "units", 'F' );
    //factory->AddSpectator( "spec2 := var1*3",  "Spectator 2", "", 'F' );

    factory->AddSpectator( "eventType",  "Event Type", "", 'I' );
    factory->AddSpectator( "photonEC",  "photon early conversion flag", "", 'I' );
    factory->AddSpectator( "thrustPrinciple",  "thrustPrinciple", "", 'D' );
    factory->AddSpectator( "mcCloseToZ",  "mcCloseToZ", "", 'I' );
    
    // Read training and test data
    // (it is also possible to use ASCII format as input -> see TMVA Users Guide)

    /*
    TString fnameSg = "./pan_signal_cut.root";
    TString fnameBg = "./pan_bg_cut.root";
    if (gSystem->AccessPathName( fnameSg ))  // file does not exist in local directory
     throw;
    if (gSystem->AccessPathName( fnameBg ))  // file does not exist in local directory
     throw;

    TFile *inputSg = TFile::Open( fnameSg );
    TFile *inputBg = TFile::Open( fnameBg );
    */


    /*
       TChain  *signal= new TChain(treeName.Data());
        //chain ;
        for (int iter_1=0;iter_1<nFiles;iter_1++){
            TString name = namePrefix + TString::Format("%d",iter_1)+nameSignalSuffix;
            signal->Add(name);
        }
    */
    
    
    TStringVec filters, signalFilesNames;
    filters.push_back(".root");
    (void) GetFileList(mySginalList, filters, signalFilesNames);
    
    const TString treeName= "sel";
    TChain  *pSignalChain= new TChain(treeName.Data());
    for (TStringVec::const_iterator iter = signalFilesNames.begin(), iterEnd = signalFilesNames.end(); iter != iterEnd; ++iter)
    {   
        pSignalChain->Add(*iter);
    }
    TChainDoubleMap signalWeightMap;
    if (!signalWeightMap.insert(TChainDoubleMap::value_type(pSignalChain, GetWeight(mySginalList, pSignalChain))).second)
        std::cout<<" Fail to insert in to signal weight map"<<mySginalList <<std::endl;
    
    TChainDoubleMap bkgWeightMap;
    if ("" != myBackgroundList)
    {
        TStringVec bkgFolderVec = TMVA::gTools().SplitString( myBackgroundList, ',' );
        for (TStringVec::const_iterator iter = bkgFolderVec.begin(), iterEnd = bkgFolderVec.end(); iter != iterEnd; ++iter)
        {   
            TStringVec bkgFilters, bkgFilesNames;
            bkgFilters.push_back(".root");
            const TString bkgFolderName(*iter);
            (void) GetFileList(bkgFolderName, bkgFilters, bkgFilesNames);
            TChain  *pBkgChain= new TChain(treeName.Data());
            for (TStringVec::const_iterator jIter = bkgFilesNames.begin(), jIterEnd = bkgFilesNames.end(); jIter != jIterEnd; ++jIter)
            { 
                pBkgChain->Add(*jIter);
            }
            if (!bkgWeightMap.insert(TChainDoubleMap::value_type(pBkgChain, GetWeight(bkgFolderName, pBkgChain))).second)
                std::cout<<" Fail to insert in to background weight map"<<bkgFolderName <<std::endl;
        }
    }

    
    std::cout << "--- TMVAClassification       : Using input file: " << mySginalList + ":" + myBackgroundList << std::endl;


    // --- Register the training and test trees

    //TTree *signal     = (TTree*)inputSg->Get("pan");
    //TTree *background = (TTree*)inputBg->Get("pan");

    // global event weights per tree (see below for setting event-wise weights)
    //Double_t signalWeight     = 1;
    //Double_t backgroundWeight = 1;

    // You can add an arbitrary number of signal or background trees

    for (TChainDoubleMap::iterator iter = signalWeightMap.begin(), iterEnd = signalWeightMap.end(); iter != iterEnd; ++iter)
    {
        factory->AddSignalTree    (iter->first,iter->second);
        factory->AddBackgroundTree    (iter->first, iter->second);
    }
    //for (TChainDoubleMap::iterator iter = bkgWeightMap.begin(), iterEnd = bkgWeightMap.end(); iter != iterEnd; ++iter)
    //{
    //    factory->AddBackgroundTree    (iter->first, iter->second);
    //}
    
    //factory->AddSignalTree    ( signal,     signalWeight     );
    //factory->AddBackgroundTree( background, backgroundWeight );

    // To give different trees for training and testing, do as follows:
    //    factory->AddSignalTree( signalTrainingTree, signalTrainWeight, "Training" );
    //    factory->AddSignalTree( signalTestTree,     signalTestWeight,  "Test" );

    // Use the following code instead of the above two or four lines to add signal and background
    // training and test events "by hand"
    // NOTE that in this case one should not give expressions (such as "var1+var2") in the input
    //      variable definition, but simply compute the expression before adding the event
    //
    //     // --- begin ----------------------------------------------------------
    //     std::vector<Double_t> vars( 4 ); // vector has size of number of input variables
    //     Float_t  treevars[4], weight;
    //
    //     // Signal
    //     for (UInt_t ivar=0; ivar<4; ivar++) signal->SetBranchAddress( Form( "var%i", ivar+1 ), &(treevars[ivar]) );
    //     for (UInt_t i=0; i<signal->GetEntries(); i++) {
    //        signal->GetEntry(i);
    //        for (UInt_t ivar=0; ivar<4; ivar++) vars[ivar] = treevars[ivar];
    //        // add training and test events; here: first half is training, second is testing
    //        // note that the weight can also be event-wise
    //        if (i < signal->GetEntries()/2.0) factory->AddSignalTrainingEvent( vars, signalWeight );
    //        else                              factory->AddSignalTestEvent    ( vars, signalWeight );
    //     }
    //
    //     // Background (has event weights)
    //     background->SetBranchAddress( "weight", &weight );
    //     for (UInt_t ivar=0; ivar<4; ivar++) background->SetBranchAddress( Form( "var%i", ivar+1 ), &(treevars[ivar]) );
    //     for (UInt_t i=0; i<background->GetEntries(); i++) {
    //        background->GetEntry(i);
    //        for (UInt_t ivar=0; ivar<4; ivar++) vars[ivar] = treevars[ivar];
    //        // add training and test events; here: first half is training, second is testing
    //        // note that the weight can also be event-wise
    //        if (i < background->GetEntries()/2) factory->AddBackgroundTrainingEvent( vars, backgroundWeight*weight );
    //        else                                factory->AddBackgroundTestEvent    ( vars, backgroundWeight*weight );
    //     }
    // --- end ------------------------------------------------------------
    //
    // --- end of tree registration

    // Set individual event weights (the variables must exist in the original TTree)
    //    for signal    : factory->SetSignalWeightExpression    ("weight1*weight2");
    //    for background: factory->SetBackgroundWeightExpression("weight1*weight2");

    //factory->SetBackgroundWeightExpression( "weight" );

    // Apply additional cuts on the signal and background samples (can be different) 
    //  !(nElectron == 1 && EEHCalRatio > 0.95) && !(nMuon == 1) && !(nPionCharge > 1) && ! (nPionCharge == 1 && nPhoton < 1)
    const TCut selectionCut = " thrustPrinciple > 0.99 && photonEC < 1 && mcCloseToZ < 1 && nPionCharge == 1 && nPhoton >0";
    TCut mycut = "TMath::Abs(eventType - 5 ) < 0.2 " + selectionCut;//"e_2>1 && d<150 && dLayer<200 && duplicated<1";//"eventType<1.0"; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
    TCut mycutb ="eventType > 0 && eventType != 5" + selectionCut; //"e_2>1 && d<150 && dLayer<200 && duplicated<1";//"!( eventType<1.0) " ; // for example: TCut mycutb = "abs(var1)<0.5";
    // same cut 


    // Tell the factory how to use the training and testing events
    //
    // If no numbers of events are given, half of the events in the tree are used
    // for training, and the other half for testing:
    //    factory->PrepareTrainingAndTestTree( mycut, "SplitMode=random:!V" );
    // To also specify the number of testing events, use:
    //    factory->PrepareTrainingAndTestTree( mycut,
    //                                         "NSigTrain=3000:NBkgTrain=3000:NSigTest=3000:NBkgTest=3000:SplitMode=Random:!V" );
    factory->PrepareTrainingAndTestTree( mycut, mycutb,
                                         "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );

    // ---- Book MVA methods
    //
    // Please lookup the various method configuration options in the corresponding cxx files, eg:
    // src/MethoCuts.cxx, etc, or here: http://tmva.sourceforge.net/optionRef.html
    // it is possible to preset ranges in the option string in which the cut optimisation should be done:
    // "...:CutRangeMin[2]=-1:CutRangeMax[2]=1"...", where [2] is the third input variable

    // Cut optimisation
    if (Use["Cuts"])
        factory->BookMethod( TMVA::Types::kCuts, "Cuts",
                             "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart" );

    if (Use["CutsD"])
        factory->BookMethod( TMVA::Types::kCuts, "CutsD",
                             "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=Decorrelate" );

    if (Use["CutsPCA"])
        factory->BookMethod( TMVA::Types::kCuts, "CutsPCA",
                             "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=PCA" );

    if (Use["CutsGA"])
        factory->BookMethod( TMVA::Types::kCuts, "CutsGA",
                             "H:!V:FitMethod=GA:CutRangeMin[0]=-10:CutRangeMax[0]=10:VarProp[1]=FMax:EffSel:Steps=30:Cycles=3:PopSize=400:SC_steps=10:SC_rate=5:SC_factor=0.95" );

    if (Use["CutsSA"])
        factory->BookMethod( TMVA::Types::kCuts, "CutsSA",
                             "!H:!V:FitMethod=SA:EffSel:MaxCalls=150000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale" );

    // Likelihood ("naive Bayes estimator")
    if (Use["Likelihood"])
        factory->BookMethod( TMVA::Types::kLikelihood, "Likelihood",
                             "H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmoothBkg[1]=10:NSmooth=1:NAvEvtPerBin=50" );

    // Decorrelated likelihood
    if (Use["LikelihoodD"])
        factory->BookMethod( TMVA::Types::kLikelihood, "LikelihoodD",
                             "!H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=Decorrelate" );

    // PCA-transformed likelihood
    if (Use["LikelihoodPCA"])
        factory->BookMethod( TMVA::Types::kLikelihood, "LikelihoodPCA",
                             "!H:!V:!TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=PCA" );

    // Use a kernel density estimator to approximate the PDFs
    if (Use["LikelihoodKDE"])
        factory->BookMethod( TMVA::Types::kLikelihood, "LikelihoodKDE",
                             "!H:!V:!TransformOutput:PDFInterpol=KDE:KDEtype=Gauss:KDEiter=Adaptive:KDEFineFactor=0.3:KDEborder=None:NAvEvtPerBin=50" );

    // Use a variable-dependent mix of splines and kernel density estimator
    if (Use["LikelihoodMIX"])
        factory->BookMethod( TMVA::Types::kLikelihood, "LikelihoodMIX",
                             "!H:!V:!TransformOutput:PDFInterpolSig[0]=KDE:PDFInterpolBkg[0]=KDE:PDFInterpolSig[1]=KDE:PDFInterpolBkg[1]=KDE:PDFInterpolSig[2]=Spline2:PDFInterpolBkg[2]=Spline2:PDFInterpolSig[3]=Spline2:PDFInterpolBkg[3]=Spline2:KDEtype=Gauss:KDEiter=Nonadaptive:KDEborder=None:NAvEvtPerBin=50" );

    // Test the multi-dimensional probability density estimator
    // here are the options strings for the MinMax and RMS methods, respectively:
    //      "!H:!V:VolumeRangeMode=MinMax:DeltaFrac=0.2:KernelEstimator=Gauss:GaussSigma=0.3" );
    //      "!H:!V:VolumeRangeMode=RMS:DeltaFrac=3:KernelEstimator=Gauss:GaussSigma=0.3" );
    if (Use["PDERS"])
        factory->BookMethod( TMVA::Types::kPDERS, "PDERS",
                             "!H:!V:NormTree=T:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600" );

    if (Use["PDERSD"])
        factory->BookMethod( TMVA::Types::kPDERS, "PDERSD",
                             "!H:!V:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600:VarTransform=Decorrelate" );

    if (Use["PDERSPCA"])
        factory->BookMethod( TMVA::Types::kPDERS, "PDERSPCA",
                             "!H:!V:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600:VarTransform=PCA" );

    // Multi-dimensional likelihood estimator using self-adapting phase-space binning
    if (Use["PDEFoam"])
        factory->BookMethod( TMVA::Types::kPDEFoam, "PDEFoam",
                             "!H:!V:SigBgSeparate=F:TailCut=0.001:VolFrac=0.0666:nActiveCells=500:nSampl=2000:nBin=5:Nmin=100:Kernel=None:Compress=T" );

    if (Use["PDEFoamBoost"])
        factory->BookMethod( TMVA::Types::kPDEFoam, "PDEFoamBoost",
                             "!H:!V:Boost_Num=30:Boost_Transform=linear:SigBgSeparate=F:MaxDepth=4:UseYesNoCell=T:DTLogic=MisClassificationError:FillFoamWithOrigWeights=F:TailCut=0:nActiveCells=500:nBin=20:Nmin=400:Kernel=None:Compress=T" );

    // K-Nearest Neighbour classifier (KNN)
    if (Use["KNN"])
        factory->BookMethod( TMVA::Types::kKNN, "KNN",
                             "H:nkNN=20:ScaleFrac=0.8:SigmaFact=1.0:Kernel=Gaus:UseKernel=F:UseWeight=T:!Trim" );

    // H-Matrix (chi2-squared) method
    if (Use["HMatrix"])
        factory->BookMethod( TMVA::Types::kHMatrix, "HMatrix", "!H:!V:VarTransform=None" );

    // Linear discriminant (same as Fisher discriminant)
    if (Use["LD"])
        factory->BookMethod( TMVA::Types::kLD, "LD", "H:!V:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10" );

    // Fisher discriminant (same as LD)
    if (Use["Fisher"])
        factory->BookMethod( TMVA::Types::kFisher, "Fisher", "H:!V:Fisher:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10" );

    // Fisher with Gauss-transformed input variables
    if (Use["FisherG"])
        factory->BookMethod( TMVA::Types::kFisher, "FisherG", "H:!V:VarTransform=Gauss" );

    // Composite classifier: ensemble (tree) of boosted Fisher classifiers
    if (Use["BoostedFisher"])
        factory->BookMethod( TMVA::Types::kFisher, "BoostedFisher",
                             "H:!V:Boost_Num=20:Boost_Transform=log:Boost_Type=AdaBoost:Boost_AdaBoostBeta=0.2:!Boost_DetailedMonitoring" );

    // Function discrimination analysis (FDA) -- test of various fitters - the recommended one is Minuit (or GA or SA)
    if (Use["FDA_MC"])
        factory->BookMethod( TMVA::Types::kFDA, "FDA_MC",
                             "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MC:SampleSize=100000:Sigma=0.1" );

    if (Use["FDA_GA"]) // can also use Simulated Annealing (SA) algorithm (see Cuts_SA options])
        factory->BookMethod( TMVA::Types::kFDA, "FDA_GA",
                             "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:PopSize=300:Cycles=3:Steps=20:Trim=True:SaveBestGen=1" );

    if (Use["FDA_SA"]) // can also use Simulated Annealing (SA) algorithm (see Cuts_SA options])
        factory->BookMethod( TMVA::Types::kFDA, "FDA_SA",
                             "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=SA:MaxCalls=15000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale" );

    if (Use["FDA_MT"])
        factory->BookMethod( TMVA::Types::kFDA, "FDA_MT",
                             "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=2:UseImprove:UseMinos:SetBatch" );

    if (Use["FDA_GAMT"])
        factory->BookMethod( TMVA::Types::kFDA, "FDA_GAMT",
                             "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:Converger=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=0:!UseImprove:!UseMinos:SetBatch:Cycles=1:PopSize=5:Steps=5:Trim" );

    if (Use["FDA_MCMT"])
        factory->BookMethod( TMVA::Types::kFDA, "FDA_MCMT",
                             "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MC:Converger=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=0:!UseImprove:!UseMinos:SetBatch:SampleSize=20" );

    // TMVA ANN: MLP (recommended ANN) -- all ANNs in TMVA are Multilayer Perceptrons
    if (Use["MLP"])
        factory->BookMethod( TMVA::Types::kMLP, "MLP", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:!UseRegulator" );

    if (Use["MLPBFGS"])
        factory->BookMethod( TMVA::Types::kMLP, "MLPBFGS", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:!UseRegulator" );

    if (Use["MLPBNN"])
        factory->BookMethod( TMVA::Types::kMLP, "MLPBNN", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:UseRegulator" ); // BFGS training with bayesian regulators

    // CF(Clermont-Ferrand)ANN
    if (Use["CFMlpANN"])
        factory->BookMethod( TMVA::Types::kCFMlpANN, "CFMlpANN", "!H:!V:NCycles=2000:HiddenLayers=N+1,N"  ); // n_cycles:#nodes:#nodes:...

    // Tmlp(Root)ANN
    if (Use["TMlpANN"])
        factory->BookMethod( TMVA::Types::kTMlpANN, "TMlpANN", "!H:!V:NCycles=200:HiddenLayers=N+1,N:LearningMethod=BFGS:ValidationFraction=0.3"  ); // n_cycles:#nodes:#nodes:...

    // Support Vector Machine
    if (Use["SVM"])
        factory->BookMethod( TMVA::Types::kSVM, "SVM", "Gamma=0.25:Tol=0.001:VarTransform=Norm" );

    // Boosted Decision Trees
    if (Use["BDTG"]) // Gradient Boost
        factory->BookMethod( TMVA::Types::kBDT, "BDTG",
                             "!H:!V:NTrees=1000:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=2" );

    if (Use["BDT"])  // Adaptive Boost
        factory->BookMethod( TMVA::Types::kBDT, "BDT",
                             "!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20:DoPreselection=True" );

    if (Use["BDTB"]) // Bagging
        factory->BookMethod( TMVA::Types::kBDT, "BDTB",
                             "!H:!V:NTrees=400:BoostType=Bagging:SeparationType=GiniIndex:nCuts=20" );

    if (Use["BDTD"]) // Decorrelation + Adaptive Boost
        factory->BookMethod( TMVA::Types::kBDT, "BDTD",
                             "!H:!V:NTrees=400:MinNodeSize=5%:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:VarTransform=Decorrelate" );

    if (Use["BDTF"])  // Allow Using Fisher discriminant in node splitting for (strong) linearly correlated variables
        factory->BookMethod( TMVA::Types::kBDT, "BDTMitFisher",
                             "!H:!V:NTrees=50:MinNodeSize=2.5%:UseFisherCuts:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20" );

    // RuleFit -- TMVA implementation of Friedman's method
    if (Use["RuleFit"])
        factory->BookMethod( TMVA::Types::kRuleFit, "RuleFit",
                             "H:!V:RuleFitModule=RFTMVA:Model=ModRuleLinear:MinImp=0.001:RuleMinDist=0.001:NTrees=20:fEventsMin=0.01:fEventsMax=0.5:GDTau=-1.0:GDTauPrec=0.01:GDStep=0.01:GDNSteps=10000:GDErrScale=1.02" );

    // For an example of the category classifier usage, see: TMVAClassificationCategory

    // --------------------------------------------------------------------------------------------------

    // ---- Now you can optimize the setting (configuration) of the MVAs using the set of training events

    // ---- STILL EXPERIMENTAL and only implemented for BDT's !
    // factory->OptimizeAllMethods("SigEffAt001","Scan");
    // factory->OptimizeAllMethods("ROCIntegral","FitGA");

    // --------------------------------------------------------------------------------------------------

    // ---- Now you can tell the factory to train, test, and evaluate the MVAs

    // Train MVAs using the set of training events
    factory->TrainAllMethods();

    // ---- Evaluate all MVAs using the set of test events
    factory->TestAllMethods();

    // ----- Evaluate and compare performance of all configured MVAs
    factory->EvaluateAllMethods();

    // --------------------------------------------------------------

    // Save the output
    outputFile->Close();

    std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
    std::cout << "==> TMVAClassification is done!" << std::endl;

    delete factory;

    // Launch the GUI for the root macros
    if (!gROOT->IsBatch()) TMVAGui( outfileName );
}
