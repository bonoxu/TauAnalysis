#include <iostream>
#include <cmath>
#include <limits>
#include <iomanip>

#include "TROOT.h"

#include <TChain.h>

#include <TApplication.h>

#include <TCut.h>

enum EVENT_TYPE
{
    E = 1,
    MU = 2,
    PION = 3,
    PION2PHOTON = 4,
    PION4PHOTON = 5,
    PION2PION = 6,
    OTHER = 0
};
typedef std::vector<TString> TStringVec;
typedef std::vector<double> DoubleVec;
typedef std::map<TString, DoubleVec> TStringDoubleVecMap;
typedef std::map<EVENT_TYPE, TString> EvtTypeNameMap;
typedef std::map<TString, TChain *> StrTChainMap;
typedef std::vector<TCut> TCutVec;
typedef std::map<TString, double> TStringDoubleMap;

//------------------------------------------------------------------------------------------------------------------------------------------

void initialiseMaps(const TString &id, const TString &tfile, const TString &ttree, const double totalNEvent, StrTChainMap &idChainMap, TStringDoubleMap &idWeightMap)
{
    TChain *pTChain= new TChain(ttree.Data());
    pTChain->Add(tfile.Data());
    idChainMap.insert(StrTChainMap::value_type(id, pTChain));
    const double weight(totalNEvent / static_cast<double>(pTChain->GetEntries())); // see https://twiki.cern.ch/twiki/bin/view/CLIC/MonteCarloSamplesForTheHiggsPaper
    idWeightMap.insert(TStringDoubleMap::value_type(id, weight));
}

//------------------------------------------------------------------------------------------------------------------------------------------

TString getName(const int id)
{
    switch (id)
    {
        case 0:
            return "Back";
        case 1:
            return "E";
        case 2:
            return "Mu";
        case 3:
            return "1Pion";
        case 4:
            return "Pion2Y";
        case 5:
            return "Pion4Y";
        case 6:
            return "3Pion";
        default:
            return "unknown";
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void getNumbers(StrTChainMap &idChainMap, TStringDoubleMap &idWeightMap, const TCutVec &selectionCuts, const TCutVec &tmvaSelectionCuts,
    const TString &tmvaRootFile)
{
    TCut currentSelectionCut;
    
    const int columnWidth(14);
    std::cout << std::setprecision(2) << std::fixed;
    for (int i = 0, iEnd = 7; i < iEnd; ++i)
    {
        std::cout<< std::setw(columnWidth)<<getName(i)<<";";
    }
    std::cout<<std::endl;
    for (StrTChainMap::const_iterator jIter = idChainMap.begin(), jIterEnd = idChainMap.end(); jIter != jIterEnd; ++jIter)
    {
        std::cout<<jIter->first<<std::endl;
        currentSelectionCut = "";
        const double weight(idWeightMap.at(jIter->first));
        //const double weight(1.f);
        for (TCutVec::const_iterator iter = selectionCuts.begin(), iterEnd = selectionCuts.end(); iter != iterEnd; ++iter)
        {
            currentSelectionCut+=*iter;
            const double nTotal(jIter->second->GetEntries(currentSelectionCut) * weight);
            for (int i = 0, iEnd = 7; i < iEnd; ++i)
            {
                const TCut eventType(TString::Format("TMath::Abs(eventType - %d ) < 0.2" , i));
                const double signal(jIter->second->GetEntries(eventType + currentSelectionCut) * weight);
                std::cout<<std::setw(columnWidth)<<TString::Format("%4.2f,%4.3f", signal , signal / nTotal)<<";";
            }
            std::cout<<*iter<<std::endl;
        }
    }

}

//------------------------------------------------------------------------------------------------------------------------------------------

int main(int argc, char* argv[]) {
    
    //TApplication* pRootapp = new TApplication("example",&argc, argv);
    
    const TString defaultTreeName("sel");
    StrTChainMap idChainMap;
    TStringDoubleMap idWeightMap;
    //initialiseMaps("improved","/r02/lc/xu/tau/PandoraBono/20160126am12/rootCustom/tauAnalysisTemplate_*.root", defaultTreeName, 428.2776, idChainMap, idWeightMap);
    //initialiseMaps("cheat","/r02/lc/xu/tau/PandoraBono/20160212am19cheat/rootCustom/tauAnalysisTemplate_*.root", defaultTreeName, 0, idChainMap, idWeightMap);
    //initialiseMaps("default","/r02/lc/xu/tau/PandoraDefault/20160213am21/rootCustom/tauAnalysisTemplate_*.root", defaultTreeName, 428.2776, idChainMap, idWeightMap);
    //initialiseMaps("split","/r02/lc/xu/tau/PandoraBonoSplit/20160216am19/rootCustom/tauAnalysisTemplate_*.root", defaultTreeName, 500, idChainMap, idWeightMap);
    //initialiseMaps("100improved","/r02/lc/xu/tau/PandoraBonoNew/100GeV/rootCustom/tauAnalysisTemplate_*.root", defaultTreeName, 52869, idChainMap, idWeightMap);
    //initialiseMaps("100default","/r02/lc/xu/tau/PandoraDefaultNew/100GeV/rootCustom/tauAnalysisTemplate_*.root", defaultTreeName, 52869, idChainMap, idWeightMap);
    initialiseMaps("200improved","/r02/lc/xu/tau/PandoraBonoNew/200GeV/rootCustom/tauAnalysisTemplate_*.root", defaultTreeName, 2844.395, idChainMap, idWeightMap);
    initialiseMaps("200default","/r02/lc/xu/tau/PandoraDefaultNew/200GeV/rootCustom/tauAnalysisTemplate_*.root", defaultTreeName, 2844.395, idChainMap, idWeightMap);

    
    TCutVec muonSelectionCuts, tmvaDummyCuts;
    TString tmvaDummyRoot("");
    muonSelectionCuts.push_back("thrustPrinciple > 0.99 && photonEC < 1 && mcCloseToZ < 1 && nMuon == 1 ");
    (void) getNumbers(idChainMap, idWeightMap, muonSelectionCuts, tmvaDummyCuts,tmvaDummyRoot);

    TCutVec electronSelectionCuts;
    electronSelectionCuts.push_back("thrustPrinciple > 0.99 && photonEC < 1 && mcCloseToZ < 1 && !(nMuon == 1) && nElectron == 1 && EEHCalRatio  > 0.95");
    (void) getNumbers(idChainMap, idWeightMap, electronSelectionCuts, tmvaDummyCuts,tmvaDummyRoot);

    TCutVec pion3SelectionCuts;
    pion3SelectionCuts.push_back("thrustPrinciple > 0.99 && photonEC < 1 && mcCloseToZ < 1 && !(nMuon == 1) && !(nElectron == 1 && EEHCalRatio  > 0.95) && (nPionCharge > 1) ");
    (void) getNumbers(idChainMap, idWeightMap, pion3SelectionCuts, tmvaDummyCuts,tmvaDummyRoot);

    TCutVec pionSelectionCuts;
    pionSelectionCuts.push_back("thrustPrinciple > 0.99 && photonEC < 1 && mcCloseToZ < 1 && !(nMuon == 1) && !(nElectron == 1 && EEHCalRatio  > 0.95) && (nPionCharge == 1) && nPhoton == 0 ");
    (void) getNumbers(idChainMap, idWeightMap, pionSelectionCuts, tmvaDummyCuts,tmvaDummyRoot);

    TCutVec pion2YSelectionCuts;
    //pion2YSelectionCuts.push_back("thrustPrinciple > 0.99 && photonEC < 1 && mcCloseToZ < 1 && !(nMuon == 1) && !(nElectron == 1 && EEHCalRatio  > 0.95) && nPionCharge == 1 ");
    pion2YSelectionCuts.push_back("thrustPrinciple > 0.99 && photonEC < 1 && mcCloseToZ < 1 && !(nMuon == 1) && !(nElectron == 1 && EEHCalRatio  > 0.95) && nPionCharge == 1 && (nPhoton == 1 || (nPhoton == 2 && (mPhoton < 0.3 && mPhoton > 0.05) ))");
    (void) getNumbers(idChainMap, idWeightMap, pion2YSelectionCuts, tmvaDummyCuts,tmvaDummyRoot);

    TCutVec pion4YSelectionCuts;
    pion4YSelectionCuts.push_back("thrustPrinciple > 0.99 && photonEC < 1 && mcCloseToZ < 1 && !(nMuon == 1) && !(nElectron == 1 && EEHCalRatio  > 0.95) && nPionCharge == 1 && nPhoton > 1 && !(nPhoton == 2 && mPhoton < 0.3 && mPhoton > 0.05 ) && mPhoton >= 0.3 ");
    (void) getNumbers(idChainMap, idWeightMap, pion4YSelectionCuts, tmvaDummyCuts,tmvaDummyRoot);


    TCutVec selectionCuts;
    selectionCuts.push_back("");
    selectionCuts.push_back("thrustPrinciple > 0.99");
    
    selectionCuts.push_back("photonEC < 1");
    selectionCuts.push_back("mcCloseToZ < 1");
    //selectionCuts.push_back("nElectron == 1");// ELECTRON  nElectron == 1 && EEHCalRatio  > 0.95
    //selectionCuts.push_back("nCharge > 0");
    //selectionCuts.push_back("nCharge > 1");
    //selectionCuts.push_back("nElectron == 1");
    //selectionCuts.push_back("nMuon >= 1");
    //selectionCuts.push_back("nMuon == 1");// MUON nMuon == 1
    //selectionCuts.push_back("nPionCharge >= 3");
    //selectionCuts.push_back("nPionCharge == 3");
    
    //selectionCuts.push_back("nPionCharge == 1");
    //selectionCuts.push_back("EEHCalRatio  > 0.05");
    selectionCuts.push_back("! (nMuon == 1 )");
    selectionCuts.push_back("! (nElectron == 1 && EEHCalRatio  > 0.95)");
    //selectionCuts.push_back("nPionCharge > 0"); // 3PION nPionCharge > 1
    selectionCuts.push_back("! (nPionCharge > 1) ");
    selectionCuts.push_back("(nPionCharge == 1) ");
    selectionCuts.push_back("nPhoton > 0");
    selectionCuts.push_back("nPhoton > 1");
    selectionCuts.push_back("nPhoton > 2");
    selectionCuts.push_back("nPhoton > 3");
    selectionCuts.push_back("nPhoton > 4");
    //selectionCuts.push_back("! (nPionCharge == 1 && nPhoton < 1) "); // PION nPionCharge == 1 && nPhoton < 1
    //selectionCuts.push_back("! (nPionCharge == 1 && nPhoton < 3) "); // PION2Y nPionCharge == 1 && nPhoton < 3
    //selectionCuts.push_back("nPionCharge > 2");
    //selectionCuts.push_back("R0 < 1");
    //selectionCuts.push_back("R0 < 0.5");
    //selectionCuts.push_back("R0 < 0.2");
    //selectionCuts.push_back("nPfo > 0");
    //selectionCuts.push_back("nPfo > 1");
    //selectionCuts.push_back("nPfo > 2");
    //selectionCuts.push_back("nPfo > 3");
    //selectionCuts.push_back("nPfo > 4");
    //selectionCuts.push_back("nPfo > 5");
    //selectionCuts.push_back("nPfo > 6");
    //selectionCuts.push_back("nCharge > 0");
    //selectionCuts.push_back("nCharge > 1");
    //selectionCuts.push_back("nCharge > 2");
    //selectionCuts.push_back("nCharge > 3");
    //selectionCuts.push_back("(nPionCharge == 1 && nPhoton > 2)");
    //selectionCuts.push_back("nPhoton > 3");
    //selectionCuts.push_back("nPhoton > 4");
    //selectionCuts.push_back("nPhoton > 2");
    //selectionCuts.push_back("nPhoton > 3");
    //selectionCuts.push_back("nPhoton > 4");
    TCutVec tmvaSelectionCuts;
    //tmvaSelectionCuts.push_back("BDT>0.2587");
    // Get BDT cuts
    TString tmvaRootR0_7("TMVA20151211R0_7_btag2_full3.root");
    
    std::cout<<"Stats for selected R=0.7 "<<tmvaRootR0_7<<std::endl;
    (void) getNumbers(idChainMap, idWeightMap, selectionCuts, tmvaSelectionCuts,tmvaRootR0_7);
    //pRootapp->Run();

    return 0;

}

