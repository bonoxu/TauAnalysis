#include <iostream>
#include <cmath>
#include <limits>
#include <iomanip>

#include "TROOT.h"

#include <TChain.h>

#include <TApplication.h>

#include <TCut.h>

enum TAU_DECAY_FINAL_STATES
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

TString getName(const TString &id)
{
    if ("6022" == id)
        return "HH";
    else if ("5527" == id)
        return "qqqqvv";
    else if ("5572" == id)
        return "qqqqll";
    else if ("5594" == id)
        return "qqqqlv";
    else if ("4034" == id)
        return "qqqq";
    else if ("3366" == id)
        return "eY:qqqqe,EPA";
    else if ("3369" == id)
        return "eY:qqqqe,BS";
    else if ("3372" == id)
        return "Ye:qqqqe,EPA";
    else if ("3375" == id)
        return "Ye:qqqqe,BS";
    else if ("3378" == id)
        return "eY:qqqqv,EPA";
    else if ("3381" == id)
        return "eY:qqqqv,BS";
    else if ("3384" == id)
        return "Ye:qqqqv,EPA";
    else if ("3387" == id)
        return "Ye:qqqqv,BS";
    else if ("3414" == id)
        return "YY:qqqq,2EPA";
    else if ("3417" == id)
        return "YY:qqqq,EPA,BS";
    else if ("3420" == id)
        return "YY:qqqq,BS,EPA";
    else if ("3423" == id)
        return "YY:qqqq,2BS";
    else if ("3243" == id)
        return "qqvv";
    else if ("3246" == id)
        return "qqll";
    else if ("3249" == id)
        return "qqlv";
    else
        return id;
}

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

TCut getSampleCut(const TString &id)
{
    return TCut("TMath::Abs(prodID -" + id + ") < 0.2");
}

//------------------------------------------------------------------------------------------------------------------------------------------

void getNumbers(StrTChainMap &idChainMap, TStringDoubleMap &idWeightMap, const TCutVec &selectionCuts, const TCutVec &tmvaSelectionCuts,
    const TString &tmvaRootFile)
{
    const TCut signalCut = "TMath::Abs(eventType -1 ) < 0.2";
    const TCut bkgCut = "TMath::Abs(eventType -1 ) > 0.2 ";
    TCut currentSelectionCut;
    
    const int columnWidth(14);
    std::cout << std::setprecision(2) << std::fixed;
    for (StrTChainMap::reverse_iterator jIter = idChainMap.rbegin(), jIterEnd = idChainMap.rend(); jIter != jIterEnd; ++jIter)
    {
        if ("6022" == jIter->first)
            std::cout<< std::setw(columnWidth)<<getName(jIter->first) + "->bbWW"<<";"<< std::setw(columnWidth)<<getName(jIter->first) + "->other"<<";";
        else
            std::cout<< std::setw(columnWidth)<<getName(jIter->first)<<";";
    }
    std::cout<<std::endl;
    
    for (TCutVec::const_iterator iter = selectionCuts.begin(), iterEnd = selectionCuts.end(); iter != iterEnd; ++iter)
    {
        currentSelectionCut+=*iter;
        double nSgl(0.f), nBgk(0.f);
        for (StrTChainMap::reverse_iterator jIter = idChainMap.rbegin(), jIterEnd = idChainMap.rend(); jIter != jIterEnd; ++jIter)
        {
            //const double weight(idWeightMap.at(jIter->first));
            const double weight(1.f);
            if ("6022" == jIter->first)
            {
                const double signal(jIter->second->GetEntries(signalCut + currentSelectionCut) * weight);
                const double bkg(jIter->second->GetEntries(bkgCut + currentSelectionCut) * weight);
                std::cout<<std::setw(columnWidth)<<signal<<";"<<
                    std::setw(columnWidth)<<bkg<<";";
                nSgl += signal;
                nBgk += bkg;
            }

            else
            {
                const double bkg(jIter->second->GetEntries(currentSelectionCut) * weight);
                std::cout<<std::setw(columnWidth)<<bkg<<";";
                nBgk += bkg;
            }
        }
        std::cout<<std::setw(columnWidth)<<"S:" + TString::Format("%f",nSgl/std::sqrt(nBgk + nSgl) ) + ";";
        std::cout<<std::setw(columnWidth)<<"R:" + TString::Format("%G",nSgl/ nBgk) + ";";
        if ("" != *iter)
        {
            std::cout<<"After selection cut "<<*iter<<std::endl;
        }
        else
            std::cout<<std::endl;
        
    }
    
    TChain *pTChain= new TChain("TestTree");
    pTChain->Add(tmvaRootFile);

    TCut currentTmvaSelectionCut;
    for (TCutVec::const_iterator iter = tmvaSelectionCuts.begin(), iterEnd = tmvaSelectionCuts.end(); iter != iterEnd; ++iter)
    {
        currentSelectionCut+=*iter;
        
        double nSgl(0.f), nBgk(0.f);
        for (TStringDoubleMap::reverse_iterator jIter = idWeightMap.rbegin(), jIterEnd = idWeightMap.rend(); jIter != jIterEnd; ++jIter)
        {
            //const double weight(2 * jIter->second);
            const double weight(1.f);
            if ("6022" == jIter->first)
            {
                const double signal(pTChain->GetEntries(signalCut + currentSelectionCut + getSampleCut(jIter->first)) * weight);
                const double bkg(pTChain->GetEntries(bkgCut + currentSelectionCut + getSampleCut(jIter->first)) * weight);
                std::cout<<std::setw(columnWidth)<<signal<<";"<<
                    std::setw(columnWidth)<<bkg<<";";
                nSgl += signal;
                nBgk += bkg;
            }

            else
            {
                const double bkg(pTChain->GetEntries(currentSelectionCut + getSampleCut(jIter->first)) * weight);
                std::cout<<std::setw(columnWidth)<<bkg<<";";
                nBgk += bkg;
            }
                
        }
        std::cout<<std::setw(columnWidth)<<"S:" + TString::Format("%f",nSgl/std::sqrt(nBgk + nSgl) ) + ";";
        std::cout<<std::setw(columnWidth)<<"R:" + TString::Format("%G",nSgl/ nBgk) + ";";
        std::cout<<"After selection cut "<<*iter<<std::endl;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

int main(int argc, char* argv[]) {
    
    //TApplication* pRootapp = new TApplication("example",&argc, argv);
    const double integratedLumi1400GeV(1500.f); //unit fb
    const double lumiCorrectionEGamma(0.75f),lumiCorrectionGammaE(0.75f), lumiCorrectionGammaGamma(0.64f);
    
    const TString defaultTreeName("sel");
    StrTChainMap idChainMap;
    TStringDoubleMap idWeightMap;
    initialiseMaps("0","/r02/lc/xu/tau/PandoraBono/20160126am12/rootCustom/tauAnalysisTemplate_*.root", defaultTreeName, 0.149 * integratedLumi1400GeV, idChainMap, idWeightMap);

    TCutVec selectionCuts;
    //selectionCuts.push_back("");

    TCutVec tmvaSelectionCuts;
    //tmvaSelectionCuts.push_back("BDT>0.2587");
    // Get BDT cuts
    TString tmvaRootR0_7("TMVA20151211R0_7_btag2_full3.root");
    
    std::cout<<"Stats for selected R=0.7 "<<tmvaRootR0_7<<std::endl;
    //(void) getNumbers(idChainMap, idWeightMap, selectionCuts, tmvaSelectionCuts,tmvaRootR0_7);
    //pRootapp->Run();

    return 0;

}

