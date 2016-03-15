#include <iostream>
#include <cmath>
#include <limits>
#include <iomanip>

#include "TROOT.h"

#include <TChain.h>

#include <TApplication.h>
#include <TColor.h>
#include <TCut.h>
#include <TStyle.h>
#include <TH1.h>
#include "TLegend.h"
#include "TSystem.h"
#include "TPad.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TF1.h"

typedef std::vector<TString> TStringVec;
typedef std::vector<double> DoubleVec;
typedef std::map<TString, DoubleVec> TStringDoubleVecMap;

typedef std::map<TString, TChain *> StrTChainMap;
typedef std::vector<TCut> TCutVec;
typedef std::map<TString, double> TStringDoubleMap;

class Parameters
{
public:
    Parameters(const TString &name, const TString &tfileName, const TString &ttreeName, const TStringDoubleMap &eventTypeTotalWeightMap, const TCut &commonCut);
    
    TChain              *m_pTChain;
    TStringDoubleMap     m_evtTypeWeightMap;
    TString              m_name;
    TCut                 m_commonCut;
};

typedef std::vector<Parameters> ParametersVec;

//------------------------------------------------------------------------------------------------------------------------------------------

Parameters::Parameters(const TString &name, const TString &tfileName, const TString &ttreeName, const TStringDoubleMap &eventTypeTotalWeightMap, const TCut &commonCut):
m_name(name),
m_commonCut(commonCut)
{
    m_pTChain = new TChain(ttreeName.Data());
    m_pTChain->Add(tfileName.Data());
    
    m_evtTypeWeightMap = eventTypeTotalWeightMap;
    for (TStringDoubleMap::const_iterator iter = eventTypeTotalWeightMap.begin(), iterEnd = eventTypeTotalWeightMap.end(); iter != iterEnd; ++iter)
    {
        // ATTN chain needs to be initialised first. This uselessEntries will always be the same as m_pTChain->GetEntries(), no matter what cut is
        if (eventTypeTotalWeightMap.begin() == iter)
        {
            const int uselessEntries(m_pTChain->GetEntries(TCut(iter->first) + commonCut));
            const int uselessEntries2(m_pTChain->GetEntries(TCut(iter->first) + commonCut));
        }
        const double totalWeight(iter->second);
        if (totalWeight < std::numeric_limits<double>::epsilon())
        {
            m_evtTypeWeightMap.at(iter->first) = 1.f;
        }
        else
        {
            m_evtTypeWeightMap.at(iter->first) = totalWeight / static_cast<double>(m_pTChain->GetEntries(TCut(iter->first) + commonCut));
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

class Hist1DParameters
{
public:
    Hist1DParameters(const TString &varName);
    Hist1DParameters(const TString &varName, const float xLow, const float xHigh, const int xBin);
    Hist1DParameters(const TString &varName, const float xLow, const float xHigh, const int xBin, const float yLow, const float yHigh);
    Hist1DParameters(const TString &varName, const float xLow, const float xHigh, const int xBin, const float yLow, const float yHigh, 
        const TString &drawOption, const TString &cuts);
    Hist1DParameters(const TString &varName, const float xLow, const float xHigh, const int xBin, const float yLow, const float yHigh, 
        const TString &drawOption, const TString &cuts, const TString &xAxis, const TString &yAxis);
    float   m_xLow;
    float   m_xHigh;
    int   m_xNBins;
    float   m_yLow;
    float   m_yHigh;
    TString m_varName;
    TString m_drawOption;
    TString m_cuts;
    TString m_xAxis;
    TString m_yAxis;
};

//------------------------------------------------------------------------------------------------------------------------------------------

typedef std::vector<Hist1DParameters> Hist1DParametersVec;

//------------------------------------------------------------------------------------------------------------------------------------------

Hist1DParameters::Hist1DParameters(const TString &varName, const float xLow, const float xHigh, const int xBin, const float yLow, const float yHigh, 
    const TString &drawOption, const TString &cuts, const TString &xAxis, const TString &yAxis):
m_xLow(xLow),
m_xHigh(xHigh),
m_xNBins(xBin),
m_yLow(yLow),
m_yHigh(yHigh),
m_varName(varName),
m_drawOption(drawOption),
m_cuts(cuts),
m_xAxis(xAxis),
m_yAxis(yAxis)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

Hist1DParameters::Hist1DParameters(const TString &varName, const float xLow, const float xHigh, const int xBin, const float yLow, const float yHigh, 
    const TString &drawOption, const TString &cuts):
m_xLow(xLow),
m_xHigh(xHigh),
m_xNBins(xBin),
m_yLow(yLow),
m_yHigh(yHigh),
m_varName(varName),
m_drawOption(drawOption),
m_cuts(cuts),
m_xAxis(varName),
m_yAxis("Entries")
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

Hist1DParameters::Hist1DParameters(const TString &varName, const float xLow, const float xHigh, const int xBin, const float yLow, const float yHigh):
m_xLow(xLow),
m_xHigh(xHigh),
m_xNBins(xBin),
m_yLow(yLow),
m_yHigh(yHigh),
m_varName(varName),
m_drawOption(""),
m_cuts(""),
m_xAxis(varName),
m_yAxis("Entries")
{
}


//------------------------------------------------------------------------------------------------------------------------------------------

Hist1DParameters::Hist1DParameters(const TString &varName, const float xLow, const float xHigh, const int xBin):
m_xLow(xLow),
m_xHigh(xHigh),
m_xNBins(xBin),
m_yLow(0),
m_yHigh(0),
m_varName(varName),
m_drawOption(""),
m_cuts(""),
m_xAxis(varName),
m_yAxis("Entries")
{
}
//------------------------------------------------------------------------------------------------------------------------------------------

Hist1DParameters::Hist1DParameters(const TString &varName):
m_xLow(0),
m_xHigh(0),
m_xNBins(20),
m_yLow(0),
m_yHigh(0),
m_varName(varName),
m_drawOption(""),
m_cuts(""),
m_xAxis(varName),
m_yAxis("Entries")
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CLICStyle()
{
  gROOT->SetStyle("Plain"); /*Default white background for all plots*/
  /* set bkg color of all to kWhite: white, but not 0*/
  gStyle->SetCanvasColor(kWhite);
  gStyle->SetFrameFillColor(kWhite);
  gStyle->SetStatColor(kWhite);
  gStyle->SetPadColor(kWhite);
  gStyle->SetFillColor(10);
  gStyle->SetTitleFillColor(kWhite);
  
  
   /* SetPaperSize wants width & height in cm: A4 is 20,26 & US is 20,24*/
   gStyle->SetPaperSize(20, 26); 
   /* No yellow border around histogram*/
   gStyle->SetDrawBorder(0);
   /* remove border of canvas*/
   gStyle->SetCanvasBorderMode(0);
   /* remove border of pads*/
   gStyle->SetPadBorderMode(0);
   gStyle->SetFrameBorderMode(0);
   gStyle->SetLegendBorderSize(0);
  
   /* default text size*/
   gStyle->SetTextSize(0.05);
   gStyle->SetTitleSize(0.07,"xyz");
   gStyle->SetLabelSize(0.06,"xyz");
   /* title offset: distance between given text and axis, here x,y,z*/
   gStyle->SetLabelOffset(0.0015,"xyz");// default 0.0015
   gStyle->SetTitleOffset(1.1,"yz");//1.1
   gStyle->SetTitleOffset(1.0,"x");

   /* Use visible font for all text*/
   int font = 42; //132
   gStyle->SetTitleFont(font);
   gStyle->SetTitleFontSize(0.06);
   gStyle->SetStatFont(font);
   gStyle->SetStatFontSize(0.07);
   gStyle->SetTextFont(font);
   gStyle->SetLabelFont(font,"xyz");
   gStyle->SetTitleFont(font,"xyz");
   gStyle->SetTitleBorderSize(0);
   gStyle->SetStatBorderSize(1);

   /* big marker points*/
   gStyle->SetMarkerStyle(1);
   gStyle->SetLineWidth(2);  
   gStyle->SetMarkerSize(1.2);
   /*set palette in 2d histogram to nice and colorful one*/
   gStyle->SetPalette(1,0); 

   /*No title for histograms*/ //turn on title
   gStyle->SetOptTitle(1);
   /* show the errors on the stat box */
   gStyle->SetOptStat(0); 
   /* show errors on fitted parameters*/
   gStyle->SetOptFit(0); 
   /* number of decimals used for errors*/
   gStyle->SetEndErrorSize(5);   

   /* set line width to 2 by default so that histograms are visible when printed small
      idea: emphasize the data, not the frame around*/
   gStyle->SetHistLineWidth(2);
   gStyle->SetFrameLineWidth(2);
   gStyle->SetFuncWidth(2);
   gStyle->SetHistLineColor(kBlack);
   gStyle->SetFuncColor(kRed);
   gStyle->SetLabelColor(kBlack,"xyz");

   //set the margins
   gStyle->SetPadBottomMargin(0.18);
   gStyle->SetPadTopMargin(0.08);
   gStyle->SetPadRightMargin(0.05);//0.15
   gStyle->SetPadLeftMargin(0.17);// 0.17
   
   //set the number of divisions to show
   gStyle->SetNdivisions(506, "xy");
   
   //turn off xy grids
   gStyle->SetPadGridX(0);
   gStyle->SetPadGridY(0);
   
   //set the tick mark style
   gStyle->SetPadTickX(1);
   gStyle->SetPadTickY(1);

   gStyle->SetCanvasDefW(800);
   gStyle->SetCanvasDefH(700);

   gROOT->ForceStyle();
}

//------------------------------------------------------------------------------------------------------------------------------------------

TString safeName(const TString &name)
{
    TObject* old = gROOT->FindObject(name);
    if (old) delete old;
    return name;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void setStyle()
{
    (void) CLICStyle();
    gStyle->SetCanvasDefW(1600);
    gStyle->SetCanvasDefH(900);
    gStyle->SetTitleColor(TColor::GetColor( 33, 33, 33),"xyz");//1
    gStyle->SetLabelColor(TColor::GetColor( 97, 97, 97),"xyz");//1
    gStyle->SetGridColor(TColor::GetColor( 97, 97, 97));//1
    gStyle->SetAxisColor(TColor::GetColor( 97, 97, 97),"xyz");
    
    gStyle->SetPadRightMargin(0.05);//0.15 default for 2D plot
    gStyle->SetPadLeftMargin(0.17);// 0.17
    gStyle->SetTitleOffset(1.1);//1.1 
    gStyle->SetLabelOffset(0.02,"yz");//0.0015

}

//------------------------------------------------------------------------------------------------------------------------------------------

void setHistStyle(TH1 *&pHist, const Hist1DParameters &hist1Dparameters, const int counter)
{
    const int arrColour[] = {TColor::GetColor(114,147,203),TColor::GetColor(225,151,76),TColor::GetColor(131,186,91),
        TColor::GetColor(211,94,96),TColor::GetColor(128,133,133),TColor::GetColor(144,103,167),
        TColor::GetColor(171,104,87),TColor::GetColor(204,194,16)};

    const int arrColourFill[] = {TColor::GetColor(57,106,177),TColor::GetColor(218,124,48),TColor::GetColor(62,150,81),
        TColor::GetColor(204,37,41),TColor::GetColor(83,81,84),TColor::GetColor(107,76,154),
        TColor::GetColor(146,36,40),TColor::GetColor(148,139,61)};
        
    pHist->SetLineColor(arrColour[(counter % sizeof(arrColour))]);
    pHist->SetFillStyle(3001 + counter);
    
    pHist->SetFillColor(arrColour[(counter % sizeof(arrColour))]);
    pHist->SetTitle("");
    pHist->GetYaxis()->SetRangeUser(hist1Dparameters.m_yLow, hist1Dparameters.m_yHigh);
    pHist->GetXaxis()->SetTitle(hist1Dparameters.m_xAxis);
    pHist->GetYaxis()->SetTitle(hist1Dparameters.m_yAxis);
}

//------------------------------------------------------------------------------------------------------------------------------------------

TString getName(const TString &tstring)
{
    if (tstring.Contains("eventType == 1"))
    {
        return "e^{-}";
    }
    else if (tstring.Contains("eventType == 2"))
    {
        return "#mu^{-}";
    }
    else if (tstring.Contains("eventType == 3"))
    {
        return "#pi^{-}";
    }
    else if (tstring.Contains("eventType == 4"))
    {
        return "#pi^{-}#pi^{0}";
    }
    else if (tstring.Contains("eventType == 5"))
    {
        return "#pi^{-}#pi^{0}#pi^{0}";
    }
    else if (tstring.Contains("eventType == 6"))
    {
        return "#pi^{-}#pi^{+}#pi^{-}";
    }
    else if (tstring.Contains("eventType == 7"))
    {
        return "#pi^{-}#pi^{+}#pi^{-}#pi^{0}";
    }
    return tstring;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void saveGraph(TCanvas *&pCanvas, const TString &direcotry, const TString &fileName)
{
    if (!gSystem->OpenDirectory(direcotry))
        gSystem->MakeDirectory(direcotry);
    pCanvas->SaveAs(direcotry + fileName+".eps");
    pCanvas->SaveAs(direcotry + fileName+".C");
    pCanvas->SaveAs(direcotry + fileName+".pdf");
    pCanvas->SaveAs(direcotry + fileName+".png");
}

//------------------------------------------------------------------------------------------------------------------------------------------

int main(int argc, char* argv[]) {
    
    TApplication* pRootapp = new TApplication("example",&argc, argv);
    
    TStringDoubleMap evtTypeTotalWeightMap;
    evtTypeTotalWeightMap.insert(TStringDoubleMap::value_type("eventType == 1", 1.f));
    evtTypeTotalWeightMap.insert(TStringDoubleMap::value_type("eventType == 2", 1.f));
    evtTypeTotalWeightMap.insert(TStringDoubleMap::value_type("eventType == 3", 1.f));
    evtTypeTotalWeightMap.insert(TStringDoubleMap::value_type("eventType == 4", 1.f));
    evtTypeTotalWeightMap.insert(TStringDoubleMap::value_type("eventType == 5", 1.f));
    evtTypeTotalWeightMap.insert(TStringDoubleMap::value_type("eventType == 6", 1.f));
    evtTypeTotalWeightMap.insert(TStringDoubleMap::value_type("eventType == 7", 1.f));
    const bool debugFlag(false);
    const TCut commonCut("photonEC < 1 && mcCloseToZ < 1");
    ParametersVec parametersVec;
    parametersVec.push_back(Parameters("100GeV_improved", "/r06/lc/xu/TautauAnalysis/PandoraBono/100GeV/rootCustom/tauAnalysisTemplate_mod1_*.root", "sel", evtTypeTotalWeightMap, commonCut));
    //parametersVec.push_back(Parameters("200GeV_improved", "/r06/lc/xu/TautauAnalysis/PandoraBono/200GeV/rootCustom/tauAnalysisTemplate_mod1_*.root", "sel", evtTypeTotalWeightMap));
    //parametersVec.push_back(Parameters("500GeV_improved", "/r06/lc/xu/TautauAnalysis/PandoraBono/500GeV/rootCustom/tauAnalysisTemplate_mod1_*.root", "sel", evtTypeTotalWeightMap));
    //parametersVec.push_back(Parameters("1000GeV_improved", "/r06/lc/xu/TautauAnalysis/PandoraBono/1000GeV/rootCustom/tauAnalysisTemplate_mod1_*.root", "sel", evtTypeTotalWeightMap));
    parametersVec.push_back(Parameters("100GeV_old", "/r06/lc/xu/TautauAnalysis/PandoraDefault/100GeV/rootCustom/tauAnalysisTemplate_mod1_*.root", "sel", evtTypeTotalWeightMap, commonCut));
    //parametersVec.push_back(Parameters("200GeV_old", "/r06/lc/xu/TautauAnalysis/PandoraDefault/200GeV/rootCustom/tauAnalysisTemplate_mod1_*.root", "sel", evtTypeTotalWeightMap));
    //parametersVec.push_back(Parameters("500GeV_old", "/r06/lc/xu/TautauAnalysis/PandoraDefault/500GeV/rootCustom/tauAnalysisTemplate_mod1_*.root", "sel", evtTypeTotalWeightMap));
    //parametersVec.push_back(Parameters("1000GeV_old", "/r06/lc/xu/TautauAnalysis/PandoraDefault/1000GeV/rootCustom/tauAnalysisTemplate_mod1_*.root", "sel", evtTypeTotalWeightMap));
    

    Hist1DParametersVec hist1DParametersVec;
    //hist1DParametersVec.push_back(Hist1DParameters("nPhoton", -0.5, 6.5, 7, 0, 0, "", "", "Number of photons", "Normalised entries" ));
    //hist1DParametersVec.push_back(Hist1DParameters("nElectron", -0.5, 6.5, 7, 0, 1, "", "", "Number of e^{-}", "Normalised entries" ));
    hist1DParametersVec.push_back(Hist1DParameters("nMuon", -0.5, 6.5, 7, 0, 1, "", "", "Number of #mu^{-}", "Normalised entries" ));
    //hist1DParametersVec.push_back(Hist1DParameters("nCharge", -0.5, 6.5, 7, 0, 0, "", "", "Number of charged PFOs", "Normalised entries" ));
    //hist1DParametersVec.push_back(Hist1DParameters("nPionCharge", -0.5, 6.5, 7, 0, 0, "", "", "Number of #pi^{+/-}", "Normalised entries" ));
    //hist1DParametersVec.push_back(Hist1DParameters("mPionRhoFit", 0.01, 0.4, 100, 0, 0.2, "", "", "M_{#pi^{0}} Fit #rho", "Normalised entries" ));
    //hist1DParametersVec.push_back(Hist1DParameters("ePhoton", 0, 500, 100, 0, 1, "", "", "E_{#Sigma#gamma}", "Normalised entries" ));
    //hist1DParametersVec.push_back(Hist1DParameters("mRhoRhoFit", 0.01, 2, 100, 0, 0.2, "", "", "mass of #rho Fit #rho", "Normalised entries" ));
    //hist1DParametersVec.push_back(Hist1DParameters("mA1A1Fit", 0.01, 2, 100, 0, 0.2, "", "", "mass of a_{1} Fit a_{1}", "Normalised entries" ));
    //hist1DParametersVec.push_back(Hist1DParameters("mVis", 0.0, 2, 200, 0, 1, "", "", "mass of visibles", "Normalised entries" ));
    //hist1DParametersVec.push_back(Hist1DParameters("mPionChargeMod", 0.0, 2, 200, 0, 1, "", "", "M_{#pi^{+/-}} (mod)", "Normalised entries" ));
    //hist1DParametersVec.push_back(Hist1DParameters("mPionCharge", 0.0, 2, 200, 0, 1, "", "", "M_{#pi^{+/-}}", "Normalised entries" ));
    //hist1DParametersVec.push_back(Hist1DParameters("mChargeMod", 0.0, 2, 200, 0, 1, "", "", "mass of charged PFOs(mod)", "Normalised entries" ));
    //hist1DParametersVec.push_back(Hist1DParameters("mCharge", 0.0, 2, 200, 0, 1, "", "", "mass of charged PFOs", "Normalised entries" ));
    //hist1DParametersVec.push_back(Hist1DParameters("mPhoton", 0, 10));
    //hist1DParametersVec.push_back(Hist1DParameters("EEHCalRatio", 0.0, 1, 100, 0, 1, "", "", "E_{ECal}/E", "Normalised entries" ));
    
    //Hist1DParameters::Hist1DParameters(const TString &varName, const float xLow, const float xHigh, const int xBin, const float yLow, const float yHigh, 
    //const TString &drawOption, const TString &cuts, const TString &xAxis, const TString &yAxis):
    
    (void) setStyle();

    TCanvas* pCanvas = new TCanvas(safeName("canvas1"), "");
    pCanvas->UseCurrentStyle();//
    
    for (ParametersVec::const_iterator hIter = parametersVec.begin(), hIterEnd = parametersVec.end(); hIter != hIterEnd; ++hIter)
    {
        const Parameters &parameters(*hIter);
        
        for (Hist1DParametersVec::const_iterator iter = hist1DParametersVec.begin(), iterEnd = hist1DParametersVec.end(); iter != iterEnd; ++iter)
        {
            int counter(0);
            pCanvas->Clear();
            const Hist1DParameters &hist1DParameters(*iter);
            const TString padSize(TString::Format("(%i,%f,%f)", hist1DParameters.m_xNBins, hist1DParameters.m_xLow, hist1DParameters.m_xHigh));
            TLegend *legend=new TLegend(0.65,0.60,0.98,0.90);
            legend->SetTextSize(0.05);

            for (TStringDoubleMap::const_iterator jIter = parameters.m_evtTypeWeightMap.begin(), jIterEnd = parameters.m_evtTypeWeightMap.end(); jIter != jIterEnd; ++jIter)
            {
                const TString &evtTypeCut(jIter->first);
                const float weight(jIter->second);
                const TString cuts(TCut(hist1DParameters.m_cuts) + TCut(evtTypeCut) + parameters.m_commonCut);
                const TString weightAndcuts(TString::Format("%f*(", weight) + cuts + ")");
                const TString drawOption(hist1DParameters.m_drawOption + (parameters.m_evtTypeWeightMap.begin() != jIter ? "same" : "") );

                parameters.m_pTChain->Draw(hist1DParameters.m_varName + ">>"+ safeName(TString::Itoa(counter, 10)) + padSize, weightAndcuts, drawOption);
                TH1 *pTH1 = ((TH1*)gPad->GetPrimitive(TString::Itoa(counter, 10)));
                (void) setHistStyle(pTH1, hist1DParameters, counter);
                if (debugFlag)
                    std::cout << parameters.m_name << " " << hist1DParameters.m_xAxis << " " << getName(evtTypeCut) << " RMS " << pTH1->GetRMS() << std::endl;
                
                legend->AddEntry(pTH1, getName(evtTypeCut), "f");
                counter++;
                
            }
            legend-> SetFillStyle(0);
            legend-> SetBorderSize(0);
            legend->SetTextColor(TColor::GetColor( 97, 97, 97));
            legend->Draw("SAME");

            (void) saveGraph(pCanvas, "/var/clus/usera/xu/ILCSOFT/TauAnalysisNew/plots/", hist1DParameters.m_varName + "_" + parameters.m_name);
        }
    }
    pRootapp->Run();
    return 0;

}

