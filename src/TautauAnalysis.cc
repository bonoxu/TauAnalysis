#include "TautauAnalysis.h"
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <limits>

#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>

#include "UTIL/LCRelationNavigator.h"

#include "TLorentzVector.h"

#include "AlgorithmHelper.h"
#include "VarName.h"
#include "PdgTable.h"
#include "MCHelper.h"
#include "RecoHelper.h"

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

using namespace lcio;
using namespace marlin;

TautauAnalysis aTautauAnalysis;

TautauAnalysis::TautauAnalysis() : Processor("TautauAnalysis") 
{
    registerOptionalParameter( "PrintFlag" ,
                               "Flag for onscreen printing",
                               m_print,
                               false);
                               
    std::string tfileNameDefault("test.root");
    registerOptionalParameter( "TFileName" ,
                               "Name of Root Tree file",
                               m_tfileName,
                               tfileNameDefault);

    std::string ttreeNameDefault("sel");
    registerOptionalParameter( "TTreeName" ,
                               "Name of Tree",
                               m_ttreeName,
                               ttreeNameDefault);

    std::string ttreeIDDefault("sel");
    registerOptionalParameter( "TTreeID" ,
                               "ID of Tree",
                               m_ttreeID,
                               ttreeIDDefault);

    std::string pfoCollectionNameDefault("PandoraPfos");
    registerOptionalParameter( "PFOCollectionName" ,
                               "PFO Collection Name",
                               m_pfoCollectionName,
                               pfoCollectionNameDefault);

    registerOptionalParameter( "MCPFOCollectionName" ,
                               "MC PFO Collection Name",
                               m_mcPfoCollectionName,
                               std::string("MCParticle"));
                               
    registerOptionalParameter("RecoMCTruthLinkCollectionName" , 
				            "Name of reco particle to mc particle truth link relation"  ,
				            m_recoMCTruthLinkCollectionName,
				            std::string("RecoMCTruthLink"));
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void TautauAnalysis::init() 
{ 
    streamlog_out(DEBUG) << "   init called  " << std::endl;
    printParameters();

    m_nRun = 0;
    m_nEvent = 0;

    m_pTTreeHelper = new TTreeHelper(m_tfileName.c_str(), "RECREATE", m_ttreeName.c_str(), m_ttreeID.c_str());
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void TautauAnalysis::processRunHeader(LCRunHeader *run)
{ 
    ++m_nRun;
} 

//-----------------------------------------------------------------------------------------------------------------------------------------

void TautauAnalysis::processEvent(LCEvent *pEvent)
{ 
    if (m_nEvent % 500 == 0)
        streamlog_out(MESSAGE) << "Run " << m_nRun << ", TautauAnalysis::processEvent: " << m_nEvent << std::endl;
    streamlog_out(DEBUG) << "Run " << m_nRun << ", TautauAnalysis::processEvent: " << m_nEvent << std::endl;
    
    // DEBUG
    if(m_nEvent==0)
    {
        for (StringVec::const_iterator name = pEvent->getCollectionNames()->begin(), nameEnd = pEvent->getCollectionNames()->end(); name != nameEnd; ++name)
        {
            streamlog_out(MESSAGE) << "TautauAnalysis::processEvent Input Collections : " << *name << "  elements : " << pEvent->getCollection(*name)->getNumberOfElements() << std::endl;
        }
    }

    ++m_nEvent;
    
    LCCollection* pMCPfoCollection(pEvent->getCollection(m_mcPfoCollectionName));
    this->AnalyseMC(pMCPfoCollection);
    
    LCCollection* pPfoCollection(pEvent->getCollection(m_pfoCollectionName));

    ReconstructedParticleVec pfoVecMinus, pfoVecPlus;
    this->GetPfosInHemisphere(pPfoCollection, pfoVecMinus, pfoVecPlus);

    // DEBUG


    
    MCHelper::MCMCMap McToTauMCMap;
    MCHelper::AddAllDaughtersToMap(m_pMCTauMinus, McToTauMCMap);
    MCHelper::AddAllDaughtersToMap(m_pMCTauPlus, McToTauMCMap);
    
    /*for (MCHelper::MCMCMap::const_iterator iter = McToTauMCMap.begin(), iterEnd = McToTauMCMap.end(); iter != iterEnd; ++iter)
    {
        std::cout<< " McToTauMCMap "<< iter->first << " : " << iter->second << std::endl;
    }*/
    
    // DEBUG
    streamlog_out(DEBUG) << " nMCPFO in tau: " << McToTauMCMap.size() << std::endl;
    
    const LCCollection* pRecoMCTruthLinkCollection(pEvent->getCollection(m_recoMCTruthLinkCollectionName));
    const LCRelationNavigator* pRecoMCNavigator = new LCRelationNavigator(pRecoMCTruthLinkCollection);

    // DEBUG
    const MCParticle* pPlusMC(MCHelper::GetMCParticle(pfoVecPlus, pRecoMCNavigator, McToTauMCMap));
    const MCParticle* pMinusMC(MCHelper::GetMCParticle(pfoVecMinus, pRecoMCNavigator, McToTauMCMap));
    streamlog_out(DEBUG) << " +ve hemi is " << ((pPlusMC == m_pMCTauPlus) ? " plus mc " : " minus mc ") << " -ve hemi is " << ((pMinusMC == m_pMCTauPlus) ? " plus mc " : " minus mc ") << std::endl;
    
    // DEBUG
    streamlog_out(DEBUG) << " +ve hemi process nPFO in plus: " << pfoVecPlus.size() << " pdg: ";
    for (ReconstructedParticleVec::const_iterator iter = pfoVecPlus.begin(), iterEnd = pfoVecPlus.end(); iter != iterEnd; ++iter)
        streamlog_out(DEBUG) << (*iter)->getType() << "," ;
    streamlog_out(DEBUG)<< std::endl;
    
    this->AnalyseHemisphere(pfoVecPlus, pPfoCollection, pPlusMC, pMCPfoCollection);
    
    streamlog_out(DEBUG) << " -ve hemi process nPFO in minus: " << pfoVecMinus.size() << " pdg: "; 
    for (ReconstructedParticleVec::const_iterator iter = pfoVecMinus.begin(), iterEnd = pfoVecMinus.end(); iter != iterEnd; ++iter)
        streamlog_out(DEBUG) << (*iter)->getType() << ",";
    streamlog_out(DEBUG)<< std::endl;
    
    this->AnalyseHemisphere(pfoVecMinus, pPfoCollection, pMinusMC, pMCPfoCollection);
    //if (pPlusMC == pMinusMC)
    //    std::cin.get();
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void TautauAnalysis::check(LCEvent * pEvent)
{ 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void TautauAnalysis::end()
{
    // End of program - write out histograms and tree
   m_pTTreeHelper->WriteTTreeToTFile();
   delete m_pTTreeHelper;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void TautauAnalysis::AnalyseMC(LCCollection* pMCPfoCollection)
{
    const MCParticle *pMcEMinus(NULL), *pMcEPlus(NULL);
    this->GetMCInitialEEAfterISR(pMCPfoCollection, pMcEMinus, pMcEPlus);
    
    // DEBUG
    const MCHelper::MCIntMap &mcOrderMap(MCHelper::GetMCOrderMap(pMCPfoCollection));
    streamlog_out(DEBUG) << "E Minus id :" << mcOrderMap.at(pMcEMinus) <<std::endl;
    
    const MCParticle *pMcTauMinus(NULL), *pMcTauPlus(NULL);
    this->GetMCTauPairs(pMcEMinus, pMcTauMinus, pMcTauPlus);
    
    // DEBUG
    streamlog_out(DEBUG) << "Tau minus id :" << mcOrderMap.at(pMcTauMinus) << " plus :" << mcOrderMap.at(pMcTauPlus) <<std::endl;

    pMcTauMinus = this->GetMCAfterFSRPythia(pMcTauMinus);
    pMcTauPlus = this->GetMCAfterFSRPythia(pMcTauPlus);

    // DEBUG
    streamlog_out(DEBUG) << "New Tau minus id :" << mcOrderMap.at(pMcTauMinus) << " plus :" << mcOrderMap.at(pMcTauPlus) <<std::endl;

    m_pMCTauMinus = pMcTauMinus;
    m_pMCTauPlus  = pMcTauPlus;

    if (m_print)
        MCHelper::PrintMC(pMCPfoCollection, 50);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void TautauAnalysis::GetPfosInHemisphere(LCCollection* pPfoCollection, EVENT::ReconstructedParticleVec &pfoVecMinus, EVENT::ReconstructedParticleVec &pfoVecPlus) const
{
    FloatVec principleThrustAxis;
    pPfoCollection->parameters().getFloatVals("principleThrustAxis",principleThrustAxis);

    const int nParticles(pPfoCollection->getNumberOfElements());
    for (int i = 0; i < nParticles; ++i)
    {
        ReconstructedParticle *pReco(dynamic_cast<ReconstructedParticle*>(pPfoCollection->getElementAt(i)));
        const float dotProduct( principleThrustAxis[0] * pReco->getMomentum()[0] + principleThrustAxis[1] * pReco->getMomentum()[1] + 
                principleThrustAxis[2] * pReco->getMomentum()[2]);

        if (dotProduct > 0) 
        {
            pfoVecPlus.push_back(pReco);
        }
        else 
        {
            pfoVecMinus.push_back(pReco);
        }
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void TautauAnalysis::AnalyseHemisphere(const ReconstructedParticleVec &pfoVec, LCCollection* pPfoCollection, const MCParticle* pMainMC, LCCollection* pMCPfoCollection)
{
    m_pTTreeHelper->SetIntVar(VarName::GetName(VarName::N_EVENT), m_nEvent);
    
    const double thrustPrinciple(pPfoCollection->parameters().getFloatVal("principleThrustValue"));
    m_pTTreeHelper->SetDoubleVar(VarName::GetName(VarName::THRUST_PRINCIPLE), thrustPrinciple);
    
    streamlog_out(DEBUG) << " analyse MC " << std::endl;
    this->AnalyseHemisphereMC(pMainMC, pMCPfoCollection);
    
    //ReconstructedParticleVec preselectedPfoVec = this->GetPreselectedParticles(pfoVec);
    this->AnalyseHemisphereReco(pfoVec, pPfoCollection);
    
    m_pTTreeHelper->FillNode();
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void TautauAnalysis::AnalyseHemisphereMC(const MCParticle* pMainMC, LCCollection* pMCPfoCollection)
{
    int nE(0), nMu(0), nPionCharge(0), nPhotonFromPion(0), nOther(0), nNuE(0), nNuMu(0), nNuTau(0), nFSRPhoton(0);
    bool photonEarlyConversion(false), particleCloseToZ(false);
    if (pMainMC)
    {
        const EVENT::MCParticleVec mcDaughterVec(pMainMC->getDaughters());
        for (EVENT::MCParticleVec::const_iterator iter = mcDaughterVec.begin(), iterEnd = mcDaughterVec.end(); iter != iterEnd; ++iter)
        {
            const MCParticle* pMCDaughter(*iter);
            if (!particleCloseToZ && this->IsParticleInZAngle(pMCDaughter->getMomentum(), 0.2))
                particleCloseToZ = true;

            const int daughterPDG(std::fabs(pMCDaughter->getPDG()));
            switch (daughterPDG)
            {
                case NU_TAU:
                    nNuTau++;
                    break;
                case E_MINUS:
                    nE++;
                    break;
                case NU_E:
                    nNuE++;
                    break;
                case MU_MINUS:
                    nMu++;
                    break;
                case NU_MU:
                    nNuMu++;
                    break;
                case PI_PLUS:
                    nPionCharge++;
                    break;
                case PI_ZERO:
                {
                    this->AnalysePionNutralMC(pMCDaughter, nPhotonFromPion, nOther, photonEarlyConversion);
                    break;
                }
                case PHOTON:
                    nFSRPhoton++;
                    break;
                case RHO_770_PLUS:
                case A1_1260_PLUS:
                case W_PLUS:
                {
                    const EVENT::MCParticleVec mcDaughterDaughterVec(pMCDaughter->getDaughters());
                    for (EVENT::MCParticleVec::const_iterator jIter = mcDaughterDaughterVec.begin(), jIterEnd = mcDaughterDaughterVec.end(); jIter != jIterEnd; ++jIter)
                    {
                        const MCParticle* pMCDaughterDaughter(*jIter);
                        const int daughterDaughterPDG(std::fabs(pMCDaughterDaughter->getPDG()));
                        if (PI_PLUS == daughterDaughterPDG)
                        {
                            nPionCharge++;
                        }
                        else if (PI_ZERO == daughterDaughterPDG)
                        {
                            this->AnalysePionNutralMC(pMCDaughterDaughter, nPhotonFromPion, nOther, photonEarlyConversion);
                        }
                        else if (PHOTON == daughterDaughterPDG)
                        {
                            nFSRPhoton++;
                        }
                        else
                        {
                            streamlog_out(DEBUG) << " daughterDaughterPDG " << daughterDaughterPDG << std::endl;
                            nOther++;
                        }
                    }
                    break;
                }
                default:
                    streamlog_out(DEBUG) << "Other daughter pdg " << daughterPDG << std::endl;
                    nOther++;
                    break;
            }
        }
    }

    // DEBUG
    streamlog_out(DEBUG) << "nNuTau, nE, nNuE, nMu, nNuMu, nPionCharge, nPhotonFromPion, nOther, nFSRPhoton: " << nNuTau << "," << nE << "," << nNuE << "," << nMu << "," << nNuMu << "," <<
        nPionCharge << "," << nPhotonFromPion << "," << nOther << "," << nFSRPhoton << " early conversion " << photonEarlyConversion << std::endl;

    // Event selection
    const int nTotal(nNuTau + nE + nNuE + nMu + nNuMu + nPionCharge + nPhotonFromPion + nOther);
    TAU_DECAY_FINAL_STATES tauDecayFinalState(OTHER);
    if (1 == nNuTau && 1 == nE && 1 == nNuE && 3 == nTotal)
    {
        tauDecayFinalState = E;
    }
    else if (1 == nNuTau && 1 == nMu && 1 == nNuMu && 3 == nTotal)
    {
        tauDecayFinalState = MU;
    }
    else if (1 == nNuTau && 1 == nPionCharge && 2 == nTotal)
    {
        tauDecayFinalState = PION;
    }
    else if (1 == nNuTau && 1 == nPionCharge && 2 == nPhotonFromPion && 4 == nTotal)
    {
        tauDecayFinalState = PION2PHOTON;
    }
    else if (1 == nNuTau && 1 == nPionCharge && 4 == nPhotonFromPion && 6 == nTotal)
    {
        tauDecayFinalState = PION4PHOTON;
    }
    else if (1 == nNuTau && 3 == nPionCharge && 4 == nTotal)
    {
        tauDecayFinalState = PION2PION;
    }
    else if (1 == nNuTau && 3 == nPionCharge && 2 == nPhotonFromPion && 6 == nTotal)
    {
        tauDecayFinalState = PION2PION2PHOTON;
    }
    // DEBUG
    streamlog_out(DEBUG) << "Tau dacay final state " << tauDecayFinalState << std::endl;
    
    m_pTTreeHelper->SetIntVar(VarName::GetName(VarName::EVENT_TYPE), tauDecayFinalState);
    m_pTTreeHelper->SetIntVar(VarName::GetName(VarName::PHOTON_EARLY_CONVERSION), photonEarlyConversion);
    m_pTTreeHelper->SetIntVar(VarName::GetName(VarName::MC_CLOSE_Z), particleCloseToZ);
    m_pTTreeHelper->SetIntVar(VarName::GetName(VarName::PHOTON_FSR), nFSRPhoton);
    streamlog_out(DEBUG) << "particleCloseToZ " << particleCloseToZ << std::endl;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void TautauAnalysis::AnalyseHemisphereReco(const EVENT::ReconstructedParticleVec &pfoVec, LCCollection* pPfoCollection)
{
    m_pTTreeHelper->SetIntVar(VarName::GetName(VarName::N_PHOTON),  RecoHelper::GetNPfo(pfoVec, PHOTON));
    m_pTTreeHelper->SetIntVar(VarName::GetName(VarName::N_E),  RecoHelper::GetNPfo(pfoVec, E_MINUS));
    m_pTTreeHelper->SetIntVar(VarName::GetName(VarName::N_MU),  RecoHelper::GetNPfo(pfoVec, MU_MINUS));
    m_pTTreeHelper->SetIntVar(VarName::GetName(VarName::N_PIONCHARGE),  RecoHelper::GetNPfo(pfoVec, PI_PLUS));
    m_pTTreeHelper->SetIntVar(VarName::GetName(VarName::N_PFO),  pfoVec.size());
    
    int nCharge(0), nNeutral(0);
    double eECal(0.f), eHCal(0.f), eYoke(0.f), eLCal(0.f), eLHCal(0.f), eEHCalRatio(0.f), r0(0.f);
    TLorentzVector visMom(0.f, 0.f, 0.f, 0.f), neutralMom(0.f, 0.f, 0.f, 0.f);
    ReconstructedParticleVec photonVec;
    for (EVENT::ReconstructedParticleVec::const_iterator iter = pfoVec.begin(), iterEnd = pfoVec.end(); iter != iterEnd; ++iter)
    {
        const ReconstructedParticle * pReco(*iter);
        visMom += TLorentzVector(pReco->getMomentum()[0], pReco->getMomentum()[1], pReco->getMomentum()[2], pReco->getEnergy());
        if (PHOTON == pReco->getType())
            photonVec.push_back(const_cast<ReconstructedParticle*>(pReco));
        if (0 != std::fabs(pReco->getCharge()))
        {
            nCharge++;
            if (!pReco->getTracks().empty())
            {
                const EVENT::Track *pTrack(pReco->getTracks()[0]);
                const double d0(std::fabs(pTrack->getD0()));
                const double z0(std::fabs(pTrack->getZ0()));
                r0 = (std::sqrt(d0*d0 + z0*z0));
            }
            
            for ( EVENT::ClusterVec::const_iterator jIter = pReco->getClusters().begin(), jIterEnd = pReco->getClusters().end(); jIter != jIterEnd; ++jIter) 
            {
                const Cluster *pCluster(*jIter);
                eECal += pCluster->getSubdetectorEnergies()[0];
                eHCal += pCluster->getSubdetectorEnergies()[1];
                eYoke += pCluster->getSubdetectorEnergies()[2];
                eLCal += pCluster->getSubdetectorEnergies()[3];
                eLHCal += pCluster->getSubdetectorEnergies()[4];
            }
            eEHCalRatio = (eECal + eHCal > std::numeric_limits<double>::epsilon()) ? eECal / (eECal + eHCal) : std::numeric_limits<double>::max();
        }
        else
        {
            nNeutral++;
            neutralMom += TLorentzVector(pReco->getMomentum()[0], pReco->getMomentum()[1], pReco->getMomentum()[2], pReco->getEnergy());
        }
    }

    std::sort(photonVec.begin(), photonVec.end(), SortRecoParticleByEnergyDescendingOrder);
    TLorentzVector photonMom(0.f, 0.f, 0.f, 0.f);
    for (EVENT::ReconstructedParticleVec::const_iterator iter = photonVec.begin(), iterEnd = photonVec.end(); iter != iterEnd; ++iter)
    {
        const ReconstructedParticle * pReco(*iter);
        streamlog_out(DEBUG) << "photon E " << pReco->getEnergy() << std::endl;
        photonMom += TLorentzVector(pReco->getMomentum()[0], pReco->getMomentum()[1], pReco->getMomentum()[2], pReco->getEnergy());
    }
    
    m_pTTreeHelper->SetIntVar(VarName::GetName(VarName::N_CHARGE),  nCharge);
    m_pTTreeHelper->SetDoubleVar(VarName::GetName(VarName::E_EHCAL_RATIO),  eEHCalRatio);
    m_pTTreeHelper->SetDoubleVar(VarName::GetName(VarName::R0), r0);
    m_pTTreeHelper->SetDoubleVar(VarName::GetName(VarName::M_VIS), visMom.M());
    m_pTTreeHelper->SetDoubleVar(VarName::GetName(VarName::E_VIS), visMom.E());
    m_pTTreeHelper->SetDoubleVar(VarName::GetName(VarName::M_PHOTON), photonMom.M());
    m_pTTreeHelper->SetDoubleVar(VarName::GetName(VarName::E_PHOTON), photonMom.E());
    m_pTTreeHelper->SetDoubleVar(VarName::GetName(VarName::M_NEUTRAL), neutralMom.M());
    m_pTTreeHelper->SetIntVar(VarName::GetName(VarName::N_NEUTRAL),  nNeutral);
    m_pTTreeHelper->SetDoubleVar(VarName::GetName(VarName::E_NEUTRAL), neutralMom.E());

}

//-----------------------------------------------------------------------------------------------------------------------------------------

void TautauAnalysis::GetMCInitialEEAfterISR(const LCCollection* pMCPfoCollection, const MCParticle*& pMCEMinus, const MCParticle*& pMCEPlus) const
{
    // ATTN assume ILD mc structure, i.e. top down
    const int nParticles(pMCPfoCollection->getNumberOfElements());
    if (nParticles < 2)
    {
        streamlog_out(DEBUG) << " TautauAnalysis::GetInitialEEAfterISR less than 2 mc pfo!!" << std::endl;
        pMCEMinus = NULL;
        pMCEPlus = NULL;
        return;
    }

    const MCParticle* pMCparticle0(dynamic_cast<MCParticle*>(pMCPfoCollection->getElementAt(0)));
    const MCParticle* pMCparticle1(dynamic_cast<MCParticle*>(pMCPfoCollection->getElementAt(1)));
    if (E_MINUS == pMCparticle0->getPDG() && E_PLUS == pMCparticle1->getPDG())
    {
        pMCEMinus = pMCparticle0;
        pMCEPlus = pMCparticle1;
    }
    else if (E_MINUS == pMCparticle1->getPDG() && E_PLUS == pMCparticle0->getPDG())
    {
        pMCEMinus = pMCparticle1;
        pMCEPlus = pMCparticle0;
    }
    else
    {
        streamlog_out(DEBUG) << " TautauAnalysis::GetInitialEEAfterISR first two pfo are not electron positrons!!" << std::endl;
        pMCEMinus = NULL;
        pMCEPlus = NULL;
        return;
    }

    // Find ISR photons
    bool notISRphotons(false);
    const MCParticle *pMCEMinusTemp(pMCEMinus), *pMCEPlusTemp(pMCEPlus);
    while (!notISRphotons)
    {
        // ATTN only check e minus as e plus should be exactly same as e plus
        const EVENT::MCParticleVec mcEMinusDaughterVec(pMCEMinusTemp->getDaughters());
        for (EVENT::MCParticleVec::const_iterator iter = mcEMinusDaughterVec.begin(), iterEnd = mcEMinusDaughterVec.end(); iter != iterEnd; ++iter)
        {
            const MCParticle *pMCDaughter(*iter);
            const int daughterPdg(pMCDaughter->getPDG());
            if (E_MINUS == daughterPdg)
            {
                pMCEMinusTemp = pMCDaughter;
            }
            else if (E_PLUS == daughterPdg)
            {
                pMCEPlusTemp = pMCDaughter;
            }
            else if (PHOTON != std::fabs(daughterPdg))
            {
                notISRphotons = true;
            }
        }
        if (notISRphotons)
        {
            pMCEMinus = pMCEMinusTemp;
            pMCEPlus = pMCEPlusTemp;
        }
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void TautauAnalysis::GetMCTauPairs(const MCParticle* pMCEMinus, const MCParticle*& pMCTauMinus, const MCParticle*& pMCTauPlus) const
{
    if (!pMCEMinus)
        return;

    const EVENT::MCParticleVec mcEMinusDaughterVec(pMCEMinus->getDaughters());
    for (EVENT::MCParticleVec::const_iterator iter = mcEMinusDaughterVec.begin(), iterEnd = mcEMinusDaughterVec.end(); iter != iterEnd; ++iter)
    {
        const MCParticle *pMCDaughter(*iter);
        const int daughterPdg(pMCDaughter->getPDG());
        if (TAU_MINUS == daughterPdg)
        {
            pMCTauMinus = pMCDaughter;
        }
        else if (TAU_PLUS == daughterPdg)
        {
            pMCTauPlus = pMCDaughter;
        }
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------

const MCParticle* TautauAnalysis::GetMCAfterFSRPythia(const MCParticle* pMC) const
{
    if (!pMC)
        return NULL;

    // Find FSR photons and pythia process 92, 94
    bool notFSRorPythia(false);
    const MCParticle *pMCTemp(pMC);
    while (!notFSRorPythia)
    {
        const EVENT::MCParticleVec mcDaughterVec(pMCTemp->getDaughters());
        for (EVENT::MCParticleVec::const_iterator iter = mcDaughterVec.begin(), iterEnd = mcDaughterVec.end(); iter != iterEnd; ++iter)
        {
            const MCParticle *pMCDaughter(*iter);
            const int daughterPdg(pMCDaughter->getPDG());
            if (pMC->getPDG() == daughterPdg || 92 == std::fabs(daughterPdg) || 94 == std::fabs(daughterPdg))
            {
                pMCTemp = pMCDaughter;
            }
            else if ((PHOTON != std::fabs(daughterPdg)) && (92 != std::fabs(daughterPdg)) && (94 != std::fabs(daughterPdg)) && (pMC->getPDG() != daughterPdg) && 
                (TAU_MINUS != std::fabs(daughterPdg)))
            {
                notFSRorPythia = true;
            }
        }
    }
    return pMCTemp;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void TautauAnalysis::AnalysePionNutralMC(const MCParticle* pMCPion, int &nPhoton, int &nOther, bool &hasPhotonEarlyConversion) const
{
    if (!pMCPion)
        return;

    const EVENT::MCParticleVec mcDaughterVec(pMCPion->getDaughters());
    for (EVENT::MCParticleVec::const_iterator iter = mcDaughterVec.begin(), iterEnd = mcDaughterVec.end(); iter != iterEnd; ++iter)
    {
        const MCParticle* pMCPionDaughter(*iter);
        const int pionDaughterPDG(std::fabs(pMCPionDaughter->getPDG()));
        if (PHOTON == pionDaughterPDG)
        {
            nPhoton++;
            streamlog_out(DEBUG) << " TautauAnalysis::AnalysePionNutralMC photon energy"  << pMCPionDaughter->getEnergy() << std::endl;
            if (pMCPionDaughter->isDecayedInTracker())
            //if (!pMCPionDaughter->getDaughters().empty())
                hasPhotonEarlyConversion = true;
        }
        else
        {
            nOther++;
        }
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------

EVENT::ReconstructedParticleVec TautauAnalysis::GetPreselectedParticles(const EVENT::ReconstructedParticleVec &pfoVec) const
{
    ReconstructedParticleVec preselectedPfoVec;
    
    const ReconstructedParticle* pMostEReco(this->GetMostEChargedParticle(pfoVec));
    
    if (!pMostEReco) 
        return pfoVec;
    
    bool hasLepton(false);
    
    const float mostERecoMag2(pMostEReco->getMomentum()[0] * pMostEReco->getMomentum()[0] + pMostEReco->getMomentum()[1] * pMostEReco->getMomentum()[1] + 
                pMostEReco->getMomentum()[2] * pMostEReco->getMomentum()[2]);
    const float mostERecoMag(mostERecoMag2 < std::numeric_limits<float>::epsilon() ? std::numeric_limits<float>::epsilon() : std::sqrt(mostERecoMag2));
    
    for (EVENT::ReconstructedParticleVec::const_iterator iter = pfoVec.begin(), iterEnd = pfoVec.end(); iter != iterEnd; ++iter)
    {
        ReconstructedParticle *pReco(*iter);
        if (E_MINUS == std::fabs(pReco->getType()) || MU_MINUS == std::fabs(pReco->getType()))
        {
            hasLepton = true;
            break;
        }
        
        const float dotProduct( pMostEReco->getMomentum()[0] * pReco->getMomentum()[0] + pMostEReco->getMomentum()[1] * pReco->getMomentum()[1] + 
                pMostEReco->getMomentum()[2] * pReco->getMomentum()[2]);
        const float momMag2(pReco->getMomentum()[0] * pReco->getMomentum()[0] + pReco->getMomentum()[1] * pReco->getMomentum()[1] + 
                    pReco->getMomentum()[2] * pReco->getMomentum()[2]);
        const float momMag(momMag2 < std::numeric_limits<float>::epsilon() ? std::numeric_limits<float>::epsilon() : std::sqrt(momMag2)); 
        
        const float cosAngle(dotProduct / momMag / mostERecoMag);
        const bool bigAngleFlag(std::acos(std::fabs(cosAngle)) > 0.2 ? true : false);
        //streamlog_out(DEBUG) << "bigAngleFlag :" << bigAngleFlag << " angle :" << std::acos(std::fabs(cosAngle)) << " type:" <<pReco->getType()<<std::endl;

        if (bigAngleFlag) continue;
        
        preselectedPfoVec.push_back(pReco);
    }
    return (hasLepton) ? pfoVec : preselectedPfoVec;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

const ReconstructedParticle* TautauAnalysis::GetMostEChargedParticle(const EVENT::ReconstructedParticleVec &pfoVec) const
{
    const ReconstructedParticle* pMostEReco(NULL);
    float mostE(std::numeric_limits<float>::epsilon());
    for (EVENT::ReconstructedParticleVec::const_iterator iter = pfoVec.begin(), iterEnd = pfoVec.end(); iter != iterEnd; ++iter)
    {
        ReconstructedParticle *pReco(*iter);
        if (0 == pReco->getCharge()) continue;
        if (mostE < pReco->getEnergy())
        {
            mostE = pReco->getEnergy();
            pMostEReco = pReco;
        }
    }
    return pMostEReco;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

bool TautauAnalysis::IsParticleInZAngle(const double *momentum, const float minZAngle) const
{
    const float pt2(momentum[0] * momentum[0] + momentum[1] * momentum[1]);
    const float pt(pt2 < std::numeric_limits<float>::epsilon() ? std::numeric_limits<float>::epsilon() : std::sqrt(pt2));
    const float angle(std::atan(pt / std::fabs(momentum[2])));
    //streamlog_out(DEBUG) << " ploar angle " << angle << std::endl;
    return ((angle > minZAngle) ? false : true);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

bool TautauAnalysis::SortRecoParticleByEnergyDescendingOrder(const ReconstructedParticle * lhs, const ReconstructedParticle * rhs)
{
    return lhs->getEnergy() > rhs->getEnergy();
}

//-----------------------------------------------------------------------------------------------------------------------------------------
float TautauAnalysis::Chi2FitRho770(const EVENT::ReconstructedParticleVec &inputPfoVec, EVENT::ReconstructedParticleVec &outputPfoVec) const
{
    EVENT::ReconstructedParticleVec pionChargedVec, photonVec;
    for (EVENT::ReconstructedParticleVec::const_iterator iter = inputPfoVec.begin(), iterEnd = inputPfoVec.end(); iter != iterEnd; ++iter)
    {
        ReconstructedParticle *pReco(*iter);
        const int pdg(std::fabs(pReco->getType()));
        if (PI_PLUS == pdg)
        {
            pionChargedVec.push_back(pReco);
        }
        else if (PHOTON == pdg)
        {
            photonVec.push_back(pReco);
        }
    }

    std::size_t n = 5;
    std::size_t k = 3;

    std::vector<int> ints;
    for (int i = 0; i < n; ints.push_back(i++));

    do
    {
       for (int i = 0; i < k; ++i)
       {
          std::cout << ints[i];
       }
       std::cout << "\n";
    }
    while(AlgorithmHelper::next_combination(ints.begin(),ints.begin() + k,ints.end()));
    
    

    for (EVENT::ReconstructedParticleVec::const_iterator iter = pionChargedVec.begin(), iterEnd = pionChargedVec.end(); iter != iterEnd; ++iter)
    {
        ReconstructedParticle *pPionCharge(*iter);
        for (EVENT::ReconstructedParticleVec::const_iterator jIter = photonVec.begin(), jIterEnd = photonVec.end(); jIter != jIterEnd; ++jIter)
        {
            ReconstructedParticle *pPhotonLhs(*jIter);
            for (EVENT::ReconstructedParticleVec::const_iterator jIter = photonVec.begin(), jIterEnd = photonVec.end(); jIter != jIterEnd; ++jIter)
            {
                ReconstructedParticle *pPhotonLhs(*jIter);
            }
        }
    }
}
