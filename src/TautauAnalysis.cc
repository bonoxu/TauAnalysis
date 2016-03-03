#include "TautauAnalysis.h"
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <numeric>
#include <functional>
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
            if (!particleCloseToZ && RecoHelper::IsParticleInZAngle(pMCDaughter->getMomentum(), 0.2))
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
    double eECal(0.f), eHCal(0.f);
    ReconstructedParticleVec chargeVec, neutralVec;
    for (EVENT::ReconstructedParticleVec::const_iterator iter = pfoVec.begin(), iterEnd = pfoVec.end(); iter != iterEnd; ++iter)
    {
        ReconstructedParticle * pReco(*iter);
        if (0 != std::fabs(pReco->getCharge()))
        {
            chargeVec.push_back(pReco);
            for ( EVENT::ClusterVec::const_iterator jIter = pReco->getClusters().begin(), jIterEnd = pReco->getClusters().end(); jIter != jIterEnd; ++jIter) 
            {
                const Cluster *pCluster(*jIter);
                eECal += pCluster->getSubdetectorEnergies()[0];
                eHCal += pCluster->getSubdetectorEnergies()[1];
            }
        }
        else
        {
            neutralVec.push_back(pReco);
        }
    }
    const double eEHCalRatio((eECal + eHCal > std::numeric_limits<double>::epsilon()) ? eECal / (eECal + eHCal) : std::numeric_limits<double>::max());

    ReconstructedParticleVec photonVec(RecoHelper::GetPfoVec(pfoVec, PHOTON));

    // debug
    std::sort(photonVec.begin(), photonVec.end(), RecoHelper::SortRecoParticleByEnergyDescendingOrder);
    for (EVENT::ReconstructedParticleVec::const_iterator iter = photonVec.begin(), iterEnd = photonVec.end(); iter != iterEnd; ++iter)
    {
        streamlog_out(DEBUG) << "photon E " << (*iter)->getEnergy() << std::endl;
    }
    
    ReconstructedParticleVec muonVec(RecoHelper::GetPfoVec(pfoVec, MU_MINUS));
    ReconstructedParticleVec electronVec(RecoHelper::GetPfoVec(pfoVec, E_MINUS));
    ReconstructedParticleVec pionChargeVec(RecoHelper::GetPfoVec(pfoVec, PI_PLUS));
    
    const TLorentzVector visMom(RecoHelper::GetMomFromRecoParticleVec(pfoVec)), chargeMom(RecoHelper::GetMomFromRecoParticleVec(chargeVec)), 
        neutralMom(RecoHelper::GetMomFromRecoParticleVec(neutralVec)), photonMom(RecoHelper::GetMomFromRecoParticleVec(photonVec)), 
        electronMom(RecoHelper::GetMomFromRecoParticleVec(electronVec)), muonMom(RecoHelper::GetMomFromRecoParticleVec(muonVec)), 
        pionChargeMom(RecoHelper::GetMomFromRecoParticleVec(pionChargeVec));
    
    const int nPhoton(RecoHelper::GetNPfo(pfoVec, PHOTON)), nPionCharge(pionChargeVec.size());
    
    // rho test
    float rhoChi2(std::numeric_limits<float>::max());
    ReconstructedParticleVec rhoFitPionChargeVec, rhoFitPhotonVec;
    if (nPionCharge > 0  && nPhoton > 1)
    {
        rhoChi2 = this->Chi2FitRho770(pfoVec, rhoFitPionChargeVec, rhoFitPhotonVec);
    }
    float rhoChi2NegLog(-std::numeric_limits<float>::max()), rhoFitPionZeroM(0.f), rhoFitRhoM(0.f), rhoFitCosStarPhoton(0.f);
    if (!rhoFitPionChargeVec.empty() && (rhoFitPhotonVec.size() > 1))
    {
        rhoChi2NegLog = - std::log(rhoChi2);
        rhoFitPionZeroM = RecoHelper::GetMomFromRecoParticleVec(rhoFitPhotonVec).M();
        rhoFitRhoM = (RecoHelper::GetMomFromRecoParticleVec(rhoFitPhotonVec) + RecoHelper::GetMomFromRecoParticleVec(rhoFitPionChargeVec)).M();
        rhoFitCosStarPhoton = RecoHelper::GetCosineInCoMFrame(RecoHelper::GetMomFromRecoParticle(rhoFitPhotonVec[0]), RecoHelper::GetMomFromRecoParticle(rhoFitPhotonVec[1]));
    }

    // a1 test
    float a1Chi2(std::numeric_limits<float>::max());
    ReconstructedParticleVec a1FitPionChargeVec, a1FitPhotonVecLhs, a1FitPhotonVecRhs;
    if (nPionCharge > 0  && nPhoton > 3)
    {
        a1Chi2 = this->Chi2FitA11260(pfoVec, a1FitPionChargeVec, a1FitPhotonVecLhs, a1FitPhotonVecRhs);
    }
    float a1Chi2NegLog(-std::numeric_limits<float>::max()), a1FitPionZeroLhsM(0.f), a1FitPionZeroRhsM(0.f), a1FitA1M(0.f), a1FitCosStarLhs(0.f), a1FitCosStarRhs(0.f);
    if (!a1FitPionChargeVec.empty() && (a1FitPhotonVecLhs.size() > 1) && (a1FitPhotonVecRhs.size() > 1))
    {
        a1Chi2NegLog = - std::log(a1Chi2);
        a1FitPionZeroLhsM = RecoHelper::GetMomFromRecoParticleVec(a1FitPhotonVecLhs).M();
        a1FitPionZeroRhsM = RecoHelper::GetMomFromRecoParticleVec(a1FitPhotonVecRhs).M();
        a1FitA1M = (RecoHelper::GetMomFromRecoParticleVec(a1FitPhotonVecLhs)  + RecoHelper::GetMomFromRecoParticleVec(a1FitPhotonVecRhs) + 
            RecoHelper::GetMomFromRecoParticleVec(a1FitPionChargeVec)).M();
        a1FitCosStarLhs = RecoHelper::GetCosineInCoMFrame(RecoHelper::GetMomFromRecoParticle(a1FitPhotonVecLhs[0]), RecoHelper::GetMomFromRecoParticle(a1FitPhotonVecLhs[1]));
        a1FitCosStarRhs = RecoHelper::GetCosineInCoMFrame(RecoHelper::GetMomFromRecoParticle(a1FitPhotonVecRhs[0]), RecoHelper::GetMomFromRecoParticle(a1FitPhotonVecRhs[1]));
    }

    m_pTTreeHelper->SetDoubleVar(VarName::GetName(VarName::E_EHCAL_RATIO),  eEHCalRatio);

    m_pTTreeHelper->SetDoubleVar(VarName::GetName(VarName::M_VIS), visMom.M());
    m_pTTreeHelper->SetDoubleVar(VarName::GetName(VarName::E_VIS), visMom.E());
    m_pTTreeHelper->SetIntVar(VarName::GetName(VarName::N_PFO),  pfoVec.size());
    
    m_pTTreeHelper->SetIntVar(VarName::GetName(VarName::N_NEUTRAL),  neutralVec.size());
    m_pTTreeHelper->SetDoubleVar(VarName::GetName(VarName::M_NEUTRAL), neutralMom.M());
    m_pTTreeHelper->SetDoubleVar(VarName::GetName(VarName::E_NEUTRAL), neutralMom.E());
    
    m_pTTreeHelper->SetIntVar(VarName::GetName(VarName::N_CHARGE),  chargeVec.size());
    m_pTTreeHelper->SetIntVar(VarName::GetName(VarName::M_CHARGE),  chargeMom.M());
    m_pTTreeHelper->SetIntVar(VarName::GetName(VarName::E_CHARGE),  chargeMom.E());
    
    m_pTTreeHelper->SetIntVar(VarName::GetName(VarName::N_MU),  muonVec.size());
    m_pTTreeHelper->SetIntVar(VarName::GetName(VarName::E_MU),  muonMom.E());
    
    m_pTTreeHelper->SetIntVar(VarName::GetName(VarName::N_E),  electronVec.size());
    m_pTTreeHelper->SetIntVar(VarName::GetName(VarName::E_E),  electronMom.E());
    
    m_pTTreeHelper->SetIntVar(VarName::GetName(VarName::N_PHOTON), nPhoton);
    m_pTTreeHelper->SetDoubleVar(VarName::GetName(VarName::M_PHOTON), photonMom.M());
    m_pTTreeHelper->SetDoubleVar(VarName::GetName(VarName::E_PHOTON), photonMom.E());

    m_pTTreeHelper->SetIntVar(VarName::GetName(VarName::N_PIONCHARGE), nPionCharge);
    m_pTTreeHelper->SetIntVar(VarName::GetName(VarName::M_PIONCHARGE), pionChargeMom.M());
    m_pTTreeHelper->SetIntVar(VarName::GetName(VarName::E_PIONCHARGE), pionChargeMom.E());

    m_pTTreeHelper->SetDoubleVar(VarName::GetName(VarName::LOGCHI_RHOFIT), rhoChi2NegLog);
    m_pTTreeHelper->SetDoubleVar(VarName::GetName(VarName::M_PION_RHOFIT), rhoFitPionZeroM);
    m_pTTreeHelper->SetDoubleVar(VarName::GetName(VarName::M_RHO_RHOFIT), rhoFitRhoM);
    m_pTTreeHelper->SetDoubleVar(VarName::GetName(VarName::COSSTAR_RHOFIT), rhoFitCosStarPhoton);
    
    m_pTTreeHelper->SetDoubleVar(VarName::GetName(VarName::LOGCHI_A1FIT), a1Chi2NegLog);
    m_pTTreeHelper->SetDoubleVar(VarName::GetName(VarName::M_PION_LHS_A1FIT), a1FitPionZeroLhsM);
    m_pTTreeHelper->SetDoubleVar(VarName::GetName(VarName::M_PION_RHS_A1FIT), a1FitPionZeroRhsM);
    m_pTTreeHelper->SetDoubleVar(VarName::GetName(VarName::M_A1_A1FIT), a1FitA1M);
    m_pTTreeHelper->SetDoubleVar(VarName::GetName(VarName::COSSTAR_LHS_A1FIT), a1FitCosStarLhs);
    m_pTTreeHelper->SetDoubleVar(VarName::GetName(VarName::COSSTAR_RHS_A1FIT), a1FitCosStarRhs);

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

float TautauAnalysis::Chi2FitRho770(const EVENT::ReconstructedParticleVec &inputPfoVec, EVENT::ReconstructedParticleVec &pionChargePfoVec,
    EVENT::ReconstructedParticleVec &photonPfoVec) const
{
    EVENT::ReconstructedParticleVec pionChargedVec(RecoHelper::GetPfoVec(inputPfoVec, PI_PLUS));
    EVENT::ReconstructedParticleVec photonVec(RecoHelper::GetPfoVec(inputPfoVec, PHOTON));

    std::vector<int> pionChargeOrder;
    for (unsigned int iter = 0; iter < pionChargedVec.size(); pionChargeOrder.push_back(iter++));
    
    std::vector<int> photonOrder;
    for (unsigned int iter = 0; iter < photonVec.size(); photonOrder.push_back(iter++));

    const int nPionChargeFit(1), nPhotonFit(2);
    float bestChi2(std::numeric_limits<float>::max());
    do
    {
        EVENT::ReconstructedParticleVec pionChargedFitVec;
        for (int iter = 0; iter < nPionChargeFit; pionChargedFitVec.push_back(pionChargedVec[pionChargeOrder[iter++]]));
        
        do
        {
            EVENT::ReconstructedParticleVec photonFitVec;
            for (int iter = 0; iter < nPhotonFit; photonFitVec.push_back(photonVec[photonOrder[iter++]]));
            
            const TLorentzVector photonMom(RecoHelper::GetMomFromRecoParticleVec(photonFitVec));
            const TLorentzVector totalMom(RecoHelper::GetMomFromRecoParticleVec(pionChargedFitVec) + photonMom);
            
            const float pionChi2((photonMom.M() / PdgTable::GetParticleMass(PI_ZERO) - 1.f) * (photonMom.M() / PdgTable::GetParticleMass(PI_ZERO) - 1.f));
            const float rhoChi2((totalMom.M() / PdgTable::GetParticleMass(RHO_770_PLUS) - 1.f) * (totalMom.M() / PdgTable::GetParticleMass(RHO_770_PLUS) - 1.f));
            const float chi2(pionChi2 + rhoChi2);
            
            if (chi2 < bestChi2)
            {
                bestChi2 = chi2;
                photonPfoVec = photonFitVec;
                pionChargePfoVec = pionChargedFitVec;
            }
        }
        while(AlgorithmHelper::next_combination(photonOrder.begin(), photonOrder.begin() + nPhotonFit, photonOrder.end()));
    }
    while(AlgorithmHelper::next_combination(pionChargeOrder.begin(), pionChargeOrder.begin() + nPionChargeFit, pionChargeOrder.end()));
    
    
    // debug
    //streamlog_out(DEBUG) << "Chi2FitRho770 " << bestChi2  << " -log " << -std::log(bestChi2) << " rho m " << 
    //    (RecoHelper::GetMomFromRecoParticleVec(photonPfoVec) + RecoHelper::GetMomFromRecoParticleVec(pionChargePfoVec)).M() << 
    //    " pion m " << RecoHelper::GetMomFromRecoParticleVec(photonPfoVec).M() << 
    //    " cos* " << RecoHelper::GetCosineInCoMFrame(RecoHelper::GetMomFromRecoParticle(photonPfoVec[0]), RecoHelper::GetMomFromRecoParticle(photonPfoVec[1]))<< std::endl;
    //std::cin.get();
    return bestChi2;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

float TautauAnalysis::Chi2FitA11260(const EVENT::ReconstructedParticleVec &inputPfoVec, EVENT::ReconstructedParticleVec &pionChargePfoVec,
    EVENT::ReconstructedParticleVec &photonPfoVecLhs, EVENT::ReconstructedParticleVec &photonPfoVecRhs) const
{
    EVENT::ReconstructedParticleVec pionChargedVec(RecoHelper::GetPfoVec(inputPfoVec, PI_PLUS));
    EVENT::ReconstructedParticleVec photonVec(RecoHelper::GetPfoVec(inputPfoVec, PHOTON));

    IntVec pionChargeOrder;
    for (unsigned int iter = 0; iter < pionChargedVec.size(); pionChargeOrder.push_back(iter++));
    
    IntVec photonOrder;
    for (unsigned int iter = 0; iter < photonVec.size(); photonOrder.push_back(iter++));

    const int nPionChargeFit(1), nPhotonFit(4), nPionFit(2);
    float bestChi2(std::numeric_limits<float>::max());
    do
    {
        EVENT::ReconstructedParticleVec pionChargedFitVec;
        for (int iter = 0; iter < nPionChargeFit; pionChargedFitVec.push_back(pionChargedVec[pionChargeOrder[iter++]]));
        
        do
        {
            IntVec pionOrder;
            for (int iter = 0; iter < nPhotonFit; pionOrder.push_back(photonOrder[iter++]));
            
            do
            {
                EVENT::ReconstructedParticleVec pionFitVecLhs;
                for (int iter = 0; iter < nPionFit; pionFitVecLhs.push_back(photonVec[pionOrder[iter++]]));
                EVENT::ReconstructedParticleVec pionFitVecRhs;
                for (int iter = nPionFit; iter < nPhotonFit; pionFitVecRhs.push_back(photonVec[pionOrder[iter++]]));
                
                const TLorentzVector pionLhsMom(RecoHelper::GetMomFromRecoParticleVec(pionFitVecLhs));
                const TLorentzVector pionRhsMom(RecoHelper::GetMomFromRecoParticleVec(pionFitVecRhs));
                const TLorentzVector totalMom(RecoHelper::GetMomFromRecoParticleVec(pionChargedFitVec) + pionLhsMom + pionRhsMom);
                
                const float pionChi2Lhs((pionLhsMom.M() / PdgTable::GetParticleMass(PI_ZERO) - 1.f) * (pionLhsMom.M() / PdgTable::GetParticleMass(PI_ZERO) - 1.f));
                const float pionChi2Rhs((pionRhsMom.M() / PdgTable::GetParticleMass(PI_ZERO) - 1.f) * (pionRhsMom.M() / PdgTable::GetParticleMass(PI_ZERO) - 1.f));
                const float rhoChi2((totalMom.M() / PdgTable::GetParticleMass(A1_1260_PLUS) - 1.f) * (totalMom.M() / PdgTable::GetParticleMass(A1_1260_PLUS) - 1.f));
                const float chi2(pionChi2Lhs + pionChi2Rhs + rhoChi2);
                
                if (pionChi2Lhs < pionChi2Rhs && chi2 < bestChi2)
                {
                    bestChi2 = chi2;
                    photonPfoVecLhs = pionFitVecLhs;
                    photonPfoVecRhs = pionFitVecRhs;
                    pionChargePfoVec = pionChargedFitVec;
                }
            
            }
            while(AlgorithmHelper::next_combination(pionOrder.begin(), pionOrder.begin() + nPionFit, pionOrder.end()));

        }
        while(AlgorithmHelper::next_combination(photonOrder.begin(), photonOrder.begin() + nPhotonFit, photonOrder.end()));
    }
    while(AlgorithmHelper::next_combination(pionChargeOrder.begin(), pionChargeOrder.begin() + nPionChargeFit, pionChargeOrder.end()));
    
    
    // debug
    //streamlog_out(DEBUG) << "Chi2FitA1 " << bestChi2  << " -log " << -std::log(bestChi2) << " a1 m " << 
    //    (RecoHelper::GetMomFromRecoParticleVec(photonPfoVecLhs) + RecoHelper::GetMomFromRecoParticleVec(pionChargePfoVec) + RecoHelper::GetMomFromRecoParticleVec(photonPfoVecRhs)).M() << 
    //    " pion m " << RecoHelper::GetMomFromRecoParticleVec(photonPfoVecLhs).M() <<  " 2 " << RecoHelper::GetMomFromRecoParticleVec(photonPfoVecRhs).M() <<
    //    " cos* " << RecoHelper::GetCosineInCoMFrame(RecoHelper::GetMomFromRecoParticle(photonPfoVecLhs[0]), RecoHelper::GetMomFromRecoParticle(photonPfoVecLhs[1]))<< std::endl;
    //std::cin.get();
    return bestChi2;
}