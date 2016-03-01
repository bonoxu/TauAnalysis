#ifndef TautauAnalysis_h
#define TautauAnalysis_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
//#include <set>
#include <map>

#include "TTreeHelper.h"

using namespace lcio;
using namespace marlin;

namespace EVENT { class LCEvent; class MCParticle;}
#include <EVENT/ReconstructedParticle.h>
namespace UTIL {class LCRelationNavigator;}

class TautauAnalysis : public Processor 
{

public:

    
	virtual Processor*  newProcessor() 
	{
		return new TautauAnalysis; 
	}

	TautauAnalysis();

	virtual void init();

	virtual void processRunHeader(LCRunHeader *run);

	virtual void processEvent(LCEvent *pEvent); 

	virtual void check(LCEvent *pEvent); 

	virtual void end();

protected:

    enum TAU_DECAY_FINAL_STATES
    {
        E = 1,
        MU = 2,
        PION = 3,
        PION2PHOTON = 4,
        PION4PHOTON = 5,
        PION2PION = 6,
        PION2PION2PHOTON = 7,
        OTHER = 0
    };
    
    void AnalyseMC(LCCollection* pMCPfoCollection);
    
    /**
     *  @brief get pfo in hemisphere, based on thrust axis
     * 
     *  @param pPfoCollection the pointer to pfo collection
     *  @param pfoVecMinus pfo vector for negative hemisphere to receive
     *  @param pfoVecPlus pfo vector for positive hemisphere to receive
     */
    void GetPfosInHemisphere(LCCollection* pPfoCollection, EVENT::ReconstructedParticleVec &pfoVecMinus, EVENT::ReconstructedParticleVec &pfoVecPlus) const;
    
    /**
     *  @brief Analyse hemisphere. find the mc event type, fill the entries in the tree, write to the tree
     * 
     *  @param pfoVec pointer to pfo vector for hemisphere
     *  @param pPfoCollection pointer to colllection of all pfo 
     *  @param pMainMC pointer to main mc of the hemisphere
     *  @param pMCPfoCollection pointer to colllection of all mc particle 
     */
    void AnalyseHemisphere(const EVENT::ReconstructedParticleVec &pfoVec, LCCollection* pPfoCollection, const MCParticle* pMainMC, LCCollection* pMCPfoCollection);
    
    /**
     *  @brief Analyse hemisphere with mc information. find the mc event type, fill the entries regarding to mc in the tree
     * 
     *  @param pMainMC pointer to main mc of the hemisphere
     *  @param pMCPfoCollection pointer to colllection of all mc particle 
     */
    void AnalyseHemisphereMC(const MCParticle* pMainMC, LCCollection* pMCPfoCollection);
    
    /**
     *  @brief Analyse hemisphere with reco information. fill the entries in the tree
     * 
     *  @param pfoVec pointer to pfo vector for hemisphere
     *  @param pPfoCollection pointer to colllection of all pfo 
     */
    void AnalyseHemisphereReco(const EVENT::ReconstructedParticleVec &pfoVec, LCCollection* pPfoCollection);
    
    /**
     *  @brief get initial e plus and e mius after IST, assuming top down structure of mc particle
     * 
     *  @param pMCPfoCollection the pointer to mc pfo collection
     *  @param pMCEMinus the pointer to the e minus mc PFO to receive
     *  @param pMCEPlus the pointer to the e plus mc PFO to receive
     */
    void GetMCInitialEEAfterISR(const LCCollection* pMCPfoCollection, const MCParticle*& pMCEMinus, const MCParticle*& pMCEPlus) const;

    /**
     *  @brief get mc tau plus and tau minus
     * 
     *  @param pMCEMinus the pointer to mc e minus
     *  @param pMCTauMinus the pointer to the tau minus mc PFO to receive
     *  @param pMCTauPlus the pointer to the tau plus mc PFO to receive
     */
    void GetMCTauPairs(const MCParticle* pMCEMinus, const MCParticle*& pMCTauMinus, const MCParticle*& pMCTauPlus) const;

    /**
     *  @brief get mc particle after fsr photons and pythia process
     * 
     *  @param pMC the pointer to mc particle
     * 
     *  @return pointer to the mc particle after fsr photons and pythia process
     */
    const MCParticle* GetMCAfterFSRPythia(const MCParticle* pMC) const;
    
    /**
     *  @brief Analyse a mc neutral pion, to the number of photon daughters, number of other daughters and whether photon daughters have early conversion
     * 
     *  @param pMCPion the pointer to mc pion
     *  @param nPhoton the number of photon daughters to receive
     *  @param nOther the number of other daughters to receive
     *  @param hasPhotonEarlyConversion the pointer to mc pion
     */
    void AnalysePionNutralMC(const MCParticle* pMCPion, int &nPhoton, int &nOther, bool &hasPhotonEarlyConversion) const;

    /**
     *  @brief Chi squared fit for rho 770
     * 
     *  @param inputPfoVec input particles
     *  @param pionChargePfoVec output charged pions to receive
     *  @param photonPfoVec output pion to receive 
     *  @return Chi squared
     */
    float Chi2FitRho770(const EVENT::ReconstructedParticleVec &inputPfoVec, EVENT::ReconstructedParticleVec &pionChargePfoVec, EVENT::ReconstructedParticleVec &photonPfoVec) const;
    
    int                     m_nRun;                             // The counter for number of runs
    int                     m_nEvent;                           // The counter for number of events
    TTreeHelper             *m_pTTreeHelper;                    // The pointer to ttree helper
    bool                    m_print;                            // The flag for whether to print
    
    std::string             m_tfileName;                        // The tfile name
    std::string             m_ttreeName;                        // The ttree name
    std::string             m_ttreeID;                          // The ttree id
    std::string             m_pfoCollectionName;                // The pfo collection name
    std::string             m_mcPfoCollectionName;              // The mc pfo collection name
    std::string             m_recoMCTruthLinkCollectionName;    // The reco to mc truth linker name
    
    const MCParticle        *m_pMCTauMinus;                     // The pointer to mc tau minus
    const MCParticle        *m_pMCTauPlus;                      // The pointer to mc tau plus
};

#endif



