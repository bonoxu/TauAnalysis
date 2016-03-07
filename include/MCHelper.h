#ifndef MCHelper_h
#define MCHelper_h 1

#include <map>
#include <set>
#include <limits>
#include <typeinfo> 
 
#include "lcio.h"
#include "EVENT/MCParticle.h"
#include <EVENT/ReconstructedParticle.h>
using namespace lcio;
//namespace EVENT {class ReconstructedParticle;}
namespace UTIL {class LCRelationNavigator;}

class MCHelper
{

public:
    typedef std::map<const MCParticle*, int> MCIntMap;
    typedef std::set<const MCParticle*> MCParticleSet;
    typedef std::map<const MCParticle*, const MCParticle*> MCMCMap;

    /**
     *  @brief print mc particle from the collection
     * 
     *  @param pMCPfoCollection pointer to mc particle collection
     *  @param maxNPfo maximum number of pfo to display
     */
    static void PrintMC(const LCCollection* pMCPfoCollection, const int maxNFfo);

    /**
     *  @brief print mc particle from the collection
     * 
     *  @param pMCPfoCollection pointer to mc particle collection
     */
    static void PrintMC(const LCCollection* pMCPfoCollection);
    

    /**
     *  @brief get mc to order map, starting from 0
     * 
     *  @param pMCPfoCollection pointer to mc particle collection
     * 
     *  @return mc to order map
     */
    static MCIntMap GetMCOrderMap(const LCCollection* pMCPfoCollection);

    /**
     *  @brief  find daughters of a mc particle
     *
     *  @para	pMC pointer to MC particle
     *  @para	mcParticleDaughtersVec mc particle daughter vector to receive
     */
    static void GetAllDaughters(const MCParticle *pMC, MCParticleSet &mcParticleDaughtersSet);

    /**
     *  @brief  add all dauthers of a mc particle to the map
     *
     *  @para	pMC pointer to MC particle
     *  @para	daughterMcToParentMap daughter mc particles to mc paritcles of interests map to receive
     */
    static void AddAllDaughtersToMap(const MCParticle *pMC, MCMCMap &daughterMcToMCMap);
    
    /**
     *  @brief  find the main contributed mc particle of the mc particle of interests of a reco pfo
     *
     *  @para	pJet pointer to MC jet
     *  @para   pRecoMCNavigator pointer to reco to mc navigator
     *  @para	daughterMcToParentMap daughter mc particles to mc paritcles of interests map to receive
     * 
     */
    static const MCParticle* GetMCParticle(const ReconstructedParticle* pReco, const LCRelationNavigator* pRecoMCNavigator, const MCMCMap &daughterMcToMCMap);
    
    /**
     *  @brief  find the main contributed mc particle of the mc particle of interests of a reco pfo
     *
     *  @para	pJet pointer to MC jet
     *  @para   pRecoMCNavigator pointer to reco to mc navigator
     *  @para	daughterMcToParentMap daughter mc particles to mc paritcles of interests map to receive
     * 
     */
    static const MCParticle* GetMCParticle(const ReconstructedParticleVec &recoVec, const LCRelationNavigator* pRecoMCNavigator, const MCMCMap &daughterMcToMCMap);
private:
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline void MCHelper::PrintMC(const LCCollection* pMCPfoCollection)
{
    return PrintMC(pMCPfoCollection, std::numeric_limits<int>::max());
}

#endif



