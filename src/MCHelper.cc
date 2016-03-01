#include "MCHelper.h"
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <algorithm>

#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/LCObject.h>
#include <EVENT/ReconstructedParticle.h>
#include <UTIL/LCRelationNavigator.h>

using namespace lcio;

//-----------------------------------------------------------------------------------------------------------------------------------------

void MCHelper::PrintMC(const LCCollection* pMCPfoCollection, const int maxNFfo)
{
    const int nParticles(pMCPfoCollection->getNumberOfElements());
    const MCIntMap &mcOrderMap(GetMCOrderMap(pMCPfoCollection));
    for(int i = 0; i < std::min(nParticles, maxNFfo); ++i) 
    {
        const MCParticle* pMCparticle(dynamic_cast<MCParticle*>(pMCPfoCollection->getElementAt(i)));
        const EVENT::MCParticleVec mcDaughterVec(pMCparticle->getDaughters());
        const EVENT::MCParticleVec mcParentVec(pMCparticle->getParents());

        const double px(pMCparticle->getMomentum()[0]);
        const double py(pMCparticle->getMomentum()[1]);
        const double pz(pMCparticle->getMomentum()[2]);
        const double p2(px * px + py * py + pz * pz);
        const double p(p2 > std::numeric_limits<double>::epsilon() ? 0.f : std::sqrt(p2));
        // print
        std::cout << i << " " << pMCparticle << " : " << pMCparticle->getPDG() << " E:" << pMCparticle->getEnergy() << " P:" << p /*<< ":"<<px<<":"<<py<<":"<<pz*/ <<
            " M:" << pMCparticle->getMass();// << " nPa:" << _parents.size() << " nD:" << _daughters.size() ;
            
        std::cout << " PaPdg ";
        for (EVENT::MCParticleVec::const_iterator jIter = mcParentVec.begin(), jIterEnd = mcParentVec.end(); jIter != jIterEnd; ++jIter)
        {
            std::cout << (*jIter)->getPDG() << ":";
        }
        std::cout << " DPdg ";
        for (EVENT::MCParticleVec::const_iterator jIter = mcDaughterVec.begin(), jIterEnd = mcDaughterVec.end(); jIter != jIterEnd; ++jIter)
        {
            std::cout << (*jIter)->getPDG() << ":";
        }
        std::cout << " Pa ";
        // find parents ID
        for (EVENT::MCParticleVec::const_iterator jIter = mcParentVec.begin(), jIterEnd = mcParentVec.end(); jIter != jIterEnd; ++jIter)
        {
            if (mcOrderMap.find(*jIter) != mcOrderMap.end())
                std::cout <<mcOrderMap.at(*jIter)<< ":";
        }

        std::cout << " D ";
        // find daughter ID
        for (EVENT::MCParticleVec::const_iterator jIter = mcDaughterVec.begin(), jIterEnd = mcDaughterVec.end(); jIter != jIterEnd; ++jIter)
        {
            if (mcOrderMap.find(*jIter) != mcOrderMap.end())
                std::cout <<mcOrderMap.at(*jIter)<< ":";
        }
        
        std::cout<<" DecayedTracker,Calo,left,stop,status "<<pMCparticle->isDecayedInTracker()<<":"<<pMCparticle->isDecayedInCalorimeter()<<
            ":"<<pMCparticle->hasLeftDetector()<<":"<<pMCparticle->isStopped()<<":"<<pMCparticle->getGeneratorStatus();
        std::cout <<std::endl;
    }
    // take a break
    std::cin.get();
}

//-----------------------------------------------------------------------------------------------------------------------------------------

MCHelper::MCIntMap MCHelper::GetMCOrderMap(const LCCollection* pMCPfoCollection)
{
    MCIntMap mcOrderMap;
    const int nParticles(pMCPfoCollection->getNumberOfElements());
    for(int i = 0; i < nParticles; ++i) 
    {
        const MCParticle* pMCparticle(dynamic_cast<MCParticle*>(pMCPfoCollection->getElementAt(i)));
        if (!mcOrderMap.insert(MCIntMap::value_type(pMCparticle, i)).second)
            std::cout << "MCHelper::GetMCOrderMap fail to insert pair "<< pMCparticle << "," << i << std::endl;
    }
    return mcOrderMap;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void MCHelper::GetAllDaughters(const MCParticle *pMC, MCParticleSet &mcParticleDaughtersSet)
{
    if (!mcParticleDaughtersSet.insert(pMC).second)
        std::cout<<"MCHelper::GetAllDaughters failed to particle to set"<<std::endl;
        
    if (!pMC->getDaughters().empty())
    {
        for (EVENT::MCParticleVec::const_iterator iter = pMC->getDaughters().begin(), iterEnd = pMC->getDaughters().end(); iter != iterEnd; ++iter)
        {
            GetAllDaughters(*iter, mcParticleDaughtersSet);
        }
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void MCHelper::AddAllDaughtersToMap(const MCParticle *pMC, MCMCMap &daughterMcToMCMap)
{
    MCParticleSet mcParticleDaughtersSet;
    GetAllDaughters(pMC, mcParticleDaughtersSet);

    for (MCParticleSet::const_iterator iter = mcParticleDaughtersSet.begin(), iterEnd = mcParticleDaughtersSet.end(); iter != iterEnd; ++iter)
    {
        if (!daughterMcToMCMap.insert(MCMCMap::value_type(*iter, pMC)).second)
            std::cout << "MCHelper::AddAllDaughtersToMap failed to add particle to " << pMC << " map" << std::endl;
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------

const MCParticle* MCHelper::GetMCParticle(const ReconstructedParticle* pReco, const LCRelationNavigator* pRecoMCNavigator, const MCMCMap &daughterMcToMCMap)
{
    const LCObjectVec &mcParticlesVec(pRecoMCNavigator->getRelatedToObjects(const_cast<ReconstructedParticle*>(pReco)));
    if (mcParticlesVec.empty()) return NULL;
    //std::cout << "MCHelper::GetMCParticle mc raw " << mcParticlesVec.at(0) << " reco " << pReco << std::endl;
    
    const MCParticle *pMC(dynamic_cast<MCParticle*>(mcParticlesVec.at(0)));
    if (daughterMcToMCMap.find(pMC) != daughterMcToMCMap.end())
    {
        return daughterMcToMCMap.at(pMC);
    }
    else
    {
        return NULL;
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------

const MCParticle* MCHelper::GetMCParticle(const ReconstructedParticleVec &recoVec, const LCRelationNavigator* pRecoMCNavigator, const MCMCMap &daughterMcToMCMap)
{
    typedef std::map<const MCParticle*, float> MCFloatMap;
    MCFloatMap mcEnergyMap;
    for (ReconstructedParticleVec::const_iterator iter = recoVec.begin(), iterEnd = recoVec.end(); iter != iterEnd; ++iter)
    {
        const ReconstructedParticle* pReco(*iter);
        const float recoEnergy(pReco->getEnergy());
        const MCParticle* pMC(GetMCParticle(pReco, pRecoMCNavigator, daughterMcToMCMap));
        //std::cout << "MCHelper::GetMCParticle to insert mc " << pMC << " reco " << pReco << " type " << pReco->getType() << " e " << recoEnergy << std::endl;
        if (mcEnergyMap.find(pMC) != mcEnergyMap.end())
        {
            mcEnergyMap.at(pMC) += recoEnergy;
        }
        else
        {
            if (!mcEnergyMap.insert(MCFloatMap::value_type(pMC, recoEnergy)).second)
                std::cout << "MCHelper::GetMCParticle failed to insert mc " << pMC << " reco " << pReco << std::endl;
        }
    }

    const MCParticle* pBestMC(NULL);
    float highestE(0.f);
    for (MCFloatMap::const_iterator iter = mcEnergyMap.begin(), iterEnd = mcEnergyMap.end(); iter != iterEnd; ++iter)
    {
        if (iter->second > highestE && iter->first)
        {
            highestE = iter->second;
            pBestMC = iter->first;
        }
        //std::cout << "MCHelper::GetMCParticle mc " << iter->first << " e " << iter->second << std::endl;
    }
    //std::cout << "MCHelper::GetMCParticle best mc " << pBestMC << " e " << highestE << std::endl;
    return pBestMC;
}
