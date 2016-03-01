#ifndef RecoHelper_h
#define RecoHelper_h 1

#include "TLorentzVector.h"

#include "lcio.h"
#include <EVENT/ReconstructedParticle.h>

using namespace lcio;

class RecoHelper
{

public:

    /**
     *  @brief Get n pfo for a given pdg code, for absolute pdgCode
     * 
     *  @param pfoVec pointer to pfo vector 
     *  @param pdgCode the pdg code of a particle
     *  
     *  @return Number of particle with that pdg code
     */
    static int GetNPfo(const EVENT::ReconstructedParticleVec &pfoVec, const int pdgCode);

    /**
     *  @brief Get n pfo for a given pdg code, if signFlag then sign matters
     * 
     *  @param pfoVec pointer to pfo vector 
     *  @param pdgCode the pdg code of a particle
     *  @param signFlag true if sign matters
     *  
     *  @return Number of particle with that pdg code
     */
    static int GetNPfo(const EVENT::ReconstructedParticleVec &pfoVec, const int pdgCode, const bool signFlag);
    
    /**
     *  @brief Get n pfo vec for a given pdg code, for absolute pdgCode
     * 
     *  @param pfoVec pointer to pfo vector 
     *  @param pdgCode the pdg code of a particle
     *  
     *  @return particle vec with that pdg code
     */
    static EVENT::ReconstructedParticleVec GetPfoVec(const EVENT::ReconstructedParticleVec &pfoVec, const int pdgCode);

    /**
     *  @brief Get n pfo vec for a given pdg code, if signFlag then sign matters
     * 
     *  @param pfoVec pointer to pfo vector 
     *  @param pdgCode the pdg code of a particle
     *  @param signFlag true if sign matters
     * 
     *  @return particle vec with that pdg code
     */
    static EVENT::ReconstructedParticleVec GetPfoVec(const EVENT::ReconstructedParticleVec &pfoVec, const int pdgCode, const bool signFlag);

    /**
     *  @brief True if particle polar angle is smaller than minZAngle
     * 
     *  @param momentum momentum
     *  @param minZAngle minimum accepted Z angle
     * 
     *  @return True if particle polar angle is smaller than minZAngle
     */
    static bool IsParticleInZAngle(const double *momentum, const float minZAngle);
    
    /**
     *  @brief Comparator sort reconstructed particle by decsending order
     * 
     *  @param lhs lhs reco particle
     *  @param rhs rhs reco particle
     * 
     *  @return True if lhs energy > rhs energy
     */
    static bool SortRecoParticleByEnergyDescendingOrder(const ReconstructedParticle *lhs, const ReconstructedParticle *rhs);
    
    /**
     *  @brief Get Momentum from reco particle
     * 
     *  @param pReco reco particle
     * 
     *  @return momentum vector
     */
    static TLorentzVector GetMomFromRecoParticle(const ReconstructedParticle *pReco);
    
    /**
     *  @brief Get Momentum from reco particle vec
     * 
     *  @param pfoVec reco particle vec
     * 
     *  @return momentum vector
     */
    static TLorentzVector GetMomFromRecoParticleVec(const EVENT::ReconstructedParticleVec &pfoVec);
    
    /**
     *  @brief get cosine of angle in com centre of mass frame
     * 
     *  @param lhs lhs reco particle
     *  @param rhs rhs reco particle
     * 
     *  @return cosine of angle in com centre of mass frame
     */
    static float GetCosineInCoMFrame(const ReconstructedParticle *lhs, const ReconstructedParticle *rhs);
    
    /**
     *  @brief get cosine of angle in com centre of mass frame
     * 
     *  @param lhs lhs reco particle momentum
     *  @param rhs rhs reco particle momentum
     * 
     *  @return cosine of angle in com centre of mass frame
     */
    static float GetCosineInCoMFrame(const TLorentzVector &lhs, const TLorentzVector &rhs);
    
private:

};

//------------------------------------------------------------------------------------------------------------------------------------------

inline int RecoHelper::GetNPfo(const EVENT::ReconstructedParticleVec &pfoVec, const int pdgCode)
{
    return RecoHelper::GetNPfo(pfoVec, pdgCode, false);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline int RecoHelper::GetNPfo(const EVENT::ReconstructedParticleVec &pfoVec, const int pdgCode, const bool signFlag)
{
    int nPfo(0);
    for (EVENT::ReconstructedParticleVec::const_iterator iter = pfoVec.begin(), iterEnd = pfoVec.end(); iter != iterEnd; ++iter)
    {
        if (signFlag ? pdgCode == (*iter)->getType() : std::fabs(pdgCode) == std::fabs((*iter)->getType()))
            nPfo++;
    }
    return nPfo;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline EVENT::ReconstructedParticleVec RecoHelper::GetPfoVec(const EVENT::ReconstructedParticleVec &pfoVec, const int pdgCode)
{
    return RecoHelper::GetPfoVec(pfoVec, pdgCode, false);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline EVENT::ReconstructedParticleVec RecoHelper::GetPfoVec(const EVENT::ReconstructedParticleVec &pfoVec, const int pdgCode, const bool signFlag)
{
    EVENT::ReconstructedParticleVec outputPfoVec;
    for (EVENT::ReconstructedParticleVec::const_iterator iter = pfoVec.begin(), iterEnd = pfoVec.end(); iter != iterEnd; ++iter)
    {
        if (signFlag ? pdgCode == (*iter)->getType() : std::fabs(pdgCode) == std::fabs((*iter)->getType()))
            outputPfoVec.push_back(*iter);
    }
    return outputPfoVec;
}


//-----------------------------------------------------------------------------------------------------------------------------------------

inline bool RecoHelper::IsParticleInZAngle(const double *momentum, const float minZAngle)
{
    const float pt2(momentum[0] * momentum[0] + momentum[1] * momentum[1]);
    const float pt(pt2 < std::numeric_limits<float>::epsilon() ? std::numeric_limits<float>::epsilon() : std::sqrt(pt2));
    const float angle(std::atan(pt / std::fabs(momentum[2])));
    return (angle < minZAngle);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

inline bool RecoHelper::SortRecoParticleByEnergyDescendingOrder(const ReconstructedParticle * lhs, const ReconstructedParticle * rhs)
{
    return lhs->getEnergy() > rhs->getEnergy();
}

//-----------------------------------------------------------------------------------------------------------------------------------------

inline TLorentzVector RecoHelper::GetMomFromRecoParticle(const ReconstructedParticle *pReco)
{
    return TLorentzVector(pReco->getMomentum()[0], pReco->getMomentum()[1], pReco->getMomentum()[2], pReco->getEnergy());
}

//-----------------------------------------------------------------------------------------------------------------------------------------

inline TLorentzVector RecoHelper::GetMomFromRecoParticleVec(const EVENT::ReconstructedParticleVec &pfoVec)
{
    TLorentzVector momentum(0.f, 0.f, 0.f, 0.f);
    for (EVENT::ReconstructedParticleVec::const_iterator iter = pfoVec.begin(), iterEnd = pfoVec.end(); iter != iterEnd; ++iter)
    {
        momentum += RecoHelper::GetMomFromRecoParticle(*iter);
    }
    return momentum;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

inline float RecoHelper::GetCosineInCoMFrame(const ReconstructedParticle *lhs, const ReconstructedParticle *rhs)
{
    return RecoHelper::GetCosineInCoMFrame(RecoHelper::GetMomFromRecoParticle(lhs), RecoHelper::GetMomFromRecoParticle(rhs));
}

//-----------------------------------------------------------------------------------------------------------------------------------------

inline float RecoHelper::GetCosineInCoMFrame(const TLorentzVector &lhs, const TLorentzVector &rhs)
{
    const TLorentzVector dijet(lhs + rhs);
    TLorentzVector decayParticleL(lhs);
    decayParticleL.Boost(-dijet.BoostVector());
    return (decayParticleL.Vect().Dot(dijet.Vect()) / (dijet.Vect().Mag() * decayParticleL.Vect().Mag()));
}

#endif
