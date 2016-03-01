#ifndef VarName_h
#define VarName_h 1

#include <string>

class VarName 
{

public:
    
    enum VARIABLE
    {
        N_EVENT,
        THRUST_PRINCIPLE,
        EVENT_TYPE,
        PHOTON_EARLY_CONVERSION,
        PHOTON_FSR,
        N_PHOTON,
        N_PIONCHARGE,
        N_MU,
        N_E,
        N_PFO,
        N_CHARGE,
        E_EHCAL_RATIO,
        R0,
        M_VIS,
        E_VIS,
        MC_CLOSE_Z,
        M_PHOTON,
        E_PHOTON,
        M_NEUTRAL,
        N_NEUTRAL,
        E_NEUTRAL,
        LOGCHI_RHOFIT,
        M_PION_RHOFIT,
        M_RHO_RHOFIT,
        COSSTAR_RHOFIT
    };
    
    static std::string GetName(VARIABLE var)
    {
        switch (var)
        {
            case N_EVENT:
                return "nEvt";
            case THRUST_PRINCIPLE:
                return "thrustPrinciple";
            case EVENT_TYPE:
                return "eventType";
            case PHOTON_EARLY_CONVERSION:
                return "photonEC";
            case PHOTON_FSR:
                return "photonFSR";
            case N_PHOTON:
                return "nPhoton";
            case N_PIONCHARGE:
                return "nPionCharge";
            case N_MU:
                return "nMuon";
            case N_E:
                return "nElectron";
            case N_PFO:
                return "nPfo";
            case N_CHARGE:
                return "nCharge";
            case E_EHCAL_RATIO:
                return "EEHCalRatio";
            case R0:
                return "R0";
            case M_VIS:
                return "mVis";
            case E_VIS:
                return "eVis";
            case MC_CLOSE_Z:
                return "mcCloseToZ";
            case M_PHOTON:
                return "mPhoton";
            case E_PHOTON:
                return "ePhoton";
            case M_NEUTRAL:
                return "mNeutral";
            case N_NEUTRAL:
                return "nNeutral";
            case E_NEUTRAL:
                return "eNeutral";
            case LOGCHI_RHOFIT:
                return "logChi2RhoFit";
            case M_PION_RHOFIT:
                return "mPionRhoFit";
            case M_RHO_RHOFIT:
                return "mRhoRhoFit";
            case COSSTAR_RHOFIT:
                return "cosStarRhoFit";
            default:
                return "unknown";
                break;
        }
    }
};

#endif
