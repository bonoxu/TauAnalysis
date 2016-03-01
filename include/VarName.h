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
        MC_CLOSE_Z,
        M_PHOTON,
        M_NEUTRAL,
        N_NEUTRAL
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
            case MC_CLOSE_Z:
                return "mcCloseToZ";
            case M_PHOTON:
                return "mPhoton";
            case M_NEUTRAL:
                return "mNeutral";
            case N_NEUTRAL:
                return "nNeutral";
            default:
                return "unknown";
                break;
        }
    }
};

#endif
