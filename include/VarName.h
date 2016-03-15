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
        MC_CLOSE_Z,
        
        E_EHCAL_RATIO,

        N_PFO,
        M_VIS,
        E_VIS,

        N_CHARGE,
        M_CHARGE,
        E_CHARGE,
        
        M_CHARGE_MOD,
        E_CHARGE_MOD,
        
        N_NEUTRAL,
        M_NEUTRAL,
        E_NEUTRAL,

        N_MU,
        E_MU,

        N_E,
        E_E,

        N_PHOTON,
        M_PHOTON,
        E_PHOTON,

        N_PIONCHARGE,
        M_PIONCHARGE,
        E_PIONCHARGE,

        M_PIONCHARGE_MOD,
        E_PIONCHARGE_MOD,
        
        LOGCHI_RHOFIT,
        M_PION_RHOFIT,
        M_RHO_RHOFIT,
        COSSTAR_RHOFIT,
        
        LOGCHI_A1FIT,
        M_PION_LHS_A1FIT,
        M_PION_RHS_A1FIT,
        M_A1_A1FIT,
        COSSTAR_LHS_A1FIT,
        COSSTAR_RHS_A1FIT
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
            case MC_CLOSE_Z:
                return "mcCloseToZ";

            case E_EHCAL_RATIO:
                return "EEHCalRatio";

            case N_PFO:
                return "nPfo";
            case M_VIS:
                return "mVis";
            case E_VIS:
                return "eVis";
                
            case N_CHARGE:
                return "nCharge";
            case M_CHARGE:
                return "mCharge";
            case E_CHARGE:
                return "eCharge";

            case M_CHARGE_MOD:
                return "mChargeMod";
            case E_CHARGE_MOD:
                return "eChargeMod";

            case N_NEUTRAL:
                return "nNeutral";
            case M_NEUTRAL:
                return "mNeutral";
            case E_NEUTRAL:
                return "eNeutral";
                
            case N_MU:
                return "nMuon";
            case E_MU:
                return "eMuon";

            case N_E:
                return "nElectron";
            case E_E:
                return "eElectron";
                
            case N_PHOTON:
                return "nPhoton";
            case M_PHOTON:
                return "mPhoton";
            case E_PHOTON:
                return "ePhoton";
                
            case N_PIONCHARGE:
                return "nPionCharge";
            case M_PIONCHARGE:
                return "mPionCharge";
            case E_PIONCHARGE:
                return "ePionCharge";

            case M_PIONCHARGE_MOD:
                return "mPionChargeMod";
            case E_PIONCHARGE_MOD:
                return "ePionChargeMod";

            case LOGCHI_RHOFIT:
                return "logChi2RhoFit";
            case M_PION_RHOFIT:
                return "mPionRhoFit";
            case M_RHO_RHOFIT:
                return "mRhoRhoFit";
            case COSSTAR_RHOFIT:
                return "cosStarRhoFit";
                
            case LOGCHI_A1FIT:
                return "logChi2A1Fit";
            case M_PION_LHS_A1FIT:
                return "mPionLhsA1Fit";
            case M_PION_RHS_A1FIT:
                return "mPionRhsA1Fit";
            case M_A1_A1FIT:
                return "mA1A1Fit";
            case COSSTAR_LHS_A1FIT:
                return "cosStarLhsA1Fit";
            case COSSTAR_RHS_A1FIT:
                return "cosStarRhsA1Fit";
        
            default:
                return "unknown";
                break;
        }
    }
};

#endif
