#ifndef Physics_Analysis_KspipiAlg_H

#define Physics_Analysis_KspipiAlg_H
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/NTuple.h"
//#include "VertexFit/ReadBeamParFromDb.h"
#include "KspipiAlg/ReadBeamInfFromDb.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "VertexFit/KalmanKinematicFit.h"
#include "VertexFit/KinematicFit.h"

class KspipiAlg : public Algorithm
{

public:
    KspipiAlg(const std::string &name, ISvcLocator *pSvcLocator);
    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();

    typedef std::vector<int> Vint;
    typedef std::vector<double> Vdou;
    typedef std::vector<HepLorentzVector> Vp4;
    typedef std::vector<WTrackParameter> Vtrack;

private:
    double ecms;
    
    // define NTuple here

    NTuple::Tuple *m_tuple0;

    NTuple::Item<int> event_temp0;
    NTuple::Item<int> run_temp0;
    NTuple::Item<double> KsdecayL_temp0;
    NTuple::Item<double> KsdecayLerr_temp0;
    NTuple::Item<double> fKsdecayL_temp0;
    NTuple::Item<double> fKsdecayLerr_temp0;
    NTuple::Item<double> truthKs_x_temp0;
    NTuple::Item<double> truthKs_y_temp0;
    NTuple::Item<double> truthKs_z_temp0;
    NTuple::Item<double> truthvertex_x_temp0;
    NTuple::Item<double> truthvertex_y_temp0;
    NTuple::Item<double> truthvertex_z_temp0;
    NTuple::Item<double> truthdecayL_temp0;
    NTuple::Item<double> truthdecayL0_temp0;
    NTuple::Item<double> truthdecayR0_temp0;
    NTuple::Item<int> cutReason_temp0;
    NTuple::Item<int> fcutReason_temp0;
    NTuple::Item<int> ncharge_temp0;
    NTuple::Item<double> mKs_temp0;
    NTuple::Item<double> fmKs_temp0;
    NTuple::Item<double> m2pi_temp0;
    NTuple::Item<double> fm2pi_temp0;
    
    // mc vectors
    NTuple::Item<int> index_temp0;
    NTuple::Array<double> V4P_truth_Ks_temp0;
    NTuple::Array<double> V_P_pi_temp0;
    NTuple::Array<double> V_M_pi_temp0;
    NTuple::Array<double> V_Ks_temp0;
    NTuple::Array<double> fV_P_pi_temp0;
    NTuple::Array<double> fV_M_pi_temp0;
    NTuple::Array<double> fV_Ks_temp0;

    // decay chain
    NTuple::Item<int> idxmc_temp0;
    NTuple::Array<int> trkidx_temp0; 
    NTuple::Array<int> pdgid_temp0;
    NTuple::Array<int> motheridx_temp0;
};

#endif
