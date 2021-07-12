//author Liu Liang
//version 0.3
//add MCtruth four momentum of n nbar pi- pi+
//compare the sum of this four particles with J/psi
//check the kinematic fit
//set a new kinematicfit, use emc/2 -Epi as nbar input parameter
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/PropertyMgr.h"
#include "GaudiKernel/Bootstrap.h"
#include "GaudiKernel/INTupleSvc.h"
#include "GaudiKernel/NTuple.h"
#include "GaudiKernel/IHistogramSvc.h"

#include "VertexFit/IVertexDbSvc.h"
#include "VertexFit/KalmanKinematicFit.h"
#include "VertexFit/KinematicFit.h"
#include "VertexFit/VertexFit.h"
#include "VertexFit/Helix.h"
#include "VertexFit/SecondVertexFit.h"

#include "EventModel/EventModel.h"
#include "EventModel/Event.h"
#include "EventModel/EventHeader.h"

#include "EvtRecEvent/EvtRecEvent.h"
#include "EvtRecEvent/EvtRecTrack.h"

#include "DstEvent/TofHitStatus.h"

#include "McTruth/McParticle.h"
#include "Identifier/TofID.h"
#include "ParticleID/ParticleID.h"
#include "TrigEvent/TrigEvent.h"
#include "TrigEvent/TrigData.h"

#include "TMath.h"
#include "TF1.h"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Geometry/Point3D.h"

#include "TaglampimppAlg/lambda.h"
#include "TString.h"

#include "McDecayModeSvc/McDecayModeSvc.h"
#include "McDecayModeSvc/IMcDecayModeSvc.h"

using CLHEP::Hep2Vector;
using CLHEP::Hep3Vector;
using CLHEP::HepLorentzVector;

#ifndef ENABLE_BACKWARDS_COMPATIBITY
//#include "MCDataCorrectionAlg/mcdatacorrection.h"
typedef HepGeom::Point3D<double> HepPoint3D;
#endif

#include <vector>

const double mpi = 0.13957;
const double mkaon = 0.493677;
const double mproton = 0.938272;
const double mpion0 = 0.1349766;
const double mpion = 0.139570;
const double mEta = 0.547853;
const double n_mass = 0.939565;
const double nbar_mass = 0.939565;
const double Sigmam_mass = 1.197436;
const double Sigma_mass = 1.18937;
const double lambda_mass = 1.115683;
const double xmass[5] = {0.000511, 0.105658, 0.139570, 0.493677, 0.938272};
const double velc = 299.792458; // tof path unit in mm
typedef std::vector<int> Vint;
typedef std::vector<double> Vdouble;
typedef std::vector<HepLorentzVector> Vp4;

static int Ncut_none = 0, Ncut_ncharge = 0,
           Ncut_pid = 0, Ncut_nshower = 0,
           Ncut_npi0 = 0, Ncut_photon = 0;
static int Ncut_pppi0pi0 = 0, Ncut_ppi0 = 0,
           Ncut_pbarpi0 = 0, Ncut_nnbarpipi,
           Ncut_mc = 0;
static int Ncut_kin = 0, Ncut_nbar = 0,
           Ncut_asigma = 0, Ncut_kin_flag1 = 0,
           Ncut_kin_flag2 = 0;
static int Ncut_mass = 0, Ncut_length = 0,
           Ncut_chisq = 0, Liu = 0;

//------------------------------------------------------------------------------------------------
Taglampimpp::Taglampimpp(const std::string &name, ISvcLocator *pSvcLocator) : Algorithm(name, pSvcLocator)
{
    declareProperty("energyTreshold", m_energyTaglampimpp = 0.025);
    declareProperty("energyTreshold1", m_energyTaglampimpp1 = 0.025);
    declareProperty("energyTreshold2", m_energyTaglampimpp2 = 0.050);
    declareProperty("gammaAngleCut", m_gammaAngleCut = 20.0);
    declareProperty("rdmSeed", m_anglename = 2);
    declareProperty("cms", cms = 3.097);
}

//-------------------------------------------------------------------------------

//MCDataCorrection *nbarcorr;
//MCDataCorrection *nbarcorr1;

StatusCode Taglampimpp::initialize()
{
    MsgStream log(msgSvc(), name());

    log << MSG::INFO << "in initialize()" << endmsg;

    StatusCode status;

    NTuplePtr nt0(ntupleSvc(), "FILE1/lambdarec");
    if (nt0)
        m_tuple0 = nt0;
    else
    {
        m_tuple0 = ntupleSvc()->book("FILE1/lambdarec", CLID_ColumnWiseTuple, "ks N-Tuple example");
        if (m_tuple0)
        {
            status = m_tuple0->addItem("mclambdacos", m_mclambdacos);
            status = m_tuple0->addItem("mclambdabarcos", m_mclambdabarcos);
            status = m_tuple0->addItem("mclambdapx", m_mclambdapx);
            status = m_tuple0->addItem("mclambdapy", m_mclambdapy);
            status = m_tuple0->addItem("mclambdapz", m_mclambdapz);
            status = m_tuple0->addItem("mclambdae", m_mclambdae);
            status = m_tuple0->addItem("mclambdabarpx", m_mclambdabarpx);
            status = m_tuple0->addItem("mclambdabarpy", m_mclambdabarpy);
            status = m_tuple0->addItem("mclambdabarpz", m_mclambdabarpz);
            status = m_tuple0->addItem("mclambdabare", m_mclambdabare);
            status = m_tuple0->addItem("mcpxpim", m_mcpxpim);
            status = m_tuple0->addItem("mcpypim", m_mcpypim);
            status = m_tuple0->addItem("mcpzpim", m_mcpzpim);
            status = m_tuple0->addItem("mcepim", m_mcepim);
            status = m_tuple0->addItem("mcpxpip", m_mcpxpip);
            status = m_tuple0->addItem("mcpypip", m_mcpypip);
            status = m_tuple0->addItem("mcpzpip", m_mcpzpip);
            status = m_tuple0->addItem("mcepip", m_mcepip);
            status = m_tuple0->addItem("mcpxpp", m_mcpxpp);
            status = m_tuple0->addItem("mcpypp", m_mcpypp);
            status = m_tuple0->addItem("mcpzpp", m_mcpzpp);
            status = m_tuple0->addItem("mcepp", m_mcepp);
            status = m_tuple0->addItem("mcpxpm", m_mcpxpm);
            status = m_tuple0->addItem("mcpypm", m_mcpypm);
            status = m_tuple0->addItem("mcpzpm", m_mcpzpm);
            status = m_tuple0->addItem("mcepm", m_mcepm);

            status = m_tuple0->addItem("nmc_pp", m_nmc_pp);
            status = m_tuple0->addItem("nmc_pm", m_nmc_pm);
            status = m_tuple0->addItem("nmc_pip", m_nmc_pip);
            status = m_tuple0->addItem("nmc_pim", m_nmc_pim);
            status = m_tuple0->addItem("nmc_lambda", m_nmc_lambda);
            status = m_tuple0->addItem("nmc_lambdabar", m_nmc_lambdabar);

            status = m_tuple0->addItem("indexmc", m_idxmc, 0, 100);
            status = m_tuple0->addIndexedItem("trkidx", m_idxmc, m_trkidx);
            status = m_tuple0->addIndexedItem("pdgid", m_idxmc, m_pdgid);
            status = m_tuple0->addIndexedItem("motheridx", m_idxmc, m_motheridx);
            status = m_tuple0->addItem("nGood", m_nGood);
            status = m_tuple0->addItem("nCharge", m_nCharge);

            status = m_tuple0->addItem("runNo", m_runNo);
            status = m_tuple0->addItem("event", m_event);
            status = m_tuple0->addItem("Rvz0", m_Rvz0);
            status = m_tuple0->addItem("Rvxy0", m_Rvxy0);

            status = m_tuple0->addItem("lambda_decayL", m_lambda_decayL);
            status = m_tuple0->addItem("lambda_decayLerr", m_lambda_decayLerr);
            status = m_tuple0->addItem("lambda_chi2", m_lambda_chi2);
            status = m_tuple0->addItem("lambda_massv2", m_lambda_massv2);
            status = m_tuple0->addItem("lambda_pv2", m_lambda_pv2);

            status = m_tuple0->addItem("lambdabar_massv2", m_lambdabar_massv2);
            status = m_tuple0->addItem("lambdabar_pv2", m_lambdabar_pv2);

            status = m_tuple0->addItem("pip_p", m_pip_p);
            status = m_tuple0->addItem("pip_pt", m_pip_pt);
            status = m_tuple0->addItem("pip_cos", m_pip_cos);

            status = m_tuple0->addItem("pim_p", m_pim_p);
            status = m_tuple0->addItem("pim_pt", m_pim_pt);
            status = m_tuple0->addItem("pim_cos", m_pim_cos);

            status = m_tuple0->addItem("pm_p", m_pm_p);
            status = m_tuple0->addItem("pm_pt", m_pm_pt);
            status = m_tuple0->addItem("pm_cos", m_pm_cos);

            status = m_tuple0->addItem("pp_p", m_pp_p);
            status = m_tuple0->addItem("pp_pt", m_pp_pt);
            status = m_tuple0->addItem("pp_cos", m_pp_cos);

            status = m_tuple0->addItem("npp", m_npp);
            status = m_tuple0->addItem("npm", m_npm);
            status = m_tuple0->addItem("npip", m_npip);
            status = m_tuple0->addItem("npim", m_npim);
            status = m_tuple0->addItem("lambdabar_x", m_lambdabar_x);
            status = m_tuple0->addItem("lambdabar_y", m_lambdabar_y);
            status = m_tuple0->addItem("lambdabar_z", m_lambdabar_z);
            status = m_tuple0->addItem("lambda_x", m_lambda_x);
            status = m_tuple0->addItem("lambda_y", m_lambda_y);
            status = m_tuple0->addItem("lambda_z", m_lambda_z);
            status = m_tuple0->addItem("mcvtx_x", m_mcvtx_x);
            status = m_tuple0->addItem("mcvtx_y", m_mcvtx_y);
            status = m_tuple0->addItem("mcvtx_z", m_mcvtx_z);
            status = m_tuple0->addItem("mcdecayL", m_mcdecayL);
        }
        else
        {
            log << MSG::ERROR << "Cannot book N-Tuple:" << long(m_tuple0) << endmsg;
        }
    }

    /*
	NTuplePtr nt1(ntupleSvc(), "FILE1/FIT");
	if ( nt1 ) m_tuple1 = nt1;
	else	{
		m_tuple1 = ntupleSvc()->book ("FILE1/FIT", CLID_ColumnWiseTuple, "ks N-Tuple example");
		if( m_tuple1 )	{
				status  = m_tuple1->addItem ("indexmc", m_idxmc, 0, 100 );
				status  = m_tuple1->addIndexedItem ("trkidx", m_idxmc, m_trkidx);
				status  = m_tuple1->addIndexedItem ("pdgid", m_idxmc, m_pdgid);
				status  = m_tuple1->addIndexedItem ("motheridx", m_idxmc, m_motheridx);
				status  = m_tuple1->addItem ("nGood", m_nGood);
				status  = m_tuple1->addItem ("nCharge", m_nCharge);

  				status  = m_tuple1->addItem ("runNo", m_runNo);
				status  = m_tuple1->addItem ("event", m_event);
				status  = m_tuple1->addItem ("Rvz0", m_Rvz0);
				status  = m_tuple1->addItem ("Rvxy0", m_Rvxy0);



				status = m_tuple1->addItem ("lambda_decayL", m_lambda_decayL);
				status = m_tuple1->addItem ("lambda_decayLerr", m_lambda_decayLerr);
				status = m_tuple1->addItem ("lambda_chi2", m_lambda_chi2);
				status = m_tuple1->addItem ("lambda_massv2", m_lambda_massv2);
				status = m_tuple1->addItem ("lambda_pv2", m_lambda_pv2);

				status = m_tuple1->addItem ("lambdabar_massv2", m_lambdabar_massv2);
				status = m_tuple1->addItem ("lambdabar_pv2", m_lambdabar_pv2);



				status = m_tuple1->addItem ("pip_p", m_pip_p);
				status = m_tuple1->addItem ("pip_pt", m_pip_pt);
				status = m_tuple1->addItem ("pip_cos", m_pip_cos);

				status = m_tuple1->addItem ("pim_p", m_pim_p);
				status = m_tuple1->addItem ("pim_pt", m_pim_pt);
				status = m_tuple1->addItem ("pim_cos", m_pim_cos);

				status = m_tuple1->addItem ("pm_p", m_pm_p);
				status = m_tuple1->addItem ("pm_pt", m_pm_pt);
				status = m_tuple1->addItem ("pm_cos", m_pm_cos);

				status = m_tuple1->addItem ("pp_p", m_pp_p);
				status = m_tuple1->addItem ("pp_pt", m_pp_pt);
				status = m_tuple1->addItem ("pp_cos", m_pp_cos);


				status = m_tuple1->addItem ("npp", m_npp);
				status = m_tuple1->addItem ("npm", m_npm);
				status = m_tuple1->addItem ("npip", m_npip);
				status = m_tuple1->addItem ("npim", m_npim);
				status = m_tuple1->addItem ("lambdabar_x", m_lambdabar_x);
				status = m_tuple1->addItem ("lambdabar_y", m_lambdabar_y);
				status = m_tuple1->addItem ("lambdabar_z", m_lambdabar_z);
				status = m_tuple1->addItem ("lambda_x", m_lambda_x);
				status = m_tuple1->addItem ("lambda_y", m_lambda_y);
				status = m_tuple1->addItem ("lambda_z", m_lambda_z);
				status = m_tuple1->addItem ("mcvtx_x", m_mcvtx_x);
				status = m_tuple1->addItem ("mcvtx_y", m_mcvtx_y);
				status = m_tuple1->addItem ("mcvtx_z", m_mcvtx_z);
				status = m_tuple1->addItem ("mcdecayL", m_mcdecayL);
	

		}
		else {
				log << MSG::ERROR << "Connot book N-tuple : " << long (m_tuple1) << endmsg;
				return StatusCode::FAILURE;
		}
	}
*/
    cout << "++++++++++++++++++++++++++++++++++++++++ : lambda initial position" << endl;

    log << MSG::INFO << "successfully return from initialize()" << endmsg;
    return StatusCode::SUCCESS;
}

//========================================================================================

StatusCode Taglampimpp::execute()
{
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "in execute()" << endreq;
    SmartDataPtr<Event::EventHeader> eventHeader(eventSvc(), "/Event/EventHeader");
    
    int runNo = eventHeader->runNumber();
    int event = eventHeader->eventNumber();
    m_runNo = eventHeader->runNumber();
    m_event = eventHeader->eventNumber();
    log << MSG::DEBUG << "run, evtnum = "
        << runNo << " , "
        << event << endreq;
    Ncut_none++;
    //	cout << "**********Ncut_none :" << Ncut_none <<endl;

    int RadTrig = 0;
    if (runNo > 0)
    {
        SmartDataPtr<TrigData> trigData(eventSvc(), EventModel::Trig::TrigData);
        if (!trigData)
        {
            cout << "Could not find Trigger Data for physics analysis" << endl;
            return StatusCode::FAILURE;
        }
        RadTrig = trigData->getTrigChannel(9);
    }

    IMcDecayModeSvc *i_svc;
    StatusCode sc_DecayModeSvc = service("McDecayModeSvc", i_svc);
    if (sc_DecayModeSvc.isFailure())
    {
        log << MSG::FATAL << "Could not load McDecayModeSvc!" << endreq;
        return sc_DecayModeSvc;
    }
    m_svc = dynamic_cast<McDecayModeSvc *>(i_svc);

    if (runNo < 0)
    {
        SmartDataPtr<Event::McParticleCol> mcParticleCol(eventSvc(), "/Event/MC/McParticleCol");
        if (!mcParticleCol)
        {
            std::cout << "Could not retrieve McParticleCol" << std::endl;
            return StatusCode::FAILURE;
        }
        else
        {
            std::vector<int> pdgid;
            std::vector<int> motheridx;
            pdgid.clear();
            motheridx.clear();
            bool JpsiDecay = false;
            int rootIndex = -1;
            int m_numParticle = 0;
            Event::McParticleCol::iterator iter_mc = mcParticleCol->begin();
            for (; iter_mc != mcParticleCol->end(); iter_mc++)
            {
                if ((*iter_mc)->primaryParticle())
                    continue;
                if (!(*iter_mc)->decayFromGenerator())
                    continue;
                if ((*iter_mc)->particleProperty() == 443)
                {
                    int mode = m_svc->extract(*iter_mc, pdgid, motheridx);
                    m_numParticle = pdgid.size();
                    for (int i = 0; i < pdgid.size(); i++)
                    {
                        m_pdgid[i] = pdgid[i];
                        m_motheridx[i] = motheridx[i];
                    }
                    m_idxmc = m_numParticle;
                }
            }
        }
    }

    HepLorentzVector mc_p4pip, mc_p4pim, mc_p4pp, mc_p4pm, mc_p4lambda, mc_p4lambdabar;
    HepLorentzVector lambdabar_fposi, lambda_fposi;
    if (runNo < 0)
    {
        SmartDataPtr<Event::McParticleCol> mcParticleCol(eventSvc(), "/Event/MC/McParticleCol");
        if (mcParticleCol)
        {
            bool vphoDecay = false;
            int rootIndex = -1;
            int nmc_pip = 0, nmc_pim = 0, nmc_nbar = 0, nmc_pp = 0, nmc_lambda = 0, nmc_pi0 = 0, nmc_pm = 0, nmc_lambdabar = 0;
            Event::McParticleCol::iterator iter_mc = mcParticleCol->begin();
            for (; iter_mc != mcParticleCol->end(); iter_mc++)
            {
                if (!(*iter_mc)->decayFromGenerator())
                    continue;
                vphoDecay = true;
                HepLorentzVector mctrue_track = (*iter_mc)->initialFourMomentum();
                if ((*iter_mc)->particleProperty() == -3122 && ((*iter_mc)->mother()).particleProperty() == 443)
                { // anti-sigma+  -3112
                    mc_p4lambdabar = mctrue_track;
                    m_mclambdabarcos = cos(mc_p4lambdabar.theta());
                    m_mclambdabare = mc_p4lambdabar.e();
                    m_mclambdabarpx = mc_p4lambdabar.px();
                    m_mclambdabarpy = mc_p4lambdabar.py();
                    m_mclambdabarpz = mc_p4lambdabar.pz();
                    lambdabar_fposi = (*iter_mc)->finalPosition();

                    nmc_lambdabar++;
                }
                if ((*iter_mc)->particleProperty() == -2212 && ((*iter_mc)->mother()).particleProperty() == -3122)
                { //anti-n0 -2112
                    mc_p4pm = mctrue_track;
                    m_mcpxpm = mc_p4pm.px();
                    m_mcpypm = mc_p4pm.py();
                    m_mcpzpm = mc_p4pm.pz();
                    m_mcepm = mc_p4pm.e();
                    nmc_pm++;
                }
                if ((*iter_mc)->particleProperty() == 2212)
                {
                    mc_p4pp = mctrue_track;
                    m_mcpxpp = mc_p4pp.px();
                    m_mcpypp = mc_p4pp.py();
                    m_mcpzpp = mc_p4pp.pz();
                    m_mcepp = mc_p4pp.e();
                    nmc_pp++;
                }
                if ((*iter_mc)->particleProperty() == 211)
                {
                    mc_p4pip = mctrue_track;
                    m_mcpxpip = mc_p4pip.px();
                    m_mcpypip = mc_p4pip.py();
                    m_mcpzpip = mc_p4pip.pz();
                    m_mcepip = mc_p4pip.e();
                    nmc_pip++;
                }
                if ((*iter_mc)->particleProperty() == -211)
                {
                    mc_p4pim = mctrue_track;
                    m_mcpxpim = mc_p4pim.px();
                    m_mcpypim = mc_p4pim.py();
                    m_mcpzpim = mc_p4pim.pz();
                    m_mcepim = mc_p4pim.e();
                    nmc_pim++;
                }
                if ((*iter_mc)->particleProperty() == 3122 && ((*iter_mc)->mother()).particleProperty() == 443)
                {
                    mc_p4lambda = mctrue_track;
                    lambda_fposi = (*iter_mc)->finalPosition();
                    m_mclambdae = mctrue_track.e();
                    m_mclambdapx = mctrue_track.px();
                    m_mclambdapy = mctrue_track.py();
                    m_mclambdapz = mctrue_track.pz();
                    m_mclambdacos = cos(mctrue_track.theta());
                    nmc_lambda++;
                }
                //					if((*iter_mc)->particleProperty()==111 && ( ((*iter_mc)->mother()).particleProperty()==-3112 || ((*iter_mc)->mother()).particleProperty()==3112)) {
                //							nmc_pi0++;
                //					}
            }
            m_nmc_pp = nmc_pp;
            m_nmc_pm = nmc_pm;
            m_nmc_pip = nmc_pip;
            m_nmc_pim = nmc_pim;
            m_nmc_lambda = nmc_lambda;
            m_nmc_lambdabar = nmc_lambdabar;
            //					m_tuple0->write();
            Ncut_mc++;
        }
    }

    SmartDataPtr<EvtRecEvent> evtRecEvent(eventSvc(), EventModel::EvtRec::EvtRecEvent);
    log << MSG::DEBUG
        << " nCharge  = " << evtRecEvent->totalCharged() << " , "
        << " nNeutral = " << evtRecEvent->totalNeutral() << " , "
        << " tottrack = " << evtRecEvent->totalTracks() << endreq;
    /*	cout 	<< " nCharge  = " << evtRecEvent->totalCharged() << " , "
			<< " nneutral = " << evtRecEvent->totalNeutral() << " , "
			<< " tottrack = " << evtRecEvent->totalTracks()  << endl; */
    SmartDataPtr<EvtRecTrackCol> evtRecTrkCol(eventSvc(), EventModel::EvtRec::EvtRecTrackCol);
    //
    //	check x0, y0, z0, r0;
    //	suggest cut: |z0|<5 && r0<1
    //
    Vint iGood, ipip, ipim, iKp, iKm, ipp, ipm;
    iGood.clear();
    ipip.clear();
    ipim.clear();
    iKp.clear();
    iKm.clear();
    ipp.clear();
    ipm.clear();

    Vp4 ppip, ppim;
    ppip.clear();
    ppim.clear();

    int nCharge = 0;

    Hep3Vector xorigin(0, 0, 0);
    HepPoint3D vx(0., 0., 0.);
    HepSymMatrix Evx(3, 0);
    //if (m_reader.isRunNumberValid(runNo))
    IVertexDbSvc *vtxsvc;
    Gaudi::svcLocator()->service("VertexDbSvc", vtxsvc);
    if (vtxsvc->isVertexValid())
    {
        double *dbv = vtxsvc->PrimaryVertex();
        double *vv = vtxsvc->SigmaPrimaryVertex();
        /*		cout << "************Vertex is:" << dbv[0] <<"	  "<<dbv[1] << "    "<<dbv[2] <<"   "<< endl; */
        xorigin.setX(dbv[0]);
        xorigin.setY(dbv[1]);
        xorigin.setZ(dbv[2]);
        vx.setX(dbv[0]);
        vx.setY(dbv[1]);
        vx.setZ(dbv[2]);
        Evx[0][0] = vv[0] * vv[0];
        Evx[1][1] = vv[1] * vv[1];
        Evx[2][2] = vv[2] * vv[2];
    }
    VertexParameter vx_db;
    vx_db.setVx(vx);
    vx_db.setEvx(Evx);

    for (int i = 0; i < evtRecEvent->totalCharged(); i++)
    {
        EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + i;
        if (!(*itTrk)->isMdcTrackValid())
            continue;
        RecMdcTrack *mdcTrk = (*itTrk)->mdcTrack();
        double pch = mdcTrk->p();
        double x0 = mdcTrk->x();
        double y0 = mdcTrk->y();
        double z0 = mdcTrk->z();
        double phi0 = mdcTrk->helix(1);
        double xv = xorigin.x();
        double yv = xorigin.y();
        double Rzy = (x0 - xv) * cos(phi0) + (y0 - yv) * sin(phi0);

        HepVector a = mdcTrk->helix();
        HepSymMatrix Ea = mdcTrk->err();
        HepPoint3D point0(0., 0., 0.); // the initial point for MDC recosntruction
        HepPoint3D IP(xorigin[0], xorigin[1], xorigin[2]);
        VFHelix helixip(point0, a, Ea);
        helixip.pivot(IP);
        HepVector vecipa = helixip.a();
        double Rvxy0 = fabs(vecipa[0]); // the nearest distance to IP in xy plane
        double Rvz0 = vecipa[3];
        double Rvphi0 = vecipa[1];
        double cost = cos(mdcTrk->theta());

        m_Rvxy0 = Rvxy0;
        m_Rvz0 = Rvz0;

        //		if(pch > 2.0) continue;

        //	if(fabs(Rvz0) >= 30.0 ) continue;
        //	if(fabs(Rvxy0) >= 10.0 ) continue;
        if (fabs(cost) >= 0.93)
        {
            cout << "without kaltrack ,runNo=" << m_runNo << " eventNo=" << m_event << endl;
            continue;
        }
        iGood.push_back(i);
        RecMdcKalTrack *mdcKalTrk = (*itTrk)->mdcKalTrack();

        if (!(*itTrk)->isMdcKalTrackValid())
        {
            //	cout<<"without kaltrack ,runNo="<<m_runNo<<" eventNo="<<m_event<<endl;
            continue;
        }
        if (mdcTrk->charge() > 0 && mdcTrk->p() > 0.5)
        {
            ipp.push_back(i);
            m_pp_pt = mdcTrk->pxy();
            m_pp_p = mdcTrk->p();
            m_pp_cos = cos(mdcTrk->theta());
        }
        if (mdcTrk->charge() > 0 && mdcTrk->p() < 0.5)
        {
            ipip.push_back(i);
            m_pip_pt = mdcTrk->pxy();
            m_pip_p = mdcTrk->p();
            m_pip_cos = cos(mdcTrk->theta());
        }
        if (mdcTrk->charge() < 0 && mdcTrk->p() > 0.5)
        {
            ipm.push_back(i);
            m_pm_pt = mdcTrk->pxy();
            m_pm_p = mdcTrk->p();
            m_pm_cos = cos(mdcTrk->theta());
        }
        if (mdcTrk->charge() < 0 && mdcTrk->p() < 0.5)
        {
            ipim.push_back(i);
            m_pim_pt = mdcTrk->pxy();
            m_pim_p = mdcTrk->p();
            m_pim_cos = cos(mdcTrk->theta());
        }
        nCharge += mdcTrk->charge();
    }
    int nGood = iGood.size();

    m_nGood = nGood;
    m_nCharge = nCharge;

    int npip = ipip.size();
    int npim = ipim.size();
    int npp = ipp.size();
    int npm = ipm.size();
    m_npip = npip;
    m_npim = npim;
    m_npp = ipp.size();
    m_npm = ipm.size();

    // add MC vertex infomation --Long LI
    m_mcvtx_x = vx.x();
    m_mcvtx_y = vx.y();
    m_mcvtx_z = vx.z();

    //	final position of MC truth (lambda &lambdabar)
    m_lambdabar_x = lambdabar_fposi.x();
    m_lambdabar_y = lambdabar_fposi.y();
    m_lambdabar_z = lambdabar_fposi.z();
    m_lambda_x = lambda_fposi.x();
    m_lambda_y = lambda_fposi.y();
    m_lambda_z = lambda_fposi.z();

    //	MC decay length
    m_mcdecayL = sqrt((lambda_fposi.x() - vx.x()) * (lambda_fposi.x() - vx.x()) + (lambda_fposi.y() - vx.y()) * (lambda_fposi.y() - vx.y()) + (lambda_fposi.z() - vx.z()) * (lambda_fposi.z() - vx.z()));

    //	if(nGood < 2) return StatusCode::SUCCESS;
    Ncut_ncharge++;

    //	if(npp != 1) return SUCCESS;
    //	if(npim != 1) return SUCCESS;
    Ncut_pid++;
    if (nGood >= 2 && npp == 1 && npim == 1)
    {
        /////////////////////////////////////////////////////////////////////
        //=======================Kalman Kinematic Fit=========================

        HepLorentzVector ecms;
        ecms.setPy(0);
        ecms.setPz(0);
        ecms.setE(cms);
        ecms.setPx(cms * sin(0.022 * 0.5));
        //	cout << "===================================ipip : " << ipip[0] <<"  ipim : " << ipim[0] << endl;
        RecMdcKalTrack *pimTrk = (*(evtRecTrkCol->begin() + ipim[0]))->mdcKalTrack();
        RecMdcKalTrack *ppTrk = (*(evtRecTrkCol->begin() + ipp[0]))->mdcKalTrack();

        WTrackParameter wvppTrk = WTrackParameter(mproton, ppTrk->getZHelixP(), ppTrk->getZErrorP());
        HepLorentzVector p4pp = wvppTrk.p();

        WTrackParameter wvpimTrk = WTrackParameter(mpion, pimTrk->getZHelix(), pimTrk->getZError());
        HepLorentzVector p4pim = wvpimTrk.p();

        HepPoint3D vx_vtx(0., 0., 0.);
        HepSymMatrix Evx_vtx(3, 0);
        double bx = 1E+6;
        double by = 1E+6;
        double bz = 1E+6;
        Evx_vtx[0][0] = bx * bx;
        Evx_vtx[1][1] = by * by;
        Evx_vtx[2][2] = bz * bz;

        double chisq_vertex1 = 999.;
        bool good_lambda = false;
        double lambda_chi2 = 999., lambda_massv2 = 999.;
        WTrackParameter wlambda_vertex;
        VertexParameter vxpar;
        vxpar.setVx(vx_vtx);
        vxpar.setEvx(Evx_vtx);
        VertexFit *vtxfit = VertexFit::instance();
        vtxfit->init();
        vtxfit->AddTrack(0, wvppTrk);
        vtxfit->AddTrack(1, wvpimTrk);
        vtxfit->AddVertex(0, vxpar, 0, 1);
        if (vtxfit->Fit(0))
        {
            vtxfit->Swim(0);
            vtxfit->BuildVirtualParticle(0);
            WTrackParameter wlambda = vtxfit->wVirtualTrack(0);
            VertexParameter vtxlambda = vtxfit->vpar(0);

            SecondVertexFit *vtxfit1 = SecondVertexFit::instance();
            vtxfit1->init();
            vtxfit1->setPrimaryVertex(vx_db);
            vtxfit1->AddTrack(0, wlambda);
            vtxfit1->setVpar(vtxlambda);
            if (vtxfit1->Fit())
            {
                HepLorentzVector p4lambda = vtxfit1->p4par();
                HepLorentzVector p4lambdabar = ecms - p4lambda;

                wlambda_vertex = vtxfit1->wpar();
                m_lambda_decayL = vtxfit1->decayLength();
                m_lambda_decayLerr = vtxfit1->decayLengthError();
                m_lambda_massv2 = p4lambda.m();
                m_lambda_pv2 = p4lambda.rho();
                m_lambdabar_massv2 = p4lambdabar.m();
                m_lambdabar_pv2 = p4lambdabar.rho();

                good_lambda = true;
                lambda_chi2 = vtxfit1->chisq();
                lambda_massv2 = p4lambda.m();
            }
        }

        m_lambda_chi2 = lambda_chi2;

        //============================= Tag lambdabar ==========================================
        Ncut_kin++;
    }
    m_tuple0->write();

    return StatusCode::SUCCESS;
}

//---------------------------------------------------------------------------------
StatusCode Taglampimpp::finalize()
{
    cout << "Selection criteria:	"
         << "Survived : sigma pnpipi  V 0 - 4" << endl;
    cout << "Total number:			" << Ncut_none << endl;
    cout << "cut mc:				" << Ncut_mc << endl;
    cout << "Charged trk:			" << Ncut_ncharge << endl;
    cout << "Pid sel:				" << Ncut_pid << endl;
    cout << "Mass sel:				" << Ncut_mass << endl;
    cout << "Length sel:				" << Ncut_length << endl;
    cout << "Chisq sel:				" << Ncut_chisq << endl;
    cout << "nshower: 				" << Ncut_nshower << endl;
    cout << "npi0	: 				" << Ncut_npi0 << endl;
    cout << "nphoton:				" << Ncut_photon << " " << Liu << endl;
    cout << "nbar :					" << Ncut_nbar << endl;
    cout << "kin :					" << Ncut_kin << endl;
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "in finalize()" << endmsg;
    return StatusCode::SUCCESS;
}
