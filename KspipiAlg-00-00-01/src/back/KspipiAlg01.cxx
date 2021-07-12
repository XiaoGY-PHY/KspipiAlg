#include "GaudiKernel/MsgStream.h" //we write this package to obtain all single channel
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/PropertyMgr.h"
#include "VertexFit/IVertexDbSvc.h"
#include "GaudiKernel/Bootstrap.h"
#include "GaudiKernel/ISvcLocator.h"
#include "DTagTool/DTagTool.h"
#include "DsDsmcmode/DsDsmcmode.h"
#include "EventModel/EventModel.h"
#include "EventModel/Event.h"
#include "EvtRecEvent/EvtRecEvent.h"
#include "EvtRecEvent/EvtRecTrack.h"
#include "DstEvent/TofHitStatus.h"
#include "EventModel/EventHeader.h"
#include "Ecmset/Ecmset.h"
#include "VertexFit/IVertexDbSvc.h"
#include "VertexFit/SecondVertexFit.h"
#include "VertexFit/VertexFit.h"
#include "TMath.h"
#include "GaudiKernel/INTupleSvc.h"
#include "GaudiKernel/NTuple.h"
#include "GaudiKernel/Bootstrap.h"
#include "GaudiKernel/IHistogramSvc.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/TwoVector.h"
using CLHEP::Hep3Vector;
using CLHEP::HepLorentzVector;
#include "CLHEP/Geometry/Point3D.h"
#ifndef ENABLE_BACKWARDS_COMPATIBILITY
typedef HepGeom::Point3D<double> HepPoint3D;
#endif
#include "KspipiAlg/KspipiAlg.h"
#include "VertexFit/KinematicFit.h"
#include "VertexFit/KalmanKinematicFit.h"
#include "VertexFit/VertexFit.h"
#include "VertexFit/Helix.h"
#include "ParticleID/ParticleID.h"
#include <vector>

const double mpi = 0.13957;

KspipiAlg::KspipiAlg(const std::string &name, ISvcLocator *pSvcLocator) : Algorithm(name, pSvcLocator)
{
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
StatusCode KspipiAlg::initialize()
{
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "in initialize()" << endmsg;
    StatusCode status;

    NTuplePtr nt0(ntupleSvc(), "FILE1/Kspipirec");
    if (nt0)
        m_tuple0 = nt0;
    else
    {
        m_tuple0 = ntupleSvc()->book("FILE1/Kspipirec", CLID_ColumnWiseTuple, "ks N-Tuple example");
        if (m_tuple0)
        {
            status = m_tuple0->addItem("event", event_temp0);
            status = m_tuple0->addItem("run", run_temp0);
            status = m_tuple0->addItem("KsdecayL", KsdecayL_temp0);
            status = m_tuple0->addItem("ksdecayLerr", KsdecayLerr_temp0);
            status = m_tuple0->addItem("mKs", mKs_temp0);

            //vectors
            status = m_tuple0->addItem("index", index_temp0, 0, 10);
            status = m_tuple0->addIndexedItem("Vpi_P", index_temp0, V_P_pi_temp0);
            status = m_tuple0->addIndexedItem("Vpi_M", index_temp0, V_M_pi_temp0);
            status = m_tuple0->addIndexedItem("VKs", index_temp0, V_Ks_temp0);

            //decay chain
            status = m_tuple0->addItem("indexmc", idxmc_temp0, 0, 100);
            status = m_tuple0->addIndexedItem("trkidx", idxmc_temp0, trkidx_temp0);
            status = m_tuple0->addIndexedItem("pdgid", idxmc_temp0, pdgid_temp0);
            status = m_tuple0->addIndexedItem("motheridx", idxmc_temp0, motheridx_temp0);
        }
        else
        {
            log << MSG::ERROR << "Cannot book N-tuple:" << long(m_tuple0) << endmsg;
            return StatusCode::FAILURE;
        }
    }
}

StatusCode KspipiAlg::execute()
{
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "in execute()" << endreq;
    SmartDataPtr<Event::EventHeader> eventHeader(eventSvc(), "/Event/EventHeader");

    int runNo = eventHeader->runNumber();
    int event = eventHeader->eventNumber();

    log << MSG::DEBUG << "run, evtnum = "
        << runNo << " , "
        << event << endreq;

    //decay chain
    int idx_mc = -9999, pdgid_mc[100], motheridx_mc[100], trkidx_mc[100];
    for (int i = 0; i < 100; i++)
    {
        pdgid_mc[i] = -9999;
        motheridx_mc[i] = -9999;
        trkidx_mc[100] = -9999;
    }
    if (runNo < 0)
    {
        SmartDataPtr<Event::McParticleCol> mcParticleCol(eventSvc(), "/Event/MC/McParticleCol");
        int n_numParticle = 1;
        if (!mcParticleCol)
        {
            return StatusCode::FAILURE;
        }
        else
        {
            Event::McParticleCol::iterator iter_mc = mcParticleCol->begin();
            pdgid_mc[0] = 80022;
            motheridx_mc[0] = -1;
            trkidx_mc[0] = 0;
            for (; iter_mc != mcParticleCol->end(); iter_mc++)
            {
                if (!(*iter_mc)->decayFromGenerator())
                    continue;
                int mcidx = ((*iter_mc)->mother()).trackIndex();
                int pdgid = (*iter_mc)->particleProperty();
                int trkidx = (*iter_mc)->trackIndex();
                int motherid = ((*iter_mc)->mother()).particleProperty();
                if (mcidx == trkidx)
                {
                    motheridx_mc[n_numParticle] = 0;
                }
                else
                {
                    motheridx_mc[n_numParticle] = mcidx + 1;
                }
                trkidx_mc[n_numParticle] = trkidx + 1;
                pdgid_mc[n_numParticle] = pdgid;
                n_numParticle++;
            }
            idx_mc = n_numParticle;
        }
    }

    SmartDataPtr<EvtRecEvent> evtRecEvent(eventSvc(), EventModel::EvtRec::EvtRecEvent);
    log << MSG::DEBUG << "ncharg, nneu, tottks = "
        << evtRecEvent->totalCharged() << " , "
        << evtRecEvent->totalNeutral() << " , "
        << evtRecEvent->totalTracks() << endreq;

    SmartDataPtr<EvtRecTrackCol> evtRecTrkCol(eventSvc(), EventModel::EvtRec::EvtRecTrackCol);
    ParticleID *pid = ParticleID::instance();

    Hep3Vector xorigin(0, 0, 0);
    HepSymMatrix xoriginEx(3, 0);
    IVertexDbSvc *vtxsvc;
    Gaudi::svcLocator()->service("VertexDbSvc", vtxsvc);
    if (vtxsvc->isVertexValid())
    {
        double *dbv = vtxsvc->PrimaryVertex();
        double *vv = vtxsvc->SigmaPrimaryVertex();
        xorigin.setX(dbv[0]);
        xorigin.setY(dbv[1]);
        xorigin.setZ(dbv[2]);
        xoriginEx[0][0] = vv[0] * vv[0];
        xoriginEx[1][1] = vv[1] * vv[1];
        xoriginEx[2][2] = vv[2] * vv[2];
    }

    HepPoint3D point0(0., 0., 0.);
    HepPoint3D IP(xorigin[0], xorigin[1], xorigin[2]);

    double mKs, KsdecayL, KsdecayLerr;
    HepLorentzVector V_P_pi;
    HepLorentzVector V_M_pi;
    HepLorentzVector V_Ks;
    ////////////////////////////////////////select Kso
    for (int i = 0; i < evtRecEvent->totalCharged(); i++)
    {
        // pi+
        EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + i;
        RecMdcKalTrack::setPidType(RecMdcKalTrack::pion);
        if (!(*itTrk)->isMdcKalTrackValid())
            continue;
        RecMdcKalTrack *mdcKalTrk1 = (*itTrk)->mdcKalTrack();
        if (mdcKalTrk1->charge() != 1)
            continue;
        HepVector a1 = mdcKalTrk1->getZHelix();
        HepSymMatrix Ea1 = mdcKalTrk1->getZError();
        VFHelix helixip3_1(point0, a1, Ea1);
        helixip3_1.pivot(IP);
        // HepVector vecipa1 = helixip3_1.a();
        // double dr1 = fabs(vecipa1[0]);
        // double drz = fabs(vecipa1[3]);
        double costheta = cos(mdcKalTrk1->theta());
        // if (drz >= 20.0)
        //     continue;
        if (fabs(costheta) >= 0.93)
            continue;

        // pi-
        for (int j = 0; j < evtRecEvent->totalCharged(); j++)
        {
            EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + j;
            if (!(*itTrk)->isMdcKalTrackValid())
                continue;
            RecMdcKalTrack::setPidType(RecMdcKalTrack::pion);
            RecMdcKalTrack *mdcKalTrk2 = (*itTrk)->mdcKalTrack();
            if (mdcKalTrk2->charge() != -1)
                continue;
            HepVector a2 = mdcKalTrk2->getZHelix();
            HepSymMatrix Ea2 = mdcKalTrk2->getZError();
            VFHelix helixip3_2(point0, a2, Ea2);
            helixip3_2.pivot(IP);
            HepVector vecipa2 = helixip3_2.a();
            // double dr1 = fabs(vecipa2[0]);
            // double drz = fabs(vecipa2[3]);
            double costheta = cos(mdcKalTrk2->theta());
            // if (drz >= 20.0)
            //     continue;
            if (fabs(costheta) >= 0.93)
                continue;

            WTrackParameter pip(mpi, mdcKalTrk1->getZHelix(), mdcKalTrk1->getZError());
            WTrackParameter pim(mpi, mdcKalTrk2->getZHelix(), mdcKalTrk2->getZError());
            HepVector pip_val = HepVector(7, 0);
            HepVector pim_val = HepVector(7, 0);
            pip_val = pip.w();
            pim_val = pim.w();
            HepLorentzVector ptrktagk0(pip_val[0] + pim_val[0], pip_val[1] + pim_val[1], pip_val[2] + pim_val[2], pip_val[3] + pim_val[3]);
            double m_xmtagk0_tem = ptrktagk0.m();

            if (fabs(ptrktagk0.m() - 0.498) > 0.1)
                continue;

            Hep3Vector ip(0, 0, 0);
            HepSymMatrix ipEx(3, 0);
            IVertexDbSvc *vtxsvc;
            Gaudi::svcLocator()->service("VertexDbSvc", vtxsvc);
            if (vtxsvc->isVertexValid())
            {
                double *dbv = vtxsvc->PrimaryVertex();
                double *vv = vtxsvc->SigmaPrimaryVertex();
                ip.setX(dbv[0]);
                ip.setY(dbv[1]);
                ip.setZ(dbv[2]);
                ipEx[0][0] = vv[0] * vv[0];
                ipEx[1][1] = vv[1] * vv[1];
                ipEx[2][2] = vv[2] * vv[2];
            }
            else
            {
                continue;
            }

            VertexParameter bs;
            bs.setVx(ip);
            bs.setEvx(ipEx);
            ///////////////////////////////////////////////////////////////////
            //  set a common vertex with huge error
            ///////////////////////////////////////////////////////////////////
            HepPoint3D vx(0., 0., 0.);
            HepSymMatrix Evx(3, 0);
            double bx = 1E+6;
            double by = 1E+6;
            double bz = 1E+6;
            Evx[0][0] = bx * bx;
            Evx[1][1] = by * by;
            Evx[2][2] = bz * bz;
            VertexParameter vxpar;
            vxpar.setVx(vx);
            vxpar.setEvx(Evx);
            ///////////////////////////////////////////////////////////////////
            //  do vertex fit
            ///////////////////////////////////////////////////////////////////
            VertexFit *vtxfit = VertexFit::instance();
            vtxfit->init();

            vtxfit->AddTrack(0, pip);
            vtxfit->AddTrack(1, pim);
            vtxfit->AddVertex(0, vxpar, 0, 1);
            if (!(vtxfit->Fit(0)))
                continue;
            vtxfit->Swim(0);
            vtxfit->BuildVirtualParticle(0);
            WTrackParameter wks_Trk = vtxfit->wVirtualTrack(0);
            double vtx_chisq = vtxfit->chisq(0);

            WTrackParameter wksp = vtxfit->wtrk(0);
            WTrackParameter wksm = vtxfit->wtrk(1);

            HepVector ks_pip = HepVector(7, 0);
            ks_pip = wksp.w();
            HepLorentzVector P_pip(ks_pip[0], ks_pip[1], ks_pip[2], ks_pip[3]);
            P_pip.boost(-0.011, 0, 0);
            V_P_pi = P_pip;

            HepVector ks_pim = HepVector(7, 0);
            ks_pim = wksm.w();
            HepLorentzVector P_pim(ks_pim[0], ks_pim[1], ks_pim[2], ks_pim[3]);
            P_pim.boost(-0.011, 0, 0);
            V_M_pi = P_pim;
            ///////////////////////////////////////////////////////////////////
            //  do second vertex fit
            ///////////////////////////////////////////////////////////////////
            SecondVertexFit *svtxfit = SecondVertexFit::instance();
            svtxfit->init();
            svtxfit->setPrimaryVertex(bs);
            svtxfit->AddTrack(0, vtxfit->wVirtualTrack(0));
            svtxfit->setVpar(vtxfit->vpar(0));
            if (!svtxfit->Fit())
                continue;
            if (svtxfit->chisq() > 200.)
                continue;
            if (svtxfit->decayLength() < 0.0)
                continue;
            double len = svtxfit->decayLength();
            double lenerr = svtxfit->decayLengthError();

            KsdecayL = len;
            KsdecayLerr = lenerr;

            if ((len / lenerr) <= 2)
                continue;
            HepLorentzVector pks_temp = svtxfit->p4par();
            //       if( fabs(pks_temp.m()-0.4977)>0.012) continue;
            if (pks_temp.m() < 0.487 || pks_temp.m() > 0.511)
                continue;
            mKs = pks_temp.m();
            pks_temp.boost(-0.011, 0, 0);
            V_Ks = pks_temp;
        }
    }

    run_temp0 = runNo;
    event_temp0 = event;
    KsdecayL_temp0 = KsdecayL;
    KsdecayLerr_temp0 = KsdecayLerr;
    mKs_temp0 = mKs;

    index_temp0 = 4;

    V_P_pi_temp0[0] = V_P_pi.px();
    V_P_pi_temp0[1] = V_P_pi.py();
    V_P_pi_temp0[2] = V_P_pi.pz();
    V_P_pi_temp0[3] = V_P_pi.e();

    V_M_pi_temp0[0] = V_M_pi.px();
    V_M_pi_temp0[1] = V_M_pi.py();
    V_M_pi_temp0[2] = V_M_pi.pz();
    V_M_pi_temp0[3] = V_M_pi.e();

    V_Ks_temp0[0] = V_Ks.px();
    V_Ks_temp0[1] = V_Ks.py();
    V_Ks_temp0[2] = V_Ks.pz();
    V_Ks_temp0[3] = V_Ks.e();

    // decaychain
    for (int ntk_mc = 0; ntk_mc < idx_mc; ntk_mc++)
    {
        pdgid_temp0[ntk_mc] = pdgid_mc[ntk_mc];
        motheridx_temp0[ntk_mc] = motheridx_mc[ntk_mc];
        trkidx_temp0[ntk_mc] = trkidx_mc[ntk_mc];
    }
    idxmc_temp0 = idx_mc;

    m_tuple0->write();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
StatusCode KspipiAlg::finalize()
{
    MsgStream log(msgSvc(), name()); /////////////////////
    log << MSG::INFO << "in finalize()" << endmsg;
    return StatusCode::SUCCESS;
}
