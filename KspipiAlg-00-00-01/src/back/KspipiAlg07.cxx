#include "GaudiKernel/MsgStream.h" //we write this package to obtain all single channel
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

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Geometry/Point3D.h"
using CLHEP::Hep3Vector;
using CLHEP::HepLorentzVector;

#include "KspipiAlg/KspipiAlg.h"

#ifndef ENABLE_BACKWARDS_COMPATIBILITY
typedef HepGeom::Point3D<double> HepPoint3D;
#endif

#include <vector>

const double mpi = 0.13957;

KspipiAlg::KspipiAlg(const std::string &name, ISvcLocator *pSvcLocator) : Algorithm(name, pSvcLocator)
{
    declareProperty("cms", ecms = 3.097);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
StatusCode KspipiAlg::initialize()
{
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "in initialize()" << endmsg;
    StatusCode status;

    NTuplePtr nt0(ntupleSvc(), "FILE1/Ksrec");
    if (nt0)
        m_tuple0 = nt0;
    else
    {
        m_tuple0 = ntupleSvc()->book("FILE1/Ksrec", CLID_ColumnWiseTuple, "ks N-Tuple example");
        if (m_tuple0)
        {
            status = m_tuple0->addItem("event", event_temp0);
            status = m_tuple0->addItem("run", run_temp0);
            status = m_tuple0->addItem("KsdecayL", KsdecayL_temp0);
            status = m_tuple0->addItem("ksdecayLerr", KsdecayLerr_temp0);
            status = m_tuple0->addItem("fKsdecayL", fKsdecayL_temp0);
            status = m_tuple0->addItem("fksdecayLerr", fKsdecayLerr_temp0);
            status = m_tuple0->addItem("truth_Ks_x", truthKs_x_temp0);
            status = m_tuple0->addItem("truth_Ks_y", truthKs_y_temp0);
            status = m_tuple0->addItem("truth_Ks_z", truthKs_z_temp0);
            status = m_tuple0->addItem("truth_Vertex_x", truthvertex_x_temp0);
            status = m_tuple0->addItem("truth_Vertex_y", truthvertex_y_temp0);
            status = m_tuple0->addItem("truth_Vertex_z", truthvertex_z_temp0);
            status = m_tuple0->addItem("truth_decayL", truthdecayL_temp0);
            status = m_tuple0->addItem("truth_decayL0", truthdecayL0_temp0);
            status = m_tuple0->addItem("cutReason", cutReason_temp0);
            status = m_tuple0->addItem("fcutReason", fcutReason_temp0);
            status = m_tuple0->addItem("ncharge", ncharge_temp0);
            status = m_tuple0->addItem("mKs", mKs_temp0);
            status = m_tuple0->addItem("fmKs", fmKs_temp0);
            status = m_tuple0->addItem("m2pi", m2pi_temp0);
            status = m_tuple0->addItem("fm2pi", fm2pi_temp0);

            //vectors
            status = m_tuple0->addItem("index", index_temp0, 0, 10);
            status = m_tuple0->addIndexedItem("Vpi_P", index_temp0, V_P_pi_temp0);
            status = m_tuple0->addIndexedItem("Vpi_M", index_temp0, V_M_pi_temp0);
            status = m_tuple0->addIndexedItem("VKs", index_temp0, V_Ks_temp0);
            status = m_tuple0->addIndexedItem("fVpi_P", index_temp0, fV_P_pi_temp0);
            status = m_tuple0->addIndexedItem("fVpi_M", index_temp0, fV_M_pi_temp0);
            status = m_tuple0->addIndexedItem("fVKs", index_temp0, fV_Ks_temp0);

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
//int NNN = 0;
StatusCode KspipiAlg::execute()
{
    //cout<<NNN<<endl;
    //NNN++;
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
            pdgid_mc[0] = 443;
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

    //*********************************************************************************************
    //truth
    HepLorentzVector Ks_fposi;
    HepLorentzVector Ks_iposi;
    if (runNo < 0)
    {
        SmartDataPtr<Event::McParticleCol> mcParticleCol(eventSvc(), "/Event/MC/McParticleCol");
        if (mcParticleCol)
        {
            Event::McParticleCol::iterator iter_mc = mcParticleCol->begin();
            for (; iter_mc != mcParticleCol->end(); iter_mc++)
            {
                // if (!(*iter_mc)->decayFromGenerator())
                //     continue;
                if ((*iter_mc)->particleProperty() == 310)
                {
                    // cout << "0000000000000000000000001" << endl;
                    //cout << ((*iter_mc)->mother()).particleProperty() << endl;
                    Ks_iposi = (*iter_mc)->initialPosition();
                    Ks_fposi = (*iter_mc)->finalPosition();
                }
                // if ((*iter_mc)->particleProperty() == 443)
                // {
                //     cout << "0000000000000000000000002" << endl;
                // }
                // if ((*iter_mc)->particleProperty() == 310 && ((*iter_mc)->mother()).particleProperty() == 443)
                // {
                //     cout << "0000000000000000000000003" << endl;
                //     Ks_fposi = (*iter_mc)->finalPosition();
                // }
            }
        }
    }
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
    // VertexParameter vx_db;
    // vx_db.setVx(vx);
    // vx_db.setEvx(Evx);

    //	MCtruth decay length
    truthKs_x_temp0 = Ks_fposi.x();
    truthKs_y_temp0 = Ks_fposi.y();
    truthKs_z_temp0 = Ks_fposi.z();
    truthvertex_x_temp0 = vx.x();
    truthvertex_y_temp0 = vx.y();
    truthvertex_z_temp0 = vx.z();

    double m_mcdecayL, m_mcdecayL0;
    m_mcdecayL = sqrt((Ks_fposi.x() - vx.x()) * (Ks_fposi.x() - vx.x()) + (Ks_fposi.y() - vx.y()) * (Ks_fposi.y() - vx.y()) + (Ks_fposi.z() - vx.z()) * (Ks_fposi.z() - vx.z()));
    m_mcdecayL0 = sqrt((Ks_fposi.x() - Ks_iposi.x()) * (Ks_fposi.x() - Ks_iposi.x()) + (Ks_fposi.y() - Ks_iposi.y()) * (Ks_fposi.y() - Ks_iposi.y()) + (Ks_fposi.z() - Ks_iposi.z()) * (Ks_fposi.z() - Ks_iposi.z()));

    //**********************************************************************************************
    //rec
    SmartDataPtr<EvtRecEvent> evtRecEvent(eventSvc(), EventModel::EvtRec::EvtRecEvent);
    log << MSG::DEBUG << "ncharg, nneu, tottks = "
        << evtRecEvent->totalCharged() << " , "
        << evtRecEvent->totalNeutral() << " , "
        << evtRecEvent->totalTracks() << endreq;

    SmartDataPtr<EvtRecTrackCol> evtRecTrkCol(eventSvc(), EventModel::EvtRec::EvtRecTrackCol);
    ParticleID *pid = ParticleID::instance();

    HepPoint3D point0(0., 0., 0.);
    HepPoint3D IP(xorigin[0], xorigin[1], xorigin[2]);

    double mKs = -1000;
    double m2pi = -1000;
    double KsdecayL = -1000;
    double KsdecayLerr = -1000;
    int cutReason = -1000;
    double masscut = 9999;
    HepLorentzVector V_P_pi;
    HepLorentzVector V_M_pi;
    HepLorentzVector V_Ks;

    int ncharge = evtRecEvent->totalCharged();

    ////////////////////////////////////////select Kso
    for (int zf = 0; zf < 2; zf++)
    {
        // only one charged track
        if (evtRecEvent->totalCharged() < 2)
        {
            cutReason = 1; // total charged tracks < 2
        }

        // same charge
        if (evtRecEvent->totalCharged() == 2)
        {
            EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + 0;
            RecMdcKalTrack *mdcKalTrk1 = (*itTrk)->mdcKalTrack();
            EvtRecTrackIterator itTrk2 = evtRecTrkCol->begin() + 1;
            RecMdcKalTrack *mdcKalTrk2 = (*itTrk2)->mdcKalTrack();
            RecMdcKalTrack::setPidType(RecMdcKalTrack::pion);
            if ((*itTrk)->isMdcKalTrackValid())
            {
                //cout << mdcKalTrk1->charge() << endl;
                if ((*itTrk2)->isMdcKalTrackValid())
                {
                    //cout << mdcKalTrk2->charge() << endl;
                    if ((mdcKalTrk1->charge()) * (mdcKalTrk2->charge()) == 1)
                    {
                        cutReason = 2; // same charge
                    }
                }
            }
        }
  
        for (int i = 0; i < evtRecEvent->totalCharged(); i++)
        {
            if (cutReason == 1 || cutReason == 2)
                continue;
            // pi+
            EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + i;
            RecMdcKalTrack::setPidType(RecMdcKalTrack::pion);
            if (!(*itTrk)->isMdcKalTrackValid())
                continue;
            RecMdcKalTrack *mdcKalTrk1 = (*itTrk)->mdcKalTrack();
            if (mdcKalTrk1->charge() != 1)
                continue;
            HepVector a1;
            HepSymMatrix Ea1;
            if (zf == 0)
            {
                a1 = mdcKalTrk1->getZHelix();
                Ea1 = mdcKalTrk1->getZError();
            }
            else if (zf == 1)
            {
                a1 = mdcKalTrk1->getFHelix();
                Ea1 = mdcKalTrk1->getFError();
            }
            VFHelix helixip3_1(point0, a1, Ea1);
            helixip3_1.pivot(IP);
            // HepVector vecipa1 = helixip3_1.a();
            // double dr1 = fabs(vecipa1[0]);
            // double drz = fabs(vecipa1[3]);
            double costheta = cos(mdcKalTrk1->theta());
            // if (drz >= 20.0)
            //     continue;
            if (fabs(costheta) >= 0.93)
            {
                if (cutReason == 0)
                {
                    continue;
                }
                cutReason = 8;
                continue;
            }

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
                HepVector a2;
                HepSymMatrix Ea2;
                if (zf == 0)
                {
                    a2 = mdcKalTrk2->getZHelix();
                    Ea2 = mdcKalTrk2->getZError();
                }
                else if (zf == 1)
                {
                    a2 = mdcKalTrk2->getFHelix();
                    Ea2 = mdcKalTrk2->getFError();
                }
                VFHelix helixip3_2(point0, a2, Ea2);
                helixip3_2.pivot(IP);
                HepVector vecipa2 = helixip3_2.a();
                // double dr1 = fabs(vecipa2[0]);
                // double drz = fabs(vecipa2[3]);
                double costheta = cos(mdcKalTrk2->theta());
                // if (drz >= 20.0)
                //     continue;
                if (fabs(costheta) >= 0.93)
                {
                    if (cutReason == 0)
                    {
                        continue;
                    }
                    cutReason = -8;
                    continue;
                }

                WTrackParameter pip(mpi, a1, Ea1);
                WTrackParameter pim(mpi, a2, Ea2);

                HepVector pip_val = HepVector(7, 0);
                HepVector pim_val = HepVector(7, 0);
                pip_val = pip.w(); //啥意思？
                pim_val = pim.w();
                HepLorentzVector ptrktagk0(pip_val[0] + pim_val[0], pip_val[1] + pim_val[1], pip_val[2] + pim_val[2], pip_val[3] + pim_val[3]);
                double m_xmtagk0_tem = ptrktagk0.m();

                //cout << "bbbbbbbbbbbbbbbb! :: " << NNN << endl;
                if (fabs(ptrktagk0.m() - 0.497611) > 0.2)
                {
                    if (cutReason == 0)
                    {
                        continue;
                    }
                    cutReason = 3; // abs(mKs - 0.498) > 0.2
                    continue;
                }
                Hep3Vector ip(0, 0, 0);
                HepSymMatrix ipEx(3, 0);
                IVertexDbSvc *vtxsvc;
                Gaudi::svcLocator()->service("VertexDbSvc", vtxsvc);
                if (vtxsvc->isVertexValid())
                {
                    double *dbv = vtxsvc->PrimaryVertex();
                    double *vv = vtxsvc->SigmaPrimaryVertex(); //err
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
                double bx = 1E+6; //
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
                {
                    if (cutReason == 0)
                    {
                        continue;
                    }
                    cutReason = 4; // first fit failed
                    continue;
                }
                //cout << "cccccccccccccccc! :: " << NNN << endl;
                vtxfit->Swim(0);
                vtxfit->BuildVirtualParticle(0);
                WTrackParameter wks_Trk = vtxfit->wVirtualTrack(0);
                double vtx_chisq = vtxfit->chisq(0);

                WTrackParameter wksp = vtxfit->wtrk(0);
                WTrackParameter wksm = vtxfit->wtrk(1);

                HepVector ks_pip = HepVector(7, 0);
                ks_pip = wksp.w();
                HepLorentzVector P_pip(ks_pip[0], ks_pip[1], ks_pip[2], ks_pip[3]);
                P_pip.boost(-0.011, 0, 0); //到质心系

                HepVector ks_pim = HepVector(7, 0);
                ks_pim = wksm.w();
                HepLorentzVector P_pim(ks_pim[0], ks_pim[1], ks_pim[2], ks_pim[3]);
                P_pim.boost(-0.011, 0, 0);

                ///////////////////////////////////////////////////////////////////
                //  do second vertex fit
                ///////////////////////////////////////////////////////////////////
                SecondVertexFit *svtxfit = SecondVertexFit::instance();
                svtxfit->init();
                svtxfit->setPrimaryVertex(bs);
                svtxfit->AddTrack(0, vtxfit->wVirtualTrack(0));
                svtxfit->setVpar(vtxfit->vpar(0));
                if (!svtxfit->Fit())
                {
                    if (cutReason == 0)
                    {
                        continue;
                    }
                    cutReason = 5; // second fit failed
                    continue;
                }
                if (svtxfit->chisq() > 200.)
                {
                    if (cutReason == 0)
                    {
                        continue;
                    }
                    cutReason = 6; // chi2 > 200
                    continue;
                }
                // if (svtxfit->decayLength() < 0.0) //397
                // {
                //     if (cutReason == 0)
                //     {
                //         continue;
                //     }
                //     cutReason = 7; // fit decayL < 0
                //     continue;
                // }
                double len = svtxfit->decayLength();
                double lenerr = svtxfit->decayLengthError();

                if ((len / lenerr) <= 2)
                {
                    if (cutReason == 0)
                    {
                        continue;
                    }
                    cutReason = 7; // (len / lenerr) <= 2
                    continue;
                }
                HepLorentzVector pks_temp = svtxfit->p4par();

                //       if( fabs(pks_temp.m()-0.4977)>0.012) continue;
                // if (pks_temp.m() < 0.487 || pks_temp.m() > 0.511)
                // {
                //     if (cutReason == 0)
                //     {
                //         continue;
                //     }
                //     cutReason = 9; // cut by tight mKs condition
                //     continue;
                // }
                //cout << "dddddddddddddddd! :: " << NNN << endl;

                if (abs(pks_temp.m() - 0.497611) < masscut)
                {
                    masscut = abs(pks_temp.m() - 0.497611);
                    m2pi = (P_pip + P_pim).m();
                    mKs = pks_temp.m();
                    pks_temp.boost(-0.011, 0, 0);
                    V_Ks = pks_temp;
                    V_P_pi = P_pip;
                    V_M_pi = P_pim;
                    KsdecayL = len;
                    KsdecayLerr = lenerr;
                    cutReason = 0;
                }
            }
        }

        if (zf == 0)
        {
            KsdecayL_temp0 = KsdecayL;
            KsdecayLerr_temp0 = KsdecayLerr;
            cutReason_temp0 = cutReason;
            mKs_temp0 = mKs;
            m2pi_temp0 = m2pi;

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
        }
        else if (zf == 1)
        {
            fKsdecayL_temp0 = KsdecayL;
            fKsdecayLerr_temp0 = KsdecayLerr;
            fcutReason_temp0 = cutReason;
            fmKs_temp0 = mKs;
            fm2pi_temp0 = m2pi;

            fV_P_pi_temp0[0] = V_P_pi.px();
            fV_P_pi_temp0[1] = V_P_pi.py();
            fV_P_pi_temp0[2] = V_P_pi.pz();
            fV_P_pi_temp0[3] = V_P_pi.e();

            fV_M_pi_temp0[0] = V_M_pi.px();
            fV_M_pi_temp0[1] = V_M_pi.py();
            fV_M_pi_temp0[2] = V_M_pi.pz();
            fV_M_pi_temp0[3] = V_M_pi.e();

            fV_Ks_temp0[0] = V_Ks.px();
            fV_Ks_temp0[1] = V_Ks.py();
            fV_Ks_temp0[2] = V_Ks.pz();
            fV_Ks_temp0[3] = V_Ks.e();
        }
        
        cutReason = -1000;
    }

    run_temp0 = runNo;
    event_temp0 = event;
    truthdecayL_temp0 = m_mcdecayL;
    truthdecayL0_temp0 = m_mcdecayL0;
    ncharge_temp0 = ncharge;
    index_temp0 = 4;

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
