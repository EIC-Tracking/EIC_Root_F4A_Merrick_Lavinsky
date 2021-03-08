
#include <iostream>

#include <TVector3.h>
#include <TH1D.h>
#include <TGraph.h>

#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrack_FastSim.h>

#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/PHTFileServer.h>
#include <fun4all/SubsysReco.h>  

#include <phool/getClass.h>

#include "TrackFastSimEval.h"

using namespace std;

// ---------------------------------------------------------------------------------------

TrackFastSimEval::TrackFastSimEval(const string &name, const string &filename, 
				   const string &trackmapname)
  : SubsysReco(name)
  , _outfile_name(filename)
  , _trackmapname(trackmapname)
  , _event(0)
  , _h1d_Delta_mom(nullptr)
  , _eta_dist(nullptr)
  , ED_Gen(nullptr)
  , Eta_res(nullptr)
{
} // TrackFastSimEval::TrackFastSimEval()

// ---------------------------------------------------------------------------------------

int TrackFastSimEval::Init(PHCompositeNode *topNode)
{
  cout << PHWHERE << " Opening file " << _outfile_name << endl;
  PHTFileServer::get().open(_outfile_name, "RECREATE");

  _h1d_Delta_mom = new TH1D("_h1d_Delta_mom", "Momentum Resolution in Hybrid Model; Momentum Resolution; Counts", 100, -0.2, 0.2);

  _eta_dist = new TH1D("_eta_dist", "Event #eta Distribution of Tracks; #eta; Counts", 100, -5, 5);

  ED_Gen = new TH1D("ED_Gen", "Detector #eta Distribution of Generated Hits; #eta; Counts", 100, -5, 5);

  Eta_res = new TH1D("Eta_res", "#eta Resolution; #eta Resolution; Counts", 100, -5, 5);

  Eta_rez = new TGraph();

  //cout << "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAaaaaaaaaaaa" << endl;
  return Fun4AllReturnCodes::EVENT_OK;
} // TrackFastSimEval::Init()

// ---------------------------------------------------------------------------------------

int TrackFastSimEval::process_event(PHCompositeNode *topNode)
{
  //n=0;
  //cout << "BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB" << endl;
  if (++_event % 100 == 0) cout << PHWHERE << "Events processed: " << _event << endl;

  auto _truth_container = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  auto _trackmap = findNode::getClass<SvtxTrackMap>(topNode, _trackmapname);
  if (!_truth_container || !_trackmap) {
    cout << PHWHERE << " Either PHG4TruthInfoContainer or SvtxTrackMap node not found on node tree"
         << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  } //if

  {
    auto range = _truth_container->GetPrimaryParticleRange();
    //cout << "1000?" << range << endl;
    //cout << g4particle->phg4hitsNames << endl;
    
    // Loop through all the truth particles;
    for (auto truth_itr = range.first; truth_itr != range.second; ++truth_itr) {
      auto g4particle = truth_itr->second;
      if (!g4particle) continue;

      TVector3 truth_mom(g4particle->get_px(), g4particle->get_py(), g4particle->get_pz());
      //TVector3 reco_mom (track->get_px(), track->get_py(), track->get_pz());
      float theta = acos(g4particle->get_pz()/truth_mom.Mag());
      float eta;
      if(theta>0){
	eta = -log(tan(theta/2));}
      else{ eta = -(-log(tan((M_PI + theta)/2)));
      }
      //std::cout<<" theta "<<theta<<" eta "<<eta<<std::endl;
      //cout << g4particle->phg4hitsNames << endl;
      ED_Gen->Fill(eta);
      
      //cout << g4particle->get_name() << endl;
      // Loop through all the reconstructed particles and find a match; how about std::map?;
      for (auto track_itr = _trackmap->begin(); track_itr != _trackmap->end(); track_itr++) {
	auto track = dynamic_cast<SvtxTrack_FastSim *>(track_itr->second);
	if (!track) {
	  std::cout << "ERROR CASTING PARTICLE!" << std::endl;
	  continue;
	} //if


	//cout << "Tracks: " << g4particle->get_truth_track_id() << endl;
	// Matching reconstructed particle found, use it and break;
	if ((track->get_truth_track_id() - g4particle->get_track_id()) == 0) {
	  TVector3 truth_mom(g4particle->get_px(), g4particle->get_py(), g4particle->get_pz());
	  TVector3 reco_mom (     track->get_px(),      track->get_py(),      track->get_pz());

	  _h1d_Delta_mom->Fill((reco_mom.Mag() - truth_mom.Mag()) / truth_mom.Mag());
	  
	  //cout << track->get_eta() << endl;
	  //if()
	  _eta_dist->Fill(track->get_eta());

	  float theta_2 = acos(g4particle->get_pz()/truth_mom.Mag());
	  float eta_2;
	  if(theta_2>0){
	    eta_2 = -log(tan(theta_2/2));}
	  else{ eta_2 = -(-log(tan((M_PI + theta_2)/2)));
	  }
	  //cout << "diff in eta: " << eta_2 - track->get_eta() << ", res: " << (eta_2 - track->get_eta())/eta_2 << endl;
	  if(eta_2 != 0){
	    Eta_res->Fill((eta_2 - track->get_eta())/eta_2);
	    // cout << n << endl;
	    Eta_rez->SetPoint(n, eta_2, (eta_2 - track->get_eta())/eta_2);
	    n++;
	    //}
	    //}	  
	    break;
	  }
	} //if
      } //for 
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
} // TrackFastSimEval::process_event()

// ---------------------------------------------------------------------------------------

int TrackFastSimEval::End(PHCompositeNode *topNode)
{
  //Eta_rez->Draw("ALP");
  PHTFileServer::get().cd(_outfile_name);
  //cout << "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC" << endl;
  _h1d_Delta_mom->Write();
  _eta_dist->Write();
  ED_Gen->Write();
  //Eta_res->Write();
  Eta_rez->Write();
  
  
  return Fun4AllReturnCodes::EVENT_OK;
} // TrackFastSimEval::End()

// ---------------------------------------------------------------------------------------
