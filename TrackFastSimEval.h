// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef _TRACK_FAST_SIM_EVAL_
#define _TRACK_FAST_SIM_EVAL_

#include <string>

#include <fun4all/SubsysReco.h>

class PHCompositeNode;
class TH1D;

class TrackFastSimEval : public SubsysReco {
 public:
  // Default constructor;
  TrackFastSimEval(const std::string& name = "TrackFastSimEval",
                   const std::string& filename = "g4eval.root",
                   const std::string& trackmapname = "SvtxTrackMap");

  // Initialization;
  int Init(PHCompositeNode*);

  // Process Event;
  int process_event(PHCompositeNode*);

  // End, write and close files;
  int End(PHCompositeNode*);

  // Change output filename;
  void set_filename(const char* file) { if (file) _outfile_name = file; };

 private:
  // Output filename;
  std::string _outfile_name;

  // Name of SvtxTrackMap collection;
  std::string _trackmapname;

  // Event counter;
  int _event;

  // 1D dp/p histogram;
  TH1D *_h1d_Delta_mom;
  
  //1D Eta Distribution Hist;
  TH1D *_eta_dist;
  
};

#endif
