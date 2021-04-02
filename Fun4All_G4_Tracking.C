#include <phgenfit/Track.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllDummyInputManager.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllNoSyncDstInputManager.h>
#include <fun4all/Fun4AllOutputManager.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/SubsysReco.h>

#include <g4detectors/PHG4DetectorSubsystem.h>
#include <g4detectors/PHG4CylinderSubsystem.h>

#include <g4histos/G4HitNtuple.h>
#include <GdmlImportDetectorSubsystem.h>      
#include <g4main/PHG4ParticleGenerator.h>
#include <g4main/PHG4ParticleGeneratorBase.h>
#include <g4main/PHG4ParticleGun.h>
#include <g4main/PHG4Reco.h>
#include <g4main/PHG4SimpleEventGenerator.h>
#include <g4main/PHG4TruthSubsystem.h>
#include <g4trackfastsim/PHG4TrackFastSim.h>
#include <g4trackfastsim/PHG4TrackFastSimEval.h>
#include <phool/recoConsts.h>
#include <g4detectors/PHG4GDMLSubsystem.h>

#include <g4tpc/PHG4TpcEndCapSubsystem.h>
#include <g4main/HepMCNodeReader.h>
#include <phhepmc/Fun4AllHepMCInputManager.h>

#include <GdmlImportDetectorSubsystem.h>      
#include <SimpleNtuple.h>                        
#include <TrackFastSimEval.h>  
//
//  EicRoot VST & MM barrel geometry creation and tracking;
//

#include <fun4all/Fun4AllServer.h>
#include <g4main/PHG4Reco.h>
#include <g4main/PHG4Detector.h>
#include <g4main/PHG4ParticleGenerator.h>
#include <g4main/PHG4TruthSubsystem.h>
#include <g4trackfastsim/PHG4TrackFastSim.h>
#include <phool/recoConsts.h>

#include <EicToyModelSubsystem.h>
#include <EicRootVstSubsystem.h>
#include <EicRootGemSubsystem.h>
#include <EicRootMuMegasSubsystem.h>
#include <EtmOrphans.h>

#include <TrackFastSimEval.h>
#include "detector_setup.h"

R__LOAD_LIBRARY(libeictoydetectors.so)
// FIXME: add to CMakeLists.txt;
R__LOAD_LIBRARY(libg4trackfastsim.so)
R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libg4detectors.so)

R__LOAD_LIBRARY(libgdmlimportdetector.so)	
R__LOAD_LIBRARY(libg4tpc.so)


// This scrip is simple, sorry: either Qt display or tracking; uncomment if want to see the geometry; 

//#define _QT_DISPLAY_


//all of these variables are controlled in detecter_setup.h
//#define _TPC_
//#define _TPC_ENDCAPS_
//#define _FORWARD_SILICON_DISKS_
//#define _BACKWARD_SILICON_DISKS_
// #define _BARREL_
// #define _USE_FORWARD_PIPES_



namespace G4TPC
{
  int n_tpc_layer_inner = 16;
  double tpc_layer_thick_inner = 1.25;  // EIC- recover default inner radius of TPC vol.                                                              
  int tpc_layer_rphi_count_inner = 1152;

  int n_tpc_layer_mid = 16;
  double tpc_layer_thick_mid = 1.25;
  int tpc_layer_rphi_count_mid = 1536;

  int n_tpc_layer_outer = 16;
  double tpc_layer_thick_outer = 1.125;  // outer later reach from 60-78 cm (instead of 80 cm), that leads to radial thickness of 1.125 cm            
  int tpc_layer_rphi_count_outer = 2304;

  double outer_radius = 78.;
  double inner_cage_radius = 20.;
  double cage_length = 211.0;  // From TPC group, gives eta = 1.1 at 78 cm                                                                            
  double n_rad_length_cage = 1.13e-02;
  double cage_thickness = 28.6 * n_rad_length_cage;  // Kapton X_0 = 28.6 cm  // mocks up Kapton + carbon fiber structure                             

  string tpcgas = "sPHENIX_TPC_Gas";  //  Ne(90%) CF4(10%) - defined in g4main/PHG4Reco.cc                                                            

}  // namespace G4TPC                                                                                                                                 

void Fun4All_G4_Tracking(int nEvents = 1000)
{
  // The fun4all Server;
  Fun4AllServer *se = Fun4AllServer::instance();

  // May want to fix the random seed for reproducibility;
  recoConsts::instance()->set_IntFlag("RANDOMSEED", 12345);

  // Particle Generator Setup;
  {
    auto gen = new PHG4ParticleGenerator();
    gen->set_name("pi+");           	  // geantino, pi-, pi+, mu-, mu+, e-., e+, proton, ... 
    gen->set_vtx(0,0,0);		  // Vertex generation range
    gen->set_mom_range(0,10);	  // Momentum generation range in GeV/c
    gen->set_eta_range(0,4);		  // Detector coverage around theta
    gen->set_phi_range(0.,2.*TMath::Pi());
    se->registerSubsystem(gen); 
  }

  // fun4all Geant4 wrapper;
  PHG4Reco* g4Reco = new PHG4Reco();

  // BeAST magnetic field;
  g4Reco->set_field_map(string(getenv("CALIBRATIONROOT")) + string("/Field/Map/mfield.4col.dat"), PHFieldConfig::kFieldBeast);

  // EicRoot media import; neither bound to EicToyModel nor to a particular EicRoot detector;
  EicGeoParData::ImportMediaFile("../../examples/eicroot/media.geo");

#ifdef _VST_
  
  // EicRoot vertex tracker; be aware: "VST" will also become a SuperDetector name;
  
  auto vst = new EicRootVstSubsystem("VST");
  {
    vst->SetGeometryType(EicGeoParData::NoStructure);
    vst->SetActive(true);

    // Barrel layers; hits belonging to these layers will be labeled internally
    // according to the sequence of these calls;
    {
      auto ibcell = new MapsMimosaAssembly();
      // See other MapsMimosaAssembly class POD entries in MapsMimosaAssembly.h;
      ibcell->SetDoubleVariable("mAssemblyBaseWidth", 17.5 * etm::mm);
  
      // Compose barrel layers; parameters are:
      //  - cell assembly type;
      //  - number of staves in this layer;
      //  - number of chips in a stave;
      //  - chip center installation radius;
      //  - additional stave slope around beam line direction; [degree];
      //  - layer rotation around beam axis "as a whole"; [degree];
      vst->AddBarrelLayer(ibcell, 1*3*12,  1*9, 0.1*3*23.4 * etm::mm, 12.0, 0.0);
      vst->AddBarrelLayer(ibcell, 2*3*12,  1*9, 0.2*3*23.4 * etm::mm, 12.0, 0.0);
      vst->AddBarrelLayer(ibcell, 3*3*12,  2*9, 0.3*3*23.4 * etm::mm, 12.0, 0.0);
      vst->AddBarrelLayer(ibcell, 4*3*12,  2*9, 0.4*3*23.4 * etm::mm, 12.0, 0.0);
    }

    g4Reco->registerSubsystem(vst);
  }
#endif

#ifdef _MMT_
  
  // EicRoot micromegas central tracker barrels;
  
  auto mmt = new EicRootMuMegasSubsystem("MMT");
  {
    mmt->SetActive(true);

    {
      auto layer = new MuMegasLayer();
      // See other MuMegasLayer class POD entries in MuMegasGeoParData.h;
      layer->SetDoubleVariable("mOuterFrameWidth", 20 * etm::mm);
      
      // Compose barrel layers; parameters are: 
      //   - layer description (obviously can mix different geometries);
      //   - length along Z;
      //   - segmentation in Z;
      //   - radius;
      //   - segmentation in phi;
      //   - Z offset from 0.0 (default);
      //   - azimuthal rotation from 0.0 (default);
      mmt->AddBarrel(layer, 600 * etm::mm, 2, 300 * etm::mm, 3, 0.0, 0.0);
      mmt->AddBarrel(layer, 600 * etm::mm, 3, 400 * etm::mm, 4, 0.0, 0.0);
      //mmt->AddBarrel(layer, 2400 * etm::mm, 3, 600 * etm::mm, 4, 0.0, 0.0);
      
      mmt->SetTransparency(50);
    }

    g4Reco->registerSubsystem(mmt);
  }
  
#endif
  
#ifdef _BEAMPIPE_
  // Central pipe dimension
  // Extracted via mechanical model: Detector chamber 3-20-20
  // directly implimenting the central Be section in G4 cylinder for max speed simulation in the detector region.
  // The jointer lip structure of the pipe R = 3.2cm x L=5mm is ignored here
  const double be_pipe_radius = 3.1000;
  const double be_pipe_thickness = 3.1762 - be_pipe_radius;  // 760 um for sPHENIX
  const double be_pipe_length_plus = 66.8;                   // +z beam pipe extend.
  const double be_pipe_length_neg = -79.8;                   // -z beam pipe extend.
  bool overlapcheck = false;
  bool absorberactive = false;
	
  // process pipe extentions?
  bool use_forward_pipes = false;
#ifdef _USE_FORWARD_PIPES_
  use_forward_pipes = true;
#endif
  const bool do_pipe_hadron_forward_extension = use_forward_pipes && true;
  const bool do_pipe_electron_forward_extension = use_forward_pipes && true;

  const double be_pipe_length = be_pipe_length_plus - be_pipe_length_neg;  // pipe length
  const double be_pipe_center = 0.5 * (be_pipe_length_plus + be_pipe_length_neg);
	
  // mid-rapidity beryillium pipe
  PHG4CylinderSubsystem* cyl = new PHG4CylinderSubsystem("VAC_BE_PIPE", 0);
  cyl->set_double_param("radius", 0.0);
  cyl->set_int_param("lengthviarapidity", 0);
  cyl->set_double_param("length", be_pipe_length);
  cyl->set_double_param("place_z", be_pipe_center);
  cyl->set_string_param("material", "G4_Galactic");
  cyl->set_double_param("thickness", be_pipe_radius);
  cyl->SuperDetector("PIPE");
  cyl->OverlapCheck(overlapcheck);
  //cyl->set_color(0, 0, 0, 0); //To make vacuum invisible
  if (absorberactive) cyl->SetActive();
  g4Reco->registerSubsystem(cyl);

  cyl = new PHG4CylinderSubsystem("BE_PIPE", 1);
  cyl->set_double_param("radius", be_pipe_radius);
  cyl->set_int_param("lengthviarapidity", 0);
  cyl->set_double_param("length", be_pipe_length);
  cyl->set_double_param("place_z", be_pipe_center);
  cyl->set_string_param("material", "G4_Be");
  cyl->set_double_param("thickness", be_pipe_thickness);
  cyl->SuperDetector("PIPE");
  cyl->OverlapCheck(overlapcheck);
  //cyl->set_color(192.0/255.0, 192.0/255.0, 192.0/255.0, 1); //To give beryllium greyish colour
  if (absorberactive) cyl->SetActive();
  g4Reco->registerSubsystem(cyl);
	
  if (do_pipe_electron_forward_extension)
    {
      PHG4GDMLSubsystem* gdml = new PHG4GDMLSubsystem("ElectronForwardEnvelope");
      gdml->set_string_param("GDMPath", string(getenv("CALIBRATIONROOT")) + "/Beam/Detector chamber 3-20-20.G4Import.gdml");
      gdml->set_string_param("TopVolName", "ElectronForwardEnvelope");
      gdml->set_int_param("skip_DST_geometry_export", 1);  // do not export extended beam pipe as it is not supported by TGeo and outside Kalman filter acceptance
      gdml->OverlapCheck(1);
      g4Reco->registerSubsystem(gdml);
    }

  if (do_pipe_hadron_forward_extension)
    {
      PHG4GDMLSubsystem* gdml = new PHG4GDMLSubsystem("HadronForwardEnvelope");
      gdml->set_string_param("GDMPath", string(getenv("CALIBRATIONROOT")) + "/Beam/Detector chamber 3-20-20.G4Import.gdml");
      gdml->set_string_param("TopVolName", "HadronForwardEnvelope");
      gdml->set_int_param("skip_DST_geometry_export", 1);  // do not export extended beam pipe as it is not supported by TGeo and outside Kalman filter acceptance
      gdml->OverlapCheck(1);
      g4Reco->registerSubsystem(gdml);
    }
#endif

#ifdef _BARREL_
  // Loading barrel from gdml file
  GdmlImportDetectorSubsystem * svtPart = new GdmlImportDetectorSubsystem();
  svtPart->set_string_param("GDMPath","Vst_GDML_.gdml"); 
  svtPart->OverlapCheck(); //Doesn't do anything, it seems.
  svtPart->AddAssemblyVolume("VST");	// Barrel. ********** THIS IS THE NAME IN THE GDML FILE, FROM EICROOT ********
  svtPart->SuperDetector("SVT");
  //svtPart->SetActive();          // this saves hits in the MimosaCore volumes
  //svtPart->SetAbsorberActive();  // this saves hits in all volumes (in the absorber node)
  g4Reco->registerSubsystem(svtPart);
#endif

#ifdef _FORWARD_SILICON_DISKS_
  //Forward disks
  GdmlImportDetectorSubsystem * fstPart = new GdmlImportDetectorSubsystem();
  fstPart->set_string_param("GDMPath","Fst_GDML_.gdml"); 
  fstPart->OverlapCheck(); //Doesn't do anything, it seems.
  fstPart->AddAssemblyVolume("FST");	// Disks. ********** THIS IS THE NAME IN THE GDML FILE, FROM EICROOT ********
  fstPart->SuperDetector("FstDisks");
  //fstPart->SetActive();          // this saves hits in the MimosaCore volumes
  //fstPart->SetAbsorberActive();  // this saves hits in all volumes (in the absorber node)
  g4Reco->registerSubsystem(fstPart);
#endif
  
  
#ifdef _BACKWARD_SILICON_DISKS_
  //Backward disks
  GdmlImportDetectorSubsystem * bstPart = new GdmlImportDetectorSubsystem();
  bstPart->set_string_param("GDMPath","Bst_GDML_.gdml"); 
  bstPart->OverlapCheck(); //Doesn't do anything, it seems.
  bstPart->AddAssemblyVolume("BST");	// Disks. ********** THIS IS THE NAME IN THE GDML FILE, FROM EICROOT ********
  bstPart->SuperDetector("BstDisks");
  //bstPart->SetActive();          // this saves hits in the MimosaCore volumes
  //bstPart->SetAbsorberActive();  // this saves hits in all volumes (in the absorber node)
  g4Reco->registerSubsystem(bstPart);
#endif
  	
#ifdef _TPC_
  // time projection chamber layers ----------------------
  //This is just cylinder subsystems basically.     
  //What layer numbers to use? This will be the detector (group) number, right? So.. We can start at 60.                                              
  PHG4CylinderSubsystem* cylTPC;
  //cylTPC->CheckOverlap();

  double radius = G4TPC::inner_cage_radius;

  // inner field cage                                                                                                                                 
  cylTPC = new PHG4CylinderSubsystem("FIELDCAGE", 60);
  cylTPC->set_double_param("radius", radius);
  cylTPC->set_double_param("length", G4TPC::cage_length);
  cylTPC->set_string_param("material", "G4_KAPTON");
  cylTPC->set_double_param("thickness", G4TPC::cage_thickness);
  cylTPC->SuperDetector("FIELDCAGE");
  cylTPC->Verbosity(0);
  g4Reco->registerSubsystem(cylTPC);


  int n_tpc_layers[3] = {16, 16, 16};
  int tpc_layer_rphi_count[3] = {1152, 1536, 2304};

  double tpc_region_thickness[3] = {20., 20., 18.};
  // Active layers of the TPC (inner layers)                                                                                                          
  int nlayer = 60;
  for (int irange = 0; irange < 3; irange++)
    {
      double tpc_layer_thickness = tpc_region_thickness[irange] / n_tpc_layers[irange];  // thickness per layer                                       
      for (int ilayer = 0; ilayer < n_tpc_layers[irange]; ilayer++)
        {
          cylTPC = new PHG4CylinderSubsystem("TPC", nlayer);
          cylTPC->set_double_param("radius", radius);
          cylTPC->set_double_param("length", G4TPC::cage_length);
          cylTPC->set_string_param("material", G4TPC::tpcgas);
          cylTPC->set_double_param("thickness", tpc_layer_thickness - 0.01);
          cylTPC->SetActive();/////////////////////////////////////////////////////////sets the volume to be active in g4eval
          cylTPC->SuperDetector("TPC");
          g4Reco->registerSubsystem(cylTPC);

          radius += tpc_layer_thickness;
          nlayer++;
        }
    }

  // outer field cage                                                                                                                                 
  cylTPC = new PHG4CylinderSubsystem("FIELDCAGE", nlayer);
  cylTPC->set_double_param("radius", radius);
  cylTPC->set_int_param("lengthviarapidity", 0);
  cylTPC->set_double_param("length", G4TPC::cage_length);
  cylTPC->set_string_param("material", "G4_KAPTON");
  cylTPC->set_double_param("thickness", G4TPC::cage_thickness);  // Kapton X_0 = 28.6 cm                                                              
  cylTPC->SuperDetector("FIELDCAGE");
  g4Reco->registerSubsystem(cylTPC);
#endif //End of TPC definition                                                                                                                        

#ifdef _TPC_ENDCAPS_
  bool AbsorberActive = false;
  
  PHG4TpcEndCapSubsystem* tpc_endcap = new PHG4TpcEndCapSubsystem("TPC_ENDCAP");
  tpc_endcap->SuperDetector("TPC_ENDCAP");

  if (AbsorberActive) {
    tpc_endcap->SetActive();
  }
  tpc_endcap->CheckOverlap();

  //  tpc_endcap->set_int_param("construction_verbosity", 2);                                                                                         

  g4Reco->registerSubsystem(tpc_endcap);
#endif

  int dz = 20;
#ifdef _GEMS_
  
  // Forward GEM tracker module(s);
  auto fgt = new EicRootGemSubsystem("FGT");
  
  {
    fgt->SetActive(true);
    fgt->CheckOverlap();
    fgt->SetTGeoGeometryCheckPrecision(0.000001 * etm::um);

    {
      auto sbs = new GemModule();
      // See other GemModule class data in GemGeoParData.h;
      sbs->SetDoubleVariable("mFrameBottomEdgeWidth", 30 * etm::mm);

      // Compose sectors; parameters are: 
      //   - layer description (obviously can mix different geometries);
      //   - azimuthal segmentation;
      //   - gas volume center radius;
      //   - Z offset from 0.0 (default);
      //   - azimuthal rotation from 0.0 (default);      
      //fgt->AddWheel(sbs, 12, 420 * etm::mm, (1000+dz) * etm::mm, 0);
      //fgt->AddWheel(sbs, 12, 420 * etm::mm, (1100+dz) * etm::mm, 0);
      fgt->AddWheel(sbs, 12, 420 * etm::mm, (1200+dz) * etm::mm, 0);
      fgt->AddWheel(sbs, 12, 420 * etm::mm, (1300+dz) * etm::mm, 0);
      //fgt->AddWheel(sbs, 12, 420 * etm::mm, 1400 * etm::mm, 0);
      //fgt->AddWheel(sbs, 12, 420 * etm::mm, 1500 * etm::mm, 0);
      //fgt->AddWheel(sbs, 12, 420 * etm::mm, -1300 * etm::mm, 0);
    }

    g4Reco->registerSubsystem(fgt);
  }
/*    
  GdmlImportDetectorSubsystem * gemPart = new GdmlImportDetectorSubsystem();
  gemPart->set_string_param("","GemGeoParData.h"); 
  gemPart->AddAssemblyVolume("GEM");	// Barrel. ********** THIS IS THE NAME IN THE GDML FILE, FROM EICROOT ********
  gemPart->SuperDetector("GEM_disk");
  gemPart->SetActive();          // this saves hits in the MimosaCore volumes
  gemPart->SetAbsorberActive();  // this saves hits in all volumes (in the absorber node)
  g4Reco->registerSubsystem(gemPart);
*/    
#endif
  
    
  // Truth information;
//  g4Reco->registerSubsystem(new PHG4TruthSubsystem());
  //se->registerSubsystem(g4Reco);

  //Register the truth hits
  PHG4TruthSubsystem *truth = new PHG4TruthSubsystem();
  g4Reco->registerSubsystem(truth);


  se->registerSubsystem(g4Reco);


#ifdef _QT_DISPLAY_
  g4Reco->InitRun(se->topNode());
  g4Reco->ApplyDisplayAction();
  g4Reco->StartGui();
#else
  
  // Ideal track finder and Kalman filter;
  {
    auto kalman = new PHG4TrackFastSim("PHG4TrackFastSim");
    
    kalman->set_use_vertex_in_fitting(false);
    //===================================
    kalman->set_sub_top_node_name("SVTX");
    kalman->set_trackmap_out_name("SvtxTrackMap");
    //======================================	
    double pixelResolution = 10e-4/sqrt(12);

#ifdef _VST_    
    // Silicon tracker hits;
    
    kalman->add_phg4hits(vst->GetG4HitName(),		// const std::string& phg4hitsNames
			 PHG4TrackFastSim::Cylinder,	// const DETECTOR_TYPE phg4dettype
			 999.,				// radial-resolution [cm] (this number is not used in cylindrical geometry)
			 // 20e-4/sqrt(12) cm = 5.8e-4 cm, to emulate 20x20 um pixels;
			 5.8e-4,			// azimuthal (arc-length) resolution [cm]
			 5.8e-4,			// longitudinal (z) resolution [cm]
			 1,				// efficiency (fraction)
			 0);				// hit noise
#endif

#ifdef _MMT_
    
    // MM tracker hits;
    
      kalman->add_phg4hits(mmt->GetG4HitName(),		// const std::string& phg4hitsNames
      PHG4TrackFastSim::Cylinder,	// const DETECTOR_TYPE phg4dettype
      999.,				// radial-resolution [cm] (this number is not used in cylindrical geometry)
      // Say 50um resolution?; [cm];
      50e-4,			        // azimuthal (arc-length) resolution [cm]
      50e-4,		        	// longitudinal (z) resolution [cm]
      1,				// efficiency (fraction)
      0);				// hit noise
    
#endif
      /*
#ifdef _BARREL_
    //SVT
    //This adds the hit nodes for EACH LAYER OF THE DETECTOR.
    for (int i = 10; i < 10+_NO_OF_BARREL_LAYERS_; i++) {
      std::string nodeName = "G4HIT_SVT_" + std::to_string(i); // ****** THIS NEEDS TO MATCH THE NAME FROM GetName() in GdmlImportDetector.cc ******
      kalman->add_phg4hits(
			   nodeName,				// const std::string& phg4hitsNames
			   PHG4TrackFastSim::Cylinder,		// const DETECTOR_TYPE phg4dettype
			   999.,					// radial-resolution [cm] (this number is not used in cylindrical geometry)
			   pixelResolution,					// azimuthal (arc-length) resolution [cm]
			   pixelResolution,					// longitudinal (z) resolution [cm]
			   1,					// efficiency (fraction)
			   0					// hit noise
			   );
    }
#endif
      */
       
#ifdef _TPC_
    kalman->add_phg4hits(
                         "G4HIT_TPC",                //      const std::string& phg4hitsNames,                                                        
                         PHG4TrackFastSim::Cylinder,  //      const DETECTOR_TYPE phg4dettype,      
	                 1,                           //      const float radres,                                                                     
                         200e-4,                      //      const float phires,                                                                     
                         500e-4,                      //      const float lonres,                                                                     
                         1,                           //      const float eff,                                                                        
                         0                            //      const float noise                                                                       
                         );
#endif
   
    /*    
#ifdef _FORWARD_SILICON_DISKS_
    //Forward silicon disks
    for (int i = 40; i < 40+_NO_OF_FORWARD_DISKS_ ; i++) {
      std::string nodeName = "G4HIT_FstDisks_" + std::to_string(i); // ****** THIS NEEDS TO MATCH THE NAME FROM GetName() in GdmlImportDetector.cc ******
      //Comes from SuperDetector name when adding it to Geant4.
      kalman->add_phg4hits(
			   nodeName,				// const std::string& phg4hitsNames
			   PHG4TrackFastSim::Vertical_Plane,		// const DETECTOR_TYPE phg4dettype
			   pixelResolution,					// radial-resolution [cm] 
			   pixelResolution,					// azimuthal (arc-length) resolution [cm]
			   999.0,					// longitudinal (z) resolution [cm](not used for vertical plane)
			   1,					// efficiency (fraction)
			   0					// hit noise
			   );
    }
#endif
    
    
#ifdef _BACKWARD_SILICON_DISKS_
    //Backward silicon disks
    for (int i = 50; i < 50+_NO_OF_BACKWARD_DISKS_ ; i++) {
      std::string nodeName = "G4HIT_BstDisks_" + std::to_string(i); // ****** THIS NEEDS TO MATCH THE NAME FROM GetName() in GdmlImportDetector.cc ******
      //Comes from SuperDetector name when adding it to Geant4.
      kalman->add_phg4hits(
			   nodeName,				// const std::string& phg4hitsNames
			   PHG4TrackFastSim::Vertical_Plane,		// const DETECTOR_TYPE phg4dettype
			   pixelResolution,					// radial-resolution [cm] 
			   pixelResolution,					// azimuthal (arc-length) resolution [cm]
			   999.0,					// longitudinal (z) resolution [cm] (not used for vertical plane)
			   1,					// efficiency (fraction)
			   0					// hit noise
			   );
    }
#endif
 	*/

#ifdef _GEMS_
    // GEM tracker hits; should work;
    kalman->add_phg4hits(fgt->GetG4HitName(),		// const std::string& phg4hitsNames
			 PHG4TrackFastSim::Vertical_Plane,	// const DETECTOR_TYPE phg4dettype
			 999.,				// radial-resolution [cm] (this number is not used in cylindrical geometry)
			 // Say 70um resolution?; [cm];
			 70e-4,			        // azimuthal (arc-length) resolution [cm]
			 70e-4,		        	// longitudinal (z) resolution [cm]
			 1,				// efficiency (fraction)
			 0);				// hit noise
    
#endif
    
    se->registerSubsystem(kalman);

    TrackFastSimEval *fast_sim_eval = new TrackFastSimEval("FastTrackingEval");
    //fast_sim_eval->set_filename(TString(outputFile)+B_label+"_FastTrackingEval.root");
    //fast_sim_eval->AddProjection(radialBH);
    //fast_sim_eval->AddProjection(forwardBH);
    //fast_sim_eval->AddProjection(backwardBH);
    se->registerSubsystem(fast_sim_eval);

    SimpleNtuple * hits = new SimpleNtuple("Hits");
    /*
#ifdef _BARREL_
    //SVT starts at 10, to allow for absorbers before
    hits->AddNode("ABSORBER_SVT",0); // hits in the passive volumes
    for (int i = 10; i < 10+_NO_OF_BARREL_LAYERS_; i++) { // hits in the  MimosaCore volumes
      std::string nodeName = "SVT_" + std::to_string(i); // ****** THIS NEEDS TO MATCH WHAT'S INPUT IN GdmlImportDetector.cc, stripped of the G4HIT, because currently SimpleNtuple adds that.
      // See SimpleNtuple.cc ******
      hits->AddNode(nodeName, i);
    }
#endif
    */
#ifdef _TPC_
    hits->AddNode("TPC", 60);
#endif
    
    /*
#ifdef _FORWARD_SILICON_DISKS_
    //Forward silicon disks
    hits->AddNode("ABSORBER_FstDisks",3); // hits in the passive volumes
    for (int i = 40; i < 40+_NO_OF_FORWARD_DISKS_; i++) { // hits in the  MimosaCore volumes
      std::string nodeName = "FstDisks_" + std::to_string(i); // ****** THIS NEEDS TO MATCH WHAT'S INPUT IN GdmlImportDetector.cc, stripped of the G4HIT, because currently SimpleNtuple adds that.
      // See SimpleNtuple.cc ******
      hits->AddNode(nodeName, i);
    }
#endif
    
    
#ifdef _BACKWARD_SILICON_DISKS_
    //Backward silicon disks
    hits->AddNode("ABSORBER_BstDisks",4); // hits in the passive volumes
    for (int i = 50; i < 50+_NO_OF_FORWARD_DISKS_; i++) { // hits in the  MimosaCore volumes
      std::string nodeName = "BstDisks_" + std::to_string(i); // ****** THIS NEEDS TO MATCH WHAT'S INPUT IN GdmlImportDetector.cc, stripped of the G4HIT, because currently SimpleNtuple adds that.
      // See SimpleNtuple.cc ******
      hits->AddNode(nodeName, i);
    }
#endif
    */
#ifdef _GEMS_
    hits->AddNode("Tracker_GEMS",70);
    //std::cout << "GEMS" << std::endl;
#endif
    /*
#ifdef _VST_
    hits->AddNode("Tracker_VST",70);
#endif

#ifdef _MMT_
    hits->AddNode("Tracker_MMT",80);
#endif
    */
    cout << "HHHHHHHHHHHHHHHH" << endl;
    se->registerSubsystem(hits);


  }
  cout << "DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDddDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD" << endl;

  /*
   const std::string dst_name = std::string(outputFile)+"_gdmlimporter.root";
  //Fun4AllDstOutputManager *out = new Fun4AllDstOutputManager("DSTOUT",TString(outputFile)+"_gdmlimporter.root");
  Fun4AllDstOutputManager *out = new Fun4AllDstOutputManager("DSTOUT",dst_name);
  out->Verbosity(10);
  //Let's not output full DST tree right now
  //se->registerOutputManager(out);
	
  
  if (particle_gen != 1 && particle_gen != 2) {
		
    Fun4AllInputManager *in2 = new Fun4AllHepMCInputManager("DSTIN");
    se->registerInputManager(in2);
    se->fileopen(in2->Name().c_str(), "../../eventFiles/ggm2ccbar_28_1_100000.dat"); //HepMC event file
  } else {
    Fun4AllInputManager *in = new Fun4AllDummyInputManager("JADE");
    se->registerInputManager(in);
    }
  */

  
  // User analysis code: just a single dp/p histogram;
  se->registerSubsystem(new TrackFastSimEval());
  cout << "EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE" << endl;
  // Run it all, eventually;
  se->run(nEvents);

  cout << "GGGGGGGGGGGGGGG" << endl;
  se->End();
#endif
  
  cout << "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF" << endl;
  delete se;
  gSystem->Exit(0);
} // Fun4All_G4_Sandbox()
