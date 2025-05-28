//
// Template for this file generated with eicmkplugin.py
//

#include <algorithm>
#include <bitset>
#include <spdlog/spdlog.h>

#include <JANA/JEventProcessorSequentialRoot.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <THashList.h>

#include <edm4hep/MCParticle.h>
#include <edm4hep/SimTrackerHit.h>
#include <edm4hep/SimCalorimeterHit.h>
//#include <edm4hep/RawCalorimeterHit.h>

//#include <edm4eic/Cluster.h>
//#include <edm4eic/RawCalorimeterHit.h>
//#include <edm4eic/ProtoCluster.h>

#include <services/geometry/dd4hep/DD4hep_service.h>

#include "../../utilities/constants.h"

using namespace std;

typedef std::tuple<double,double,double> MyHit;

struct TrackClass {
  double charge;
  double X0;
  double Y0;
  double slopeX;
  double slopeY;
  double theta;
  double phi;
  double Chi2;
};

struct TreeTrackClass {
  vector<double> X0_e;
  vector<double> Y0_e;
  vector<double> slopeX_e;
  vector<double> slopeY_e;

  vector<double> X0_p;
  vector<double> Y0_p;
  vector<double> slopeX_p;
  vector<double> slopeY_p;
};

class analyzeLumiHits: public JEventProcessorSequentialRoot {
  private:

    static const int maxModules = 3;
    static const int maxSectors = 2;
   
    int Nhits_min = 1; // 15

    // spectrometer dimensions/placements in cm
    // DXY, and DZ stand for FULL widths
    double LumiSpecMag_Z = -5600;
    double LumiSpecMag_DZ = 78;
    double LumiSpecCAL_Z = -6500;
    double LumiSpecCALTower_DZ = 20;
    double LumiSpecCAL_DXY = 20;
    double LumiSpecCAL_FiveSigma = 6.9;
    double LumiConverter_Z = LumiSpecMag_Z + LumiSpecMag_DZ/2.0;
    double LumiSpecMagEnd_Z = LumiSpecMag_Z - LumiSpecMag_DZ/2.0;
    double LumiConverterCut_DXY = 6;

    double pT = 0.117; // GeV. 0.3*B(T)*dZ(m)
    // cyclotron radius = speed / cyclotron frequency -> p/(q*B) = E/(c*q*B) in ultrarelativistic limit
    double RmagPreFactor = 667.079; // (J/GeV)/(c * q * B), multiply this by E in GeV to get R in cm

    double LumiSpecTracker_Z1;
    double LumiSpecTracker_Z2;
    double LumiSpecTracker_Z3;
    vector<double> Tracker_Zs;
    double Tracker_meanZ;

    double Tracker_pixelSize = 0.005; // cm
    //maximal reduced chi^2 for tracks
    double max_chi2ndf = 0.01;
    double Tracker_sigma = 0.39; // cm from reconstructed photon origins decaying into 2 electrons.
    
    double Einput;
    int Ntrackers;

    // All histograms stored in a THashList
    THashList *gHistList;
   
    TTree *tree_MC;
    TTree *tree_CAL;
    TTree *tree_Trk;

    //MC Hits 
    double Egen;
    double XVertex;
    double YVertex;
    double ZVertex;

    //CALorimeter hits
    static const int CALNhits_max = 99999; //static and const while initilising array size is must
    int CALNhits;
    double CALE_hit[CALNhits_max];
    double CALx_hit[CALNhits_max];
    double CALy_hit[CALNhits_max];
    double CALz_hit[CALNhits_max];
    int CALsec_id[CALNhits_max];
    int CALblock[CALNhits_max];
    int CALlayer_id[CALNhits_max];
    int CALmod_id[CALNhits_max];
    int CALfiber_x_id[CALNhits_max];
    int CALfiber_y_id[CALNhits_max];
    
    //TRacKer hits
    static const int TrkNhits_max = 20;
    int TrkNhits;
    //double TrkE_hit[TrkNhits_max]; Tracker_hits doesn't store value of energy
    double Trkx_hit[TrkNhits_max];
    double Trky_hit[TrkNhits_max];
    double Trkz_hit[TrkNhits_max];
    int Trksec_id[TrkNhits_max];
    int Trkmod_id[TrkNhits_max];
    int Trktype[TrkNhits_max];

    TreeTrackClass treeTracks;
    
    // Data objects we will need from JANA e.g.
    PrefetchT<edm4hep::MCParticle> MCParticles		= {this, "MCParticles"};
    PrefetchT<edm4hep::SimCalorimeterHit> CAL_hits      = {this, "EcalLumiSpecHits"};
    //PrefetchT<edm4hep::RawCalorimeterHit> CAL_adc       = {this, "EcalLumiSpecRawHits"};
    //PrefetchT<edm4eic::CalorimeterHit> CAL_rechits      = {this, "EcalLumiSpecRecHits"};
    //PrefetchT<edm4eic::ProtoCluster> CAL_protoClusters  = {this, "EcalLumiSpecIslandProtoClusters"};
    //PrefetchT<edm4eic::Cluster> CAL_clusters            = {this, "EcalLumiSpecClusters"};

    PrefetchT<edm4hep::SimTrackerHit> Tracker_hits      = {this, "LumiSpecTrackerHits"};

  public:
    analyzeLumiHits() { SetTypeName(NAME_OF_THIS); }
    //Main functions
    void InitWithGlobalRootLock() override;
    void ProcessSequential(const std::shared_ptr<const JEvent>& event) override;
    void FinishWithGlobalRootLock() override;
    //User Added Functions
    bool PixelOverlap( MyHit hit, vector<MyHit> trackSet );
    void AssembleTracks( vector<TrackClass> *tracks, vector<vector<MyHit>> hitSet );
    double TrackerErec( double slopeY );
    double DeltaYmagnet( double E, double charge );
    double XatConverter( TrackClass track );
    double YatConverter( TrackClass track );
    double GetPairMass( TrackClass top, TrackClass bot );
    void PrintTrackInfo( vector<TrackClass> topTracks, vector<TrackClass> botTracks );
  protected:

    std::shared_ptr<DD4hep_service> m_geoSvc;
};


