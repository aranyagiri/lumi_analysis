
#include "analyzeLumiHits.h"
#include <services/rootfile/RootFile_service.h>

// The following just makes this a JANA plugin
extern "C" {
  void InitPlugin(JApplication *app) {
    InitJANAPlugin(app);
    app->Add(new analyzeLumiHits);
  }
}

//-------------------------------------------
// InitWithGlobalRootLock
//-------------------------------------------
void analyzeLumiHits::InitWithGlobalRootLock(){

  auto rootfile_svc = GetApplication()->GetService<RootFile_service>();
  auto rootfile = rootfile_svc->GetHistFile();
  
  //MC branch
  tree_MC = new TTree("tree_MC","tree_MC");
  tree_MC->Branch("Egen", &Egen);
  tree_MC->Branch("XVertex", &XVertex);
  tree_MC->Branch("YVertex", &YVertex);
  tree_MC->Branch("ZVertex", &ZVertex);
  
  //Calorimeter branches
  tree_CAL = new TTree("tree_CAL","tree_CAL");
  tree_CAL->Branch("CALNhits",&CALNhits);
  tree_CAL->Branch("CALE", CALE_hit, "CALE_hit[CALNhits]/D");
  tree_CAL->Branch("CALx", CALx_hit, "CALx_hit[CALNhits]/D");
  tree_CAL->Branch("CALy", CALy_hit, "CALy_hits[CALNhits]/D");
  tree_CAL->Branch("CALz", CALz_hit, "CALz_hits[CALNhits]/D"); 
  tree_CAL->Branch("CALsec_id", CALsec_id, "CALsec_id[CALNhits]/I");
  tree_CAL->Branch("CALlayer_id", CALlayer_id, "CALlayer_id[CALNhits]/I");
  tree_CAL->Branch("CALblock", CALblock, "CALblock[CALNhits]/I");
  tree_CAL->Branch("CALmod_id", CALmod_id, "CALmod_id[CALNhits]/I");
  tree_CAL->Branch("CALfiber_x_id", CALfiber_x_id, "CALfiber_x_id[CALNhits]/I");
  tree_CAL->Branch("CALfiber_y_id", CALfiber_y_id, "CALfiber_y_id[CALNhits]/I");

  
  //Trcaker branches
  tree_Trk = new TTree("tree_Trk","tree_Trk");
  tree_Trk->Branch("TrkNhits",&TrkNhits);
  //tree_Trk->Branch("TrkE", TrkE_hit, "TrkE_hit[Nhits]/D");
  tree_Trk->Branch("Trkx", Trkx_hit, "Trkx_hit[Nhits]/D");
  tree_Trk->Branch("Trky", Trky_hit, "Trky_hits[Nhits]/D");
  tree_Trk->Branch("Trkz", Trkz_hit, "Trkz_hits[Nhits]/D"); 
  tree_Trk->Branch("Trksec_id", Trksec_id, "Trksec_id[Nhits]/I");
  tree_Trk->Branch("Trkmod_id", Trkmod_id, "Trkmod_id[Nhits]/I");
  tree_Trk->Branch("Trktype", Trktype, "Trktype[Nhits]/I");

}
//-------------------------------------------
// ProcessSequential
//-------------------------------------------
void analyzeLumiHits::ProcessSequential(const std::shared_ptr<const JEvent>& event) {

  //cout<<"New Event"<<endl;

  auto app = GetApplication();
  m_geoSvc = app->template GetService<DD4hep_service>();

  //----------------MC Analysis--------------------------------------------- 
  for( auto particle : MCParticles() ) {
     
    edm4hep::Vector3f p = particle->getMomentum();
    edm4hep::Vector3d v = particle->getVertex(); // Units of mm !!

    if( particle->getPDG() == +11 && particle->getGeneratorStatus() == 1 ) {
      Egen = particle->getEnergy();
      XVertex = v.x;
      YVertex = v.y;
      ZVertex = v.z;

    }
    else if( particle->getPDG() == -11 && particle->getGeneratorStatus() == 1 ) {
      Egen = particle->getEnergy();
      XVertex = v.x;
      YVertex = v.y;
      ZVertex = v.z;
    }
    else {}
  }

  tree_MC->Fill();

  //-------------------------------------------------------------------------
  
  Einput = 0;
  Ntrackers = 0;
  app->GetParameter("analyzeLumiHits:Egen", Einput);
  app->GetParameter("analyzeLumiHits:Ntrackers", Ntrackers);
  
  ///////////////////////////////////////////////////////////////////////
  			//Calroimeter//
  //////////////////////////////////////////////////////////////////////
  map<string, int> CALfield_idx_Map{ {"sector", 0}, {"layer", 0}, {"module", 0}, {"block", 0}, {"fiber_x", 0}, {"fiber_y", 0}};

  int hit_num = 0;

  for( auto hit : CAL_hits() ) {

	  const auto id = hit->getCellID();
	  auto id_dec   = m_geoSvc->detector()->readout( "EcalLumiSpecHits" ).idSpec().decoder();

	  // field of readout fields
	  vector<dd4hep::BitFieldElement> hitFields = id_dec->fields();

	  // try to find the expected fields and store field index
	  for( auto field : hitFields ) {
		  if( CALfield_idx_Map.find( field.name() ) != CALfield_idx_Map.end() ) {
			  CALfield_idx_Map[ field.name() ] = id_dec->index( field.name() );
		  }
	  }

	  // look up sector,module,fiber... id of this hit 
	  CALsec_id[hit_num] 		= (int) id_dec->get( id, CALfield_idx_Map["sector"] ); // Top (0) and Bottom (1)
	  CALblock[hit_num]  		= (int) id_dec->get( id, CALfield_idx_Map["block"] );
	  CALlayer_id[hit_num]		= (int) id_dec->get( id, CALfield_idx_Map["layer"] );
	  CALmod_id[hit_num] 		= (int) id_dec->get( id, CALfield_idx_Map["module"] ); // 10x10 Matrix of bars
	  CALfiber_x_id[hit_num]  	= (int) id_dec->get( id, CALfield_idx_Map["fiber_x"] );
	  CALfiber_y_id[hit_num] 	= (int) id_dec->get( id, CALfield_idx_Map["fiber_y"] );

	  edm4hep::Vector3f vec = hit->getPosition();// mm
	  double x = vec.x;
	  double y = vec.y;
	  double z = vec.z;

	  CALE_hit[hit_num] 	= hit->getEnergy();
	  CALx_hit[hit_num] 	= x;
	  CALy_hit[hit_num] 	= y;
	  CALz_hit[hit_num] 	= z;

	  hit_num++;

  }//CAL hits close

  CALNhits = hit_num;
  tree_CAL->Fill();

  ///////////////////////////////////////////////////////////////////////
				//Tracker//
  ///////////////////////////////////////////////////////////////////////
  map<string, int> Trackerfield_idx_Map{ {"sector", 0}, {"module", 0}};

	hit_num = 0;
	for( auto hit : Tracker_hits() ){

		const auto id = hit->getCellID();
		auto id_dec   = m_geoSvc->detector()->readout( "LumiSpecTrackerHits" ).idSpec().decoder();

		bool primary = (hit->getQuality() == 0) ? true : false; // 1 << 30 if produced by secondary (edm4hep docs)

		// field of readout fields
		vector<dd4hep::BitFieldElement> hitFields = id_dec->fields();

		// try to find the expected fields and store field index
		for( auto field : hitFields ) {
			if( Trackerfield_idx_Map.find( field.name() ) != Trackerfield_idx_Map.end() ) {
				Trackerfield_idx_Map[ field.name() ] = id_dec->index( field.name() );
			}
		}

		// look up sector,module,fiber... id of this hit 
		int sec 	= (int) id_dec->get( id, Trackerfield_idx_Map["sector"] ); // top or bottom
		int mod 	= (int) id_dec->get( id, Trackerfield_idx_Map["module"] ); // layer, closet to furthest from IP

		//edm4hep::Vector3f vec = hit->getPosition();// mm
		double x = hit->x();
		double y = hit->y();
		double z = hit->z();

		//TrkE_hit[hit_num] 	= hit->getEnergy(); Tracker container doesn't contain energy information.
		Trkx_hit[hit_num] 	= x;
		Trky_hit[hit_num] 	= y;
		Trkz_hit[hit_num] 	= z;
		Trksec_id[hit_num]	= sec;
		Trkmod_id[hit_num] 	= mod;
		Trktype[hit_num] 	= (int) primary;
		hit_num++;

	}//tracker hits

	TrkNhits = hit_num;
  	tree_Trk->Fill();

} // End of the Sequential Process Function

//-------------------------------------------------------------------------
double analyzeLumiHits::TrackerErec( double slopeY ) {

  double sinTheta = fabs( sin( atan(slopeY) ) );
 
  if( sinTheta == 0 ) { return 0.0; }

  double E = pT / sinTheta;
  
  return E;
}

//-------------------------------------------------------------------------
bool analyzeLumiHits::PixelOverlap( MyHit hit, vector<MyHit> trackSet ) {

  for( auto el : trackSet ) {
    double delta = pow( std::get<0>(hit) - std::get<0>(el), 2);
    delta += pow( std::get<1>(hit) - std::get<1>(el), 2);

    if( sqrt(delta) < Tracker_pixelSize ) { return true; }
  }

  return false;
}

//-------------------------------------------------------------------------
void analyzeLumiHits::AssembleTracks( vector<TrackClass> *tracks, vector<vector<MyHit>> hitSet ) {

  // units of hits are in cm

  for( auto hit1 : hitSet[0] ) { // tracker hit closest to IP

    for( auto hit2 : hitSet[1] ) {

      for( auto hit3 : hitSet[2] ) {

        double sum_x = std::get<0>(hit1) + std::get<0>(hit2) + std::get<0>(hit3);
        double sum_y = std::get<1>(hit1) + std::get<1>(hit2) + std::get<1>(hit3);
        double sum_z = std::get<2>(hit1) + std::get<2>(hit2) + std::get<2>(hit3);
        double sum_zx = std::get<2>(hit1)*std::get<0>(hit1) + std::get<2>(hit2)*std::get<0>(hit2) + std::get<2>(hit3)*std::get<0>(hit3);
        double sum_zy = std::get<2>(hit1)*std::get<1>(hit1) + std::get<2>(hit2)*std::get<1>(hit2) + std::get<2>(hit3)*std::get<1>(hit3);
        double sum_zz = pow(std::get<2>(hit1), 2) + pow(std::get<2>(hit2), 2) + pow(std::get<2>(hit3), 2);

        // Least squares regression algorithm: assumes equal Gaussian errors with each hit point
        // y = m_y * z + y0
        // x = m_x * z + x0
        // m_y = Sum( (z_i - <z>)*(y_i - <y>) ) / Sum( (z_i - <z>)^2 )
        //     = ( N*Sum(z_i*y_i) - Sum(z_i)*Sum(y_i) ) / ( N*Sum(z_i^2) - Sum(z_i)^2 )
        // y0  = <y> - m_y*<z>
        //     = ( Sum(y_i) - m_y*Sum(z_i) ) / N

        TrackClass track;
        if( sum_y > 0 ) { track.charge = -1; }// electrons go to top CAL (B in +x direction)
        else { track.charge = +1; }
        track.slopeX = (3. * sum_zx - sum_z * sum_x) / (3. * sum_zz - sum_z * sum_z);
        track.slopeY = (3. * sum_zy - sum_z * sum_y) / (3. * sum_zz - sum_z * sum_z);
        track.X0 = (sum_x - track.slopeX * sum_z) / 3.0;
        track.Y0 = (sum_y - track.slopeY * sum_z) / 3.0;
        // calculate polar (theta) and azimuthal angles (phi)
        double delta_x = std::get<0>(hit3) - std::get<0>(hit1);
        double delta_y = std::get<1>(hit3) - std::get<1>(hit1);
        double delta_z = std::get<2>(hit3) - std::get<2>(hit1);
        track.phi = atan2( delta_y, delta_x );
        track.theta = TMath::Pi() + atan( sqrt( pow(delta_x, 2) + pow(delta_y, 2) ) / delta_z );

        // Chi2
        double Chi2 = 0;
        auto hit_group = {hit1, hit2, hit3};
        for( auto hit : hit_group ) {
          double DeltaX = pow(std::get<0>(hit) - (track.X0 + track.slopeX * std::get<2>(hit)), 2);
          double DeltaY = pow(std::get<1>(hit) - (track.Y0 + track.slopeY * std::get<2>(hit)), 2);
          Chi2 += DeltaX + DeltaY;
        }
        track.Chi2 = Chi2;

        tracks->push_back( track );
      }
    }
  }
}

//-------------------------------------------------------------------------
// returns y deflection of ultrarel electrons in magnet region assuming primordial py = 0
double analyzeLumiHits::DeltaYmagnet( double E, double charge ) {

  //return 0;

  double R = E * RmagPreFactor; // cyclotron radius of curvature
  double dy = R - sqrt( R*R - pow(LumiSpecMag_DZ,2) );

  // electrons go to the top CAL (B in +x direction), positrons to bottom CAL
  return -charge * dy; // cm
}

//-------------------------------------------------------------------------
double analyzeLumiHits::XatConverter( TrackClass track ) {

  double x_c = (track.slopeX * LumiConverter_Z + track.X0);
  
  return x_c;
}

//-------------------------------------------------------------------------
double analyzeLumiHits::YatConverter( TrackClass track ) {
  double E = TrackerErec( track.slopeY );
  double y_c = (track.slopeY * LumiSpecMagEnd_Z + track.Y0) - DeltaYmagnet( E, track.charge );
  
  return y_c;
}

//-------------------------------------------
// FinishWithGlobalRootLock
//-------------------------------------------
void analyzeLumiHits::FinishWithGlobalRootLock() {

  // Do any final calculations here.
}

//-------------------------------------------------------------------------
void analyzeLumiHits::PrintTrackInfo( vector<TrackClass> topTracks, vector<TrackClass> botTracks) {
  cout<<"N good TopTracks: "<<topTracks.size()<<endl;
  for( auto el : topTracks ) {
    cout<<"theta: "<<el.theta<<"    phi: "<<el.phi<<"   slopeX: "<<el.slopeX<<"   slopeY: "<<el.slopeY<<"    Chi2: "<<el.Chi2<<endl;
  }

  cout<<"N good BotTracks: "<<botTracks.size()<<endl;
  for( auto el : botTracks ) {
    cout<<"theta: "<<el.theta<<"    phi: "<<el.phi<<"   slopeX: "<<el.slopeX<<"   slopeY: "<<el.slopeY<<"    Chi2: "<<el.Chi2<<endl;
  }
}

//-------------------------------------------------------------------------
double analyzeLumiHits::GetPairMass( TrackClass top, TrackClass bot ) {
  
  double Etop = TrackerErec( top.slopeY );
  double Ebot = TrackerErec( bot.slopeY );
  double p_top = sqrt( pow(Etop,2) - pow(constants::mass_electron,2) );
  double p_bot = sqrt( pow(Ebot,2) - pow(constants::mass_electron,2) );

  TLorentzVector electron_p;
  TLorentzVector positron_p;
  
  electron_p.SetXYZM( 
      p_top*sin(top.theta)*cos(top.phi), 
      p_top*sin(top.theta)*sin(top.phi), 
      p_top*cos(top.theta), 
      constants::mass_electron );
  // subtract py induced by B field
  double pYZ_top = sqrt( pow(electron_p.Y(), 2) + pow(electron_p.Z(), 2) ); // conserved
  electron_p.SetY( electron_p.Y() - pT );
  electron_p.SetZ( sqrt( pow(pYZ_top, 2) - pow(electron_p.Y(), 2) ) );

  positron_p.SetXYZM( 
      p_bot*sin(bot.theta)*cos(bot.phi), 
      p_bot*sin(bot.theta)*sin(bot.phi), 
      p_bot*cos(bot.theta), 
      constants::mass_electron );
  // subtract py induced by B field
  double pYZ_bot = sqrt( pow(positron_p.Y(), 2) + pow(positron_p.Z(), 2) ); // conserved
  positron_p.SetY( positron_p.Y() + pT );
  positron_p.SetZ( sqrt( pow(pYZ_bot, 2) - pow(positron_p.Y(), 2) ) );
  
  //electron_p.Print();
  //positron_p.Print();

  return (electron_p + positron_p).M();
}
