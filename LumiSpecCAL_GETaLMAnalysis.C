#include "LumiSpecCAL.h"

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TLorentzVector.h"
#include "TDirectory.h"
#include "TStyle.h"

#include "DD4hep/IDDescriptor.h"
#include "DD4hep/BitFieldCoder.h"
#include <DD4hep/Detector.h>

void LumiSpecCAL_GETaLMAnalysis(
    const TString& inputFile = "/gpfs/mnt/gpfs02/eic/aranya/data_UserDef/recFiles/test_dataset",
    const TString& outName   = "LumiSpecCAL_Histos.root",
    bool use_campaign        = false)
{
  //gStyle->SetOptStat(0);

  TChain* tree = new TChain("events");

  if (use_campaign) {
    std::ifstream in(inputFile.Data());
    if (!in.is_open()) {
      std::cerr << "Error: cannot open campaign file list: " << inputFile << std::endl;
      delete tree;
      return;
    }

    std::string file;
    while (in >> file) tree->Add(file.c_str());
    in.close();
  }
  else {
    tree->Add(inputFile.Data());
  }

  const Long64_t nEntries = tree->GetEntries();
  std::cout << "Running input: " << inputFile << std::endl;
  std::cout << "Analyzing " << nEntries << " events" << std::endl;

  if (nEntries <= 0) {
    std::cerr << "Error: chain has 0 entries. Check input." << std::endl;
    delete tree;
    return;
  }

  TFile* fout = TFile::Open(outName, "RECREATE");
  if (!fout || fout->IsZombie()) {
    std::cerr << "Error: cannot create output file: " << outName << std::endl;
    delete tree;
    return;
  }

  // --------------------------------------------------
  // Output histogram
  // --------------------------------------------------
  InitLumiSpecCALHistograms();

  // --------------------------------------------------
  // Input Tree Reader
  // --------------------------------------------------
  TTreeReader tr(tree);

  TTreeReaderArray<Int_t>    gen_status(tr, "MCParticles.generatorStatus");
  TTreeReaderArray<Int_t>    gen_pid   (tr, "MCParticles.PDG");
  TTreeReaderArray<Double_t> gen_px    (tr, "MCParticles.momentum.x");
  TTreeReaderArray<Double_t> gen_py    (tr, "MCParticles.momentum.y");
  TTreeReaderArray<Double_t> gen_pz    (tr, "MCParticles.momentum.z");
  TTreeReaderArray<Double_t> gen_mass  (tr, "MCParticles.mass");
  TTreeReaderArray<Float_t>  gen_charge(tr, "MCParticles.charge");
  TTreeReaderArray<Double_t> gen_vx    (tr, "MCParticles.vertex.x");
  TTreeReaderArray<Double_t> gen_vy    (tr, "MCParticles.vertex.y");
  TTreeReaderArray<Double_t> gen_vz    (tr, "MCParticles.vertex.z");

  TTreeReaderArray<Double_t> gen_endx  (tr, "MCParticles.endpoint.x");
  TTreeReaderArray<Double_t> gen_endy  (tr, "MCParticles.endpoint.y");
  TTreeReaderArray<Double_t> gen_endz  (tr, "MCParticles.endpoint.z");
  TTreeReaderArray<Double_t> gen_endpx (tr, "MCParticles.momentumAtEndpoint.x");
  TTreeReaderArray<Double_t> gen_endpy (tr, "MCParticles.momentumAtEndpoint.y");
  TTreeReaderArray<Double_t> gen_endpz (tr, "MCParticles.momentumAtEndpoint.z");

  //CALs hits from GEANT
  TTreeReaderArray<Float_t> cal_hits_E        (tr, "EcalLumiSpecHits.energy");
  TTreeReaderArray<ULong_t> cal_hits_cellID   (tr, "EcalLumiSpecHits.cellID");

  //CALs EICrecon reconstructed Hits
  TTreeReaderArray<Float_t> cal_rec_E         (tr, "EcalLumiSpecRecHits.energy");
  TTreeReaderArray<ULong_t> cal_rec_cellID    (tr, "EcalLumiSpecRecHits.cellID");

  //Trackers hits from GEANT
  TTreeReaderArray<Float_t> trk_hits_eDep     (tr, "LumiSpecTrackerHits.eDep");
  TTreeReaderArray<Float_t> trk_hits_time     (tr, "LumiSpecTrackerHits.time");
  TTreeReaderArray<Double_t> trk_hits_pos_x   (tr, "LumiSpecTrackerHits.position.x");
  TTreeReaderArray<Double_t> trk_hits_pos_y   (tr, "LumiSpecTrackerHits.position.y");
  TTreeReaderArray<Double_t> trk_hits_pos_z   (tr, "LumiSpecTrackerHits.position.z");
  TTreeReaderArray<Int_t> trk_hits_quality    (tr, "LumiSpecTrackerHits.quality");
  TTreeReaderArray<ULong_t> trk_hits_cellID   (tr, "LumiSpecTrackerHits.cellID");
  
 //calorimeter ID's ------------------------------------------------
  dd4hep::IDDescriptor cal_idspec(
    "EcalLumiSpecRecHits",
    "system:8,sector:8,layer:8,module:8,block:8,fiber_x:4,fiber_y:4"
  );

  dd4hep::BitFieldCoder cal_decoder(cal_idspec.fieldDescription());

  //tracker ID's ----------------------------------------------------
  dd4hep::IDDescriptor trk_idspec(
    "LumiSpecTrackerHits",
    "system:8,sector:8,module:8,x:32:-16,y:-16"
  );

  dd4hep::BitFieldCoder trk_decoder(trk_idspec.fieldDescription());

  TLorentzVector gen_moment;
  Long64_t counter = 0;
  
  // --------------------------------------------------
  // Event loop
  // --------------------------------------------------
  while (tr.Next()) { 

    if (counter % 100 == 0) {
      std::cout << "Analyzing event " << counter << std::endl;
    }

    event_id = counter;

    //GEN
    gen_pdg_vec.clear();
    gen_status_vec.clear();
    gen_charge_vec.clear();
    gen_px_vec.clear();
    gen_py_vec.clear();
    gen_pz_vec.clear();
    gen_mass_vec.clear();
    gen_energy_vec.clear();
    gen_vx_vec.clear();
    gen_vy_vec.clear();
    gen_vz_vec.clear();

    total_gen_energy_all_status1 = 0.0;
    total_gen_energy_all_status3 = 0.0;

    evt_genE_photon_status3   = 0.0;
    evt_genE_electron_status3 = 0.0;
    evt_genE_positron_status3 = 0.0;

    evt_genE_photon_status1   = 0.0;
    evt_genE_electron_status1 = 0.0;
    evt_genE_positron_status1 = 0.0;

    evt_genE_photon_status0   = 0.0;
    evt_genE_electron_status0 = 0.0;
    evt_genE_positron_status0 = 0.0;

    //CALs
    cal_rec_energy_vec.clear();
    cal_rec_cellID_vec.clear();
    cal_rec_system_vec.clear();
    cal_rec_sector_vec.clear();
    cal_rec_layer_vec.clear();
    cal_rec_module_vec.clear();
    cal_rec_block_vec.clear();
    cal_rec_fiberx_vec.clear();
    cal_rec_fibery_vec.clear();

    total_hit_energy             = 0.0;
    total_rec_energy             = 0.0;
    total_rec_energy_top         = 0.0;
    total_rec_energy_bot         = 0.0;

    trk_hits_cellID_vec.clear();
    trk_hits_system_vec.clear();
    trk_hits_sector_vec.clear();
    trk_hits_module_vec.clear();
    trk_hits_x_id_vec.clear();
    trk_hits_y_id_vec.clear();

    trk_hits_glob_x_vec.clear();
    trk_hits_glob_y_vec.clear();
    trk_hits_glob_z_vec.clear();

    trk_hits_energy_vec.clear();
    trk_hits_quality_vec.clear();
    trk_hits_time_vec.clear();

    // --------------------------------------------------
    // MC particles
    // --------------------------------------------------
    const Int_t gen_mult = gen_status.GetSize();

    for (Int_t igen = 0; igen < gen_mult; ++igen) {

      gen_moment.SetXYZM(gen_px[igen], gen_py[igen], gen_pz[igen], gen_mass[igen]);

      const Int_t status = gen_status[igen];
      const Int_t pdg    = gen_pid[igen];
      const double E     = gen_moment.E();

      gen_pdg_vec.push_back(pdg);
      gen_status_vec.push_back(status);
      gen_charge_vec.push_back(gen_charge[igen]);
      gen_px_vec.push_back(static_cast<float>(gen_px[igen]));
      gen_py_vec.push_back(static_cast<float>(gen_py[igen]));
      gen_pz_vec.push_back(static_cast<float>(gen_pz[igen]));
      gen_mass_vec.push_back(static_cast<float>(gen_mass[igen]));
      gen_energy_vec.push_back(static_cast<float>(E));
      gen_vx_vec.push_back(static_cast<float>(gen_vx[igen]));
      gen_vy_vec.push_back(static_cast<float>(gen_vy[igen]));
      gen_vz_vec.push_back(static_cast<float>(gen_vz[igen]));

      if (status == 3) {
        total_gen_energy_all_status3 += E;

        hMCGenStatus3_Vx->Fill(gen_vx[igen]);
        hMCGenStatus3_Vy->Fill(gen_vy[igen]);
        hMCGenStatus3_Vz->Fill(gen_vz[igen]);
        hMCGenStatus3_Px->Fill(gen_px[igen]);
        hMCGenStatus3_Py->Fill(gen_py[igen]);
        hMCGenStatus3_Pz->Fill(gen_pz[igen]);
        hMCGenStatus3_E->Fill(E);

        if (pdg == 22)  evt_genE_photon_status3   += E;
        if (pdg == 11)  evt_genE_electron_status3 += E;
        if (pdg == -11) evt_genE_positron_status3 += E;
      }
      else if (status == 1) {
        total_gen_energy_all_status1 += E;

        hMCGenStatus1_Vx->Fill(gen_vx[igen]);
        hMCGenStatus1_Vy->Fill(gen_vy[igen]);
        hMCGenStatus1_Vz->Fill(gen_vz[igen]);
        hMCGenStatus1_Px->Fill(gen_px[igen]);
        hMCGenStatus1_Py->Fill(gen_py[igen]);
        hMCGenStatus1_Pz->Fill(gen_pz[igen]);
        hMCGenStatus1_E->Fill(E);

        if (pdg == 22)  evt_genE_photon_status1   += E;
        if (pdg == 11)  evt_genE_electron_status1 += E;
        if (pdg == -11) evt_genE_positron_status1 += E;
      }
      else if (status == 0) {
        hMCGenStatus0_Vx->Fill(gen_vx[igen]);
        hMCGenStatus0_Vy->Fill(gen_vy[igen]);
        hMCGenStatus0_Vz->Fill(gen_vz[igen]);
        hMCGenStatus0_Px->Fill(gen_px[igen]);
        hMCGenStatus0_Py->Fill(gen_py[igen]);
        hMCGenStatus0_Pz->Fill(gen_pz[igen]);
        hMCGenStatus0_E->Fill(E);

        if (pdg == 22)  evt_genE_photon_status0   += E;
        if (pdg == 11)  evt_genE_electron_status0 += E;
        if (pdg == -11) evt_genE_positron_status0 += E;
      }
    }

    hEvtPhotonE_Status3->Fill(evt_genE_photon_status3);
    hEvtElectronE_Status3->Fill(evt_genE_electron_status3);
    hEvtPositronE_Status3->Fill(evt_genE_positron_status3);

    hEvtPhotonE_Status1->Fill(evt_genE_photon_status1);
    hEvtElectronE_Status1->Fill(evt_genE_electron_status1);
    hEvtPositronE_Status1->Fill(evt_genE_positron_status1);

    hEvtPhotonE_Status0->Fill(evt_genE_photon_status0);
    hEvtElectronE_Status0->Fill(evt_genE_electron_status0);
    hEvtPositronE_Status0->Fill(evt_genE_positron_status0);

    // --------------------------------------------------
    // CAL Rec Hits
    // --------------------------------------------------
    const Int_t cal_rec_mult = cal_rec_E.GetSize();

    for (Int_t irec = 0; irec < cal_rec_mult; ++irec) {

      const uint64_t rec_id = static_cast<uint64_t>(cal_rec_cellID[irec]);

      const int rec_system  = cal_decoder.get(rec_id, "system");
      const int rec_sector  = cal_decoder.get(rec_id, "sector");
      const int rec_layer   = cal_decoder.get(rec_id, "layer");
      const int rec_module  = cal_decoder.get(rec_id, "module");
      const int rec_block   = cal_decoder.get(rec_id, "block");
      const int rec_fiber_x = cal_decoder.get(rec_id, "fiber_x");
      const int rec_fiber_y = cal_decoder.get(rec_id, "fiber_y");

      cal_rec_energy_vec.push_back(cal_rec_E[irec]);
      cal_rec_cellID_vec.push_back(static_cast<unsigned long long>(rec_id));
      cal_rec_system_vec.push_back(rec_system);
      cal_rec_sector_vec.push_back(rec_sector);
      cal_rec_layer_vec.push_back(rec_layer);
      cal_rec_module_vec.push_back(rec_module);
      cal_rec_block_vec.push_back(rec_block);
      cal_rec_fiberx_vec.push_back(rec_fiber_x);
      cal_rec_fibery_vec.push_back(rec_fiber_y);

      hRecHitSystem->Fill(rec_system);
      hRecHitSector->Fill(rec_sector);
      hRecHitLayer->Fill(rec_layer);
      hRecHitModule->Fill(rec_module);
      hRecHitBlock->Fill(rec_block);
      hRecHitFiberX->Fill(rec_fiber_x);
      hRecHitFiberY->Fill(rec_fiber_y);

      if (rec_sector == top) {total_rec_energy_top += cal_rec_E[irec];}
      else if (rec_sector == bot) {total_rec_energy_bot += cal_rec_E[irec];}
      total_rec_energy += cal_rec_E[irec];

    }

    hEvtHitE->Fill(total_hit_energy);
    hEvtRecE_top->Fill(total_rec_energy_top);
    hEvtRecE_bot->Fill(total_rec_energy_bot);
    hEvtRecE_tot->Fill(total_rec_energy);

    // --------------------------------------------------
    // Event-level acceptance decision
    // --------------------------------------------------
    const bool acceptedPSCAL  = (total_rec_energy_top > 0.1 && total_rec_energy_bot > 0.1);
    const bool acceptedTop    = (total_rec_energy_top > 0.1);
    const bool acceptedBottom = (total_rec_energy_bot > 0.1);

    // --------------------------------------------------
    // CAL Sampling fraction and acceptance with only CALs
    // --------------------------------------------------
    if (evt_genE_photon_status3 > 0.0) {

      const double sf_hits = total_hit_energy / evt_genE_photon_status3;
      const double sf_rec  = total_rec_energy / evt_genE_photon_status3;

      if (total_hit_energy > 0.2) hSamplingFraction_hits->Fill(sf_hits);
      if (total_rec_energy > 0.2) hSamplingFraction_rec->Fill(sf_rec);

      hGenPSCAL->Fill(evt_genE_photon_status3);
      hGenUpCAL->Fill(evt_genE_photon_status3);
      hGenDwCAL->Fill(evt_genE_photon_status3);

      if (acceptedPSCAL)  hAccPSCAL->Fill(evt_genE_photon_status3);
      if (acceptedTop)    hAccUpCAL->Fill(evt_genE_photon_status3);
      if (acceptedBottom) hAccDwCAL->Fill(evt_genE_photon_status3);
    }

    // --------------------------------------------------
    // ------------   TRACKER Hits   --------------------
    // --------------------------------------------------

    //There are two modules of trackers (closer to CAL-0, further from CAL-1)
    std::vector <TVector3> trk_pos[NSector][NTracker]; 
    const Int_t trk_hits_mult = trk_hits_cellID.GetSize();

    for (Int_t irec = 0; irec < trk_hits_mult; ++irec) {

      bool primary = trk_hits_quality[irec] ==0 ;
      if(!primary) continue; //only take into account of primary tracks for the moment.

      const uint64_t rec_id = static_cast<uint64_t>(trk_hits_cellID[irec]);

      // tracker id's
      const int rec_system  = trk_decoder.get(rec_id, "system");
      const int rec_sector  = trk_decoder.get(rec_id, "sector");
      const int rec_module  = trk_decoder.get(rec_id, "module");
      const int rec_x_id    = trk_decoder.get(rec_id, "x");
      const int rec_y_id    = trk_decoder.get(rec_id, "y");

      trk_hits_cellID_vec.push_back( static_cast<unsigned long long>(rec_id));
      trk_hits_system_vec.push_back( rec_system);
      trk_hits_sector_vec.push_back( rec_sector);
      trk_hits_module_vec.push_back( rec_module);
      trk_hits_x_id_vec.push_back(   rec_x_id);
      trk_hits_y_id_vec.push_back(   rec_y_id);

      hTrkRecHitSystem->Fill( rec_system);
      hTrkRecHitSector->Fill( rec_sector);
      hTrkRecHitModule->Fill( rec_module);
      hTrkRecHitXId->Fill(    rec_x_id);
      hTrkRecHitYId->Fill(    rec_y_id);

      hTrkRecHitGlobalX->Fill(  trk_hits_pos_x[irec] );
      hTrkRecHitGlobalY->Fill(  trk_hits_pos_y[irec] );
      hTrkRecHitGlobalZ->Fill(  trk_hits_pos_z[irec] );
      hTrkRecHitGlobalXY->Fill( trk_hits_pos_x[irec], trk_hits_pos_y[irec]);
      hTrkRecHitGlobalXZ->Fill( trk_hits_pos_x[irec], trk_hits_pos_z[irec]);

      trk_hits_energy_vec.push_back(   trk_hits_eDep[irec]);
      //trk_hits_quality_vec.push_back(  trk_hits_quality[irec]);
      trk_hits_time_vec.push_back(     trk_hits_time[irec]);

      hTrkRecHitEnergy->Fill(   trk_hits_eDep[irec]);
      hTrkRecHitTime->Fill(     trk_hits_time[irec]);
      hTrkRecHitPrimary->Fill(  primary);

      //primary track hits
      if (primary && rec_sector >= 0 && rec_sector < NSector && rec_module >= 0 && rec_module < NTracker) {
        trk_pos[rec_sector][rec_module].emplace_back( trk_hits_pos_x[irec], trk_hits_pos_y[irec], trk_hits_pos_z[irec] );
      }

    }// --------- reconsted trk hits closed -------------

    // --------------------------------------------------
    // --------  MC Photon endpoint loss maps  ----------
    // --------------------------------------------------

    for (Int_t igen = 0; igen < gen_mult; ++igen) {

      const Int_t status = gen_status[igen];
      const Int_t pdg    = gen_pid[igen];
      const double ez  = gen_endz[igen];
      const double ex  = gen_endx[igen];
      const double ey  = gen_endy[igen];
      const double epx = gen_endpx[igen];
      const double epy = gen_endpy[igen];
      const double epz = gen_endpz[igen];
      const double gvx = gen_vx[igen];
      const double gvy = gen_vy[igen];
      const double gvz = gen_vz[igen];

      gen_moment.SetXYZM(gen_px[igen], gen_py[igen], gen_pz[igen], gen_mass[igen]);
      const double gE  = gen_moment.E();

      //+++++++++++++++ for standalone simualtion only +++++++++++++++++++
      // Acceptance File has : photon (3) , e+ and e- as (1)
      // Resolution File has : e+ /e- as (1)
      if (status == 1 && std::abs(pdg) == 11) {

        TVector3 centroidPos[NSector] = {}; //even though only one sec is analyzed
        double centroidEnergy[NSector] = {}; //it receives values from a generalised function.

        TVector3 slopePos[NSector] = {};
        double slopeEnergy[NSector] = {};

        double recVx  = -999.99;
        double recVy  = -999.99;
        double diffx  = -999.99;
        double diffy  = -999.99;
        double recE   = -999.99;

        hLeptonEndX->Fill(ex);
        hLeptonEndY->Fill(ey);
        hLeptonEndZ->Fill(ez);
        hLeptonPDG_vs_EndY->Fill(pdg, ey);

        int sec = NSector;
        if( (ey > 67.0) && (ez < -63990) ){
          sec =top;
        }
        else if((ey < -67.0) && (ez < -63990)){
          sec =bot;
        }

        if (trk_pos[sec][0].size() == 1 && sec <NSector) {

            const double trkHitx = trk_pos[sec][0].front().X() / 10.0;
            const double trkHity = trk_pos[sec][0].front().Y() / 10.0;

            // -------- 1. Centroid Method  ---------
            double w0 = 4.0;
            AccumulateWeightedCentroid(&cal_rec_energy_vec, &cal_rec_sector_vec, &cal_rec_layer_vec,
                        &cal_rec_module_vec, &cal_rec_block_vec, w0, centroidPos, centroidEnergy );

            if(centroidEnergy[sec] >0.0){

              recVx = centroidPos[sec].X();
              recVy = centroidPos[sec].Y();

              diffx = recVx - trkHitx;
              diffy = recVy - trkHity;

              recE = centroidEnergy[sec];

              //==================
              // 1D histograms
              //==================
              hDiffElectronX_Ctd[sec]->Fill(diffx);
              hDiffElectronY_Ctd[sec]->Fill(diffy);

              hRecElectronVx_Ctd[sec]->Fill(recVx);
              hRecElectronVy_Ctd[sec]->Fill(recVy);

              hGenElectronVx_Ctd[sec]->Fill(trkHitx);
              hGenElectronVy_Ctd[sec]->Fill(trkHity);
              //==================
              // 2D histograms
              //==================
              hDiffElectronX_vs_gE_Ctd[sec]->Fill(gE, diffx);
              hDiffElectronY_vs_gE_Ctd[sec]->Fill(gE, diffy);

              hRecElectronVx_vs_gE_Ctd[sec]->Fill(gE, recVx);
              hRecElectronVy_vs_gE_Ctd[sec]->Fill(gE, recVy);

              hTrkHitx_vs_gE_Ctd[sec]->Fill(gE, trkHitx);
              hTrkHity_vs_gE_Ctd[sec]->Fill(gE, trkHity);

              hRecEnergy_vs_gE_Ctd[sec]->Fill(gE, recE);

            }

            // ---------- 2. Slope Method ---------
            AccumulateWeightedSlope(&cal_rec_energy_vec, &cal_rec_sector_vec, &cal_rec_layer_vec,
                                    &cal_rec_module_vec, &cal_rec_block_vec, trk_pos, slopePos, slopeEnergy);

            if(slopeEnergy[sec] >0.0){

              recVx = slopePos[sec].X();
              recVy = slopePos[sec].Y();

              diffx = recVx - trkHitx;
              diffy = recVy - trkHity;

              recE = slopeEnergy[sec];
              //==================
              // 1D histograms
              //==================
              hDiffElectronX_Slope[sec]->Fill(diffx);
              hDiffElectronY_Slope[sec]->Fill(diffy);

              hRecElectronVx_Slope[sec]->Fill(recVx);
              hRecElectronVy_Slope[sec]->Fill(recVy);

              hGenElectronVx_Slope[sec]->Fill(trkHitx);
              hGenElectronVy_Slope[sec]->Fill(trkHity);

              //==================
              // 2D histograms
              //==================
              hDiffElectronX_vs_gE_Slope[sec]->Fill(gE, diffx);
              hDiffElectronY_vs_gE_Slope[sec]->Fill(gE, diffy);

              hRecElectronVx_vs_gE_Slope[sec]->Fill(gE, recVx);
              hRecElectronVy_vs_gE_Slope[sec]->Fill(gE, recVy);

              hTrkHitx_vs_gE_Slope[sec]->Fill(gE, trkHitx);
              hTrkHity_vs_gE_Slope[sec]->Fill(gE, trkHity);
              hRecEnergy_vs_gE_Slope[sec]->Fill(gE, recE);
            }
        }
    }//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      //++++++++++++++++ for complete simulation ++++++++++++++++++++++++++
      if (status != 3 || pdg != 22) continue; // only photons shall pass this criteria.
      if (ez < -70000.0 || ez > 0.0) continue; // which are lost in far backward region.

      hPhotonEndX_vs_Z_all->Fill(ez, ex);
      hPhotonEndY_vs_Z_all->Fill(ez, ey);
      hPhotonEndPx_vs_Z_all->Fill(ez, epx);
      hPhotonEndPy_vs_Z_all->Fill(ez, epy);
      hPhotonEndPz_vs_Z_all->Fill(ez, epz);

      if (acceptedPSCAL) {
        hPhotonAcceptedZ->Fill(ez);

        hPhotonEndX_vs_Z_acc->Fill(ez, ex);
        hPhotonEndY_vs_Z_acc->Fill(ez, ey);
        hPhotonEndPx_vs_Z_acc->Fill(ez, epx);
        hPhotonEndPy_vs_Z_acc->Fill(ez, epy);
        hPhotonEndPz_vs_Z_acc->Fill(ez, epz);

      }
      else {
        hPhotonLostZ->Fill(ez);

        hPhotonEndX_vs_Z_lost->Fill(ez, ex);
        hPhotonEndY_vs_Z_lost->Fill(ez, ey);
        hPhotonEndPx_vs_Z_lost->Fill(ez, epx);
        hPhotonEndPy_vs_Z_lost->Fill(ez, epy);
        hPhotonEndPz_vs_Z_lost->Fill(ez, epz);
      }

      //detector wise
      for (const auto& det : detectorRanges) {

        if (ez < det.zmin || ez > det.zmax) continue;

        //pair - accepted
        if (acceptedPSCAL && det.saveAccepted) {
          auto& h = hAccDet[det.key]; 

          h.hEndY_vs_EndX->Fill(ex, ey); 
          h.hGenY_vs_GenX->Fill(gvx, gvy); 
          h.hGenZ->Fill(gvz); 
          h.hGenE->Fill(gE); 

          //---------- reconstruct position only for thin convereter ------
          if (det.key == "ThinFilm"){ 
        
            double  w0    = 4.0; 
            TVector3 centroidPos[NSector] = {};
            double centroidEnergy[NSector] = {};
            double vxRec = 0.0;
            double vyRec = 0.0;

            //centroid method using CALs ------
            AccumulateWeightedCentroid( &cal_rec_energy_vec, &cal_rec_sector_vec, &cal_rec_layer_vec,
                                        &cal_rec_module_vec, &cal_rec_block_vec, w0, centroidPos, centroidEnergy );
                                      
            if ( (centroidEnergy[top] >0.0) && (centroidEnergy[bot] >0.0) ) {
              //this mode ensures the denomenator is not zero.
              vxRec = ( centroidPos[top].X() + centroidPos[bot].X() ) / 2.0 ;
              vyRec = ( centroidPos[top].Y()*centroidEnergy[top] + centroidPos[bot].Y()*centroidEnergy[bot] )/ ( centroidEnergy[top] + centroidEnergy[bot] );

              //reconstructed - lost photon coordinates (cm = mm/10)
              hDiffX_Ctd->Fill( vxRec - (ex/10.0)  );
              hDiffY_Ctd->Fill( vyRec - (ey/10.0) ); 

              hRecPhotonAtThinFilmVx_Ctd->Fill( vxRec ); 
              hRecPhotonAtThinFilmVy_Ctd->Fill( vyRec ); 

              hGenPhotonAtThinFilmVx_Ctd->Fill( (ex/10.0) ); 
              hGenPhotonAtThinFilmVy_Ctd->Fill( (ey/10.0) );
            }

            TVector3 slopePos[NSector] = {};
            double slopeEnergy[NSector] = {};
            vxRec = 0.0;
            vyRec = 0.0;

            //Slope Method using CALs ------
            AccumulateWeightedSlope(&cal_rec_energy_vec, &cal_rec_sector_vec, &cal_rec_layer_vec,
                                    &cal_rec_module_vec, &cal_rec_block_vec, trk_pos, slopePos, slopeEnergy);
                                      
            if ( (slopeEnergy[top] >0.0) && (slopeEnergy[bot] >0.0) ) {
              //this mode ensures the denomenator is not zero.
              vxRec = ( slopePos[top].X() + slopePos[bot].X() ) / 2.0 ;
              vyRec = ( slopePos[top].Y()*slopeEnergy[top] + slopePos[bot].Y()*slopeEnergy[bot] )/ ( slopeEnergy[top] + slopeEnergy[bot] );

              hDiffX_Slope->Fill( vxRec - (ex/10.0)); 
              hDiffY_Slope->Fill( vyRec - (ey/10.0)); 

              hRecPhotonAtThinFilmVx_Slope->Fill( vxRec); 
              hRecPhotonAtThinFilmVy_Slope->Fill( vyRec); 

              hGenPhotonAtThinFilmVx_Slope->Fill( (ex/10.0) ); 
              hGenPhotonAtThinFilmVy_Slope->Fill( (ey/10.0) );
            }

            // ------  Based on trackers for x-check -----
            if (  (trk_pos[top][0].size() ==1) && (trk_pos[bot][0].size() ==1) && 
                  (slopeEnergy[top] >0.0) && (slopeEnergy[bot] >0.0) ) {
              
              //this mode ensures the denomenator is not zero.
              vxRec = ( trk_pos[top][0].front().X() + trk_pos[bot][0].front().X() ) / (2*10.0) ; //convert it into cm
              vyRec = ( trk_pos[top][0].front().Y()*slopeEnergy[top] + trk_pos[bot][0].front().Y()*slopeEnergy[bot] )/ ( (slopeEnergy[top] + slopeEnergy[bot])*10.0 ); //convert it into cm

              hDiffX_Trk->Fill( vxRec - (ex/10.0)); 
              hDiffY_Trk->Fill( vyRec - (ey/10.0)); 

              hRecPhotonAtThinFilmVx_Trk->Fill( vxRec); 
              hRecPhotonAtThinFilmVy_Trk->Fill( vyRec); 

              hGenPhotonAtThinFilmVx_Trk->Fill( (ex/10.0) ); 
              hGenPhotonAtThinFilmVy_Trk->Fill( (ey/10.0) );
            }
          } //--- thin film close 
          //------------------------------------------------------------------------
        }

        //pair - NOT accepted
        if (!acceptedPSCAL && det.saveLost) {
          auto& h = hLostDet[det.key];

          h.hEndY_vs_EndX->Fill(ex, ey);
          h.hGenY_vs_GenX->Fill(gvx, gvy);
          h.hGenZ->Fill(gvz);
          h.hGenE->Fill(gE);
        }
      }
    } // gen particle loop
    counter++;
  }

  // Build acceptance histograms
  hAcceptancePSCAL = (TH1D*)hAccPSCAL->Clone("hAcceptancePSCAL");
  hAcceptanceUpCAL = (TH1D*)hAccUpCAL->Clone("hAcceptanceUpCAL");
  hAcceptanceDwCAL = (TH1D*)hAccDwCAL->Clone("hAcceptanceDwCAL");

  hAcceptancePSCAL->SetTitle("PSCAL Acceptance;E_{#gamma}^{gen} [GeV];Acceptance");
  hAcceptanceUpCAL->SetTitle("Top CAL Acceptance;E_{#gamma}^{gen} [GeV];Acceptance");
  hAcceptanceDwCAL->SetTitle("Bottom CAL Acceptance;E_{#gamma}^{gen} [GeV];Acceptance");

  hAcceptancePSCAL->Divide(hAccPSCAL, hGenPSCAL, 1.0, 1.0, "B");
  hAcceptanceUpCAL->Divide(hAccUpCAL, hGenUpCAL, 1.0, 1.0, "B");
  hAcceptanceDwCAL->Divide(hAccDwCAL, hGenDwCAL, 1.0, 1.0, "B");

  // --------------------------------------------------
  // Write output
  // --------------------------------------------------
  SaveInFile(fout);

  delete tree;

  std::cout << "Saved histograms to: " << outName << std::endl;
  std::cout << "Main loss plots are in PhotonLossMap/" << std::endl;
  std::cout << "Check hPhotonEndX_vs_Z_lost, hPhotonEndY_vs_Z_lost, and hPhotonLostZ." << std::endl;
}

//_________________________________________________________________________________________________________
PhotonDetectorHistSet MakePhotonDetectorHistSet(
  const TString& prefix,
  const DetectorZRange& det,
  const TString& sample
)
{
  PhotonDetectorHistSet h;

  h.hEndY_vs_EndX = new TH2D(
    Form("%s_%s_hEndY_vs_EndX", prefix.Data(), det.key.Data()),
    Form("%s photons in %s;endpoint x [mm];endpoint y [mm]",
         sample.Data(), det.label.Data()),
    600, -300.0, 300.0,
    600, -300.0, 300.0
  );

  h.hGenY_vs_GenX = new TH2D(
    Form("%s_%s_hGenY_vs_GenX", prefix.Data(), det.key.Data()),
    Form("%s photons in %s;gen x [mm];gen y [mm]",
         sample.Data(), det.label.Data()),
    600, -300.0, 300.0,
    600, -300.0, 300.0
  );

  h.hGenZ = new TH1D(
    Form("%s_%s_hGenZ", prefix.Data(), det.key.Data()),
    Form("%s photons in %s;gen z [mm];Photons",
         sample.Data(), det.label.Data()),
    71000, -70000.0, 1000.0
  );

  h.hGenE = new TH1D(
    Form("%s_%s_hGenE", prefix.Data(), det.key.Data()),
    Form("%s photons in %s;E_{#gamma}^{gen} [GeV];Photons",
         sample.Data(), det.label.Data()),
    500, 0.0, 50.0
  );

  return h;
}

//_________________________________________________________________________________________________________
void AccumulateWeightedSlope(const std::vector<float>*  rec_energy,
                                       const std::vector<int>* rec_sector,
                                       const std::vector<int>* rec_layer,
                                       const std::vector<int>* rec_module,
                                       const std::vector<int>* rec_block,
                                       const std::vector<TVector3> tracker_pos[NSector][NTracker],
                                       TVector3 posRec[NSector],
                                       double energyRec[NSector])
{
  double vxRec = 0.0 ;
  double vyRec = 0.0 ;

  //SECTOR WISE
  double vx[NSector] = {0.0};
  double vy[NSector] = {0.0};

  if (!rec_energy || !rec_sector || !rec_layer || !rec_module || !rec_block) return;

  const size_t nHits = rec_energy->size();
  if (rec_sector->size() != nHits || rec_layer->size() != nHits ||
      rec_module->size() != nHits || rec_block->size() != nHits) {
    return;
  }

  //This is for 4x4 mm SiPMs
  int countx[NSector][NLayer][NModule][NBlockX][NBlockY]   = {};
  int county[NSector][NLayer][NModule][NBlockX][NBlockY]   = {};

  double x_energy[NSector][NLayer][NModule][NBlockX][NBlockY] = {};
  double y_energy[NSector][NLayer][NModule][NBlockX][NBlockY] = {};

  double x_pos[NSector][NLayer][NModule][NBlockX][NBlockY]    = {};
  double y_pos[NSector][NLayer][NModule][NBlockX][NBlockY]    = {};
  double z_pos[NSector][NLayer][NModule][NBlockX][NBlockY]    = {};

  double total_energy[NSector] = {};

  // rec hit loop
  for (size_t ih = 0; ih < nHits; ++ih) {
    
    const int sec   = rec_sector->at(ih);
    const int lay   = rec_layer->at(ih);
    const int mod   = rec_module->at(ih);
    const int blk   = rec_block->at(ih);
    const double e  = rec_energy->at(ih);

    if (sec < 0 || sec >= NSector)  continue;
    if (lay < 0 || lay >= NLayer)   continue;
    if (mod < 0 || mod >= NModule)  continue;
    if (e <= 0.0) continue;

    const int blockx = blk % NBlockX;
    const int blocky = blk / NBlockX;
    if (blockx < 0 || blockx >= NBlockX) continue;
    if (blocky < 0 || blocky >= NBlockY) continue;

    if (lay % 2 == 0) {

      const double y = cal_y_low_pos[sec] + mod * mod_x_size + (lay_coat_size + mod_coat_size) + (blockx * block_xy_size + block_xy_size / 2.0);
      y_pos[sec][lay][mod][blockx][blocky]      = y;
      county[sec][lay][mod][blockx][blocky]    += 1.0;
      y_energy[sec][lay][mod][blockx][blocky]  += e;

    } else {

      const double x = cal_x_low_pos + mod * mod_x_size + (lay_coat_size + mod_coat_size) + (blockx * block_xy_size + block_xy_size / 2.0);
      x_pos[sec][lay][mod][blockx][blocky]      = x;
      countx[sec][lay][mod][blockx][blocky]    += 1.0;
      x_energy[sec][lay][mod][blockx][blocky]  += e;

    } //layer orientation if-else

    const double z = cal_z_front_pos + lay*(lay_y_size) + (lay_coat_size + mod_coat_size) + (blocky*block_xy_size + block_xy_size/2.0) ;
    z_pos[sec][lay][mod][blockx][blocky] = z;

    total_energy[sec] += e;

  } // rec-hits loop closed

  //calculate slope and intercept for each sector
  for(int isec = 0; isec < NSector; ++isec){

    //if electron hits only one cal
    if(total_energy[isec] <= 0) continue;

    double total_YE = 0.0;
		double total_XE = 0.0;

		double sliceZ[NSlice]={0.0};
		
		double sliceXE[NSlice]={0.0};
		double sliceYE[NSlice]={0.0};

		double sliceXDotE[NSlice]={0.0};
		double sliceYDotE[NSlice]={0.0};

    //calculating the x,y,z positions
		for(int ilay =0; ilay<NLayer; ilay++){
      for(int imod =0; imod<NModule; imod++){
				for(int iblocx =0; iblocx<NBlockX; iblocx++){
					for(int iblocy =0; iblocy<NBlockY; iblocy++){

						int slice_number = ilay*2 + iblocy;
						sliceZ[slice_number] = z_pos[isec][ilay][imod][iblocx][iblocy];

						if(ilay%2 == 0 && county[isec][ilay][imod][iblocx][iblocy] >0)//parallel to x
						{
							sliceYDotE[slice_number]  += y_pos[isec][ilay][imod][iblocx][iblocy]*y_energy[isec][ilay][imod][iblocx][iblocy];
							sliceYE[slice_number]     += y_energy[isec][ilay][imod][iblocx][iblocy];
							total_YE                  += y_energy[isec][ilay][imod][iblocx][iblocy];

						}
						if(ilay%2 != 0 && countx[isec][ilay][imod][iblocx][iblocy] >0)//parallel to y
						{
							sliceXDotE[slice_number]  += x_pos[isec][ilay][imod][iblocx][iblocy]*x_energy[isec][ilay][imod][iblocx][iblocy];
							sliceXE[slice_number]     += x_energy[isec][ilay][imod][iblocx][iblocy];
							total_XE                  += x_energy[isec][ilay][imod][iblocx][iblocy];

						}

					}//blocky
				}//blockx
      }//module
		}//layer

    //Calculation of m(slope) and c(x/y-intercept) using Least Square Fit.
		double Xrec 	= 0.0;
		double Yrec	= 0.0;
		double sum_x	= 0.0;
		double sum_y	= 0.0;

		double slope_x = 0.0;
		double slope_y = 0.0;
		double x_int = 0.0;
		double y_int = 0.0;

		double sum_zx = 0.0;	
		double sum_zy	= 0.0;
		double count_x = 0;
		double count_y = 0;

		//for x-z angle
		double E_wt = 0.0;
		double Zrec = 0.0;
		double sum_z	= 0.0;
		double sum_zz	= 0.0;

    for(int islice=0; islice<NSlice; islice++){

			if( sliceXE[islice]>0.0 ){

				Xrec 	= sliceXDotE[islice]/sliceXE[islice];
				E_wt 	= sliceXE[islice]/total_XE; 
				Zrec 	= sliceZ[islice]; //cm

				sum_x 	+= E_wt*Xrec;
				sum_z 	+= E_wt*Zrec;
				sum_zx 	+= E_wt*Xrec*Zrec;
				sum_zz 	+= E_wt*Zrec*Zrec;

				count_x	+= E_wt;
			}//>0 check	
		}//klayer-close

    if(count_x <= 0) continue;

    double denominator = count_x * sum_zz - sum_z * sum_z;
    if (std::abs(denominator) > 1.0e-12) {
      slope_x = (count_x * sum_zx - sum_z * sum_x) / denominator;
    }else{
      slope_x = 0;
    }
    
		x_int 	= (sum_x - slope_x*sum_z)/count_x; //cm

		//for y-z angle
		E_wt = 0.0;
		Zrec = 0.0;
		sum_z	= 0.0;
		sum_zz = 0.0;

		for(int islice=0; islice<NSlice; islice++){

			if( sliceYE[islice]>0.0 ){
				Yrec 	= sliceYDotE[islice]/sliceYE[islice];
				E_wt 	= sliceYE[islice]/total_YE; 
				Zrec 	= sliceZ[islice]; //cm
				
				sum_y 	+= E_wt*Yrec;
				sum_z 	+= E_wt*Zrec;
				sum_zy 	+= E_wt*Yrec*Zrec;
				sum_zz 	+= E_wt*Zrec*Zrec;

				count_y += E_wt;
			}//>0 check
		}//klayer-close

    if(count_y <= 0) continue;

		denominator = count_y * sum_zz - sum_z * sum_z;
    if (std::abs(denominator) > 1.0e-12) {
      slope_y = (count_y * sum_zy - sum_z * sum_y) / denominator;
    }else{
      slope_y = 0;
    }

		y_int	= (sum_y - slope_y*sum_z)/count_y; //cm

    //hit location in 2nd tracker(module - 0) z position.
    
    if (tracker_pos[isec][0].empty()) continue;

    vx[isec] = slope_x*(tracker_pos[isec][0].front().Z()/10 ) + x_int; //each tracker should have 1 hit ideally
    vy[isec] = slope_y*(tracker_pos[isec][0].front().Z()/10 ) + y_int;

    posRec[isec].SetXYZ(vx[isec], vy[isec], 0.0);

    energyRec[isec] = total_energy[isec];

  }//sector close

  return; 
}

//_________________________________________________________________________________________________________
void AccumulateWeightedCentroid(const std::vector<float>* rec_energy,
                                       const std::vector<int>* rec_sector,
                                       const std::vector<int>* rec_layer,
                                       const std::vector<int>* rec_module,
                                       const std::vector<int>* rec_block,
                                       double w0,
                                       TVector3 posRec[NSector],
                                       double energyRec[NSector] )
{
  //SECTOR WISE
  double vx[NSector] = {0.0};
  double vy[NSector] = {0.0};

  if (!rec_energy || !rec_sector || !rec_layer || !rec_module || !rec_block) return;

  const size_t nHits = rec_energy->size();
  if (rec_sector->size() != nHits || rec_layer->size() != nHits ||
      rec_module->size() != nHits || rec_block->size() != nHits) {
    return;
  }

  //this is for bigger 8x8 mm SIPMs
  double countx[NSector][NLayer][NModule][NBlockX / 2]    = {};
  double county[NSector][NLayer][NModule][NBlockX / 2]    = {};
  double xpos[NSector][NLayer][NModule][NBlockX / 2]      = {};
  double ypos[NSector][NLayer][NModule][NBlockX / 2]      = {};
  double x_energy[NSector][NLayer][NModule][NBlockX / 2]  = {};
  double y_energy[NSector][NLayer][NModule][NBlockX / 2]  = {};

  //This is for 4x4 mm SiPMs
  double countx2[NSector][NLayer][NModule][NBlockX][NBlockY]    = {};
  double county2[NSector][NLayer][NModule][NBlockX][NBlockY]    = {};
  double xpos2[NSector][NLayer][NModule][NBlockX][NBlockY]      = {};
  double ypos2[NSector][NLayer][NModule][NBlockX][NBlockY]      = {};
  double x_energy2[NSector][NLayer][NModule][NBlockX][NBlockY]  = {};
  double y_energy2[NSector][NLayer][NModule][NBlockX][NBlockY]  = {};

  double total_x_energy[NSector] = {};
  double total_y_energy[NSector] = {};

  for (size_t ih = 0; ih < nHits; ++ih) {
    const int sec = rec_sector->at(ih);
    const int lay = rec_layer->at(ih);
    const int mod = rec_module->at(ih);
    const int blk = rec_block->at(ih);
    const double e = rec_energy->at(ih);

    if (sec < 0 || sec >= NSector) continue;
    if (lay < 0 || lay >= NLayer) continue;
    if (mod < 0 || mod >= NModule) continue;
    if (e <= 0.0) continue;

    //use bigger SiPMs in the farther end of CAL
    if (lay > 13) {
      const int dummyx = blk % NBlockX;
      const int blockx = dummyx / 2;
      if (blockx < 0 || blockx >= (NBlockX / 2)) continue;

      if (lay % 2 == 0) {

        const double y = cal_y_low_pos[sec] + mod * mod_x_size + (lay_coat_size + mod_coat_size) + (blockx * (block_xy_size * 2.0) + block_xy_size);
        ypos[sec][lay][mod][blockx]     = y;
        county[sec][lay][mod][blockx]   += 1.0;
        y_energy[sec][lay][mod][blockx] += e;
        total_y_energy[sec] += e;

      } else {

        const double x = cal_x_low_pos + mod * mod_x_size + (lay_coat_size + mod_coat_size) + (blockx * (block_xy_size * 2.0) + block_xy_size);
        xpos[sec][lay][mod][blockx]     = x;
        countx[sec][lay][mod][blockx]   += 1.0;
        x_energy[sec][lay][mod][blockx] += e;
        total_x_energy[sec] += e;

      }
    } else {

      const int blockx = blk % NBlockX;
      const int blocky = blk / NBlockX;
      if (blockx < 0 || blockx >= NBlockX) continue;
      if (blocky < 0 || blocky >= NBlockY) continue;

      if (lay % 2 == 0) {

        const double y = cal_y_low_pos[sec] + mod * mod_x_size + (lay_coat_size + mod_coat_size) + (blockx * block_xy_size + block_xy_size / 2.0);
        ypos2[sec][lay][mod][blockx][blocky]      = y;
        county2[sec][lay][mod][blockx][blocky]    += 1.0;
        y_energy2[sec][lay][mod][blockx][blocky]  += e;
        total_y_energy[sec] += e;

      } else {

        const double x = cal_x_low_pos + mod * mod_x_size + (lay_coat_size + mod_coat_size) + (blockx * block_xy_size + block_xy_size / 2.0);
        xpos2[sec][lay][mod][blockx][blocky]      = x;
        countx2[sec][lay][mod][blockx][blocky]    += 1.0;
        x_energy2[sec][lay][mod][blockx][blocky]  += e;
        total_x_energy[sec] += e;

      } //layer orientation if-else
    } //layers front-back if-else
  } // rec-hits loop closed

  double xpos_dot_energy[NSector]   = { 0.0 };
  double xpos_total_weight[NSector] = { 0.0 };
  double ypos_dot_energy[NSector]   = { 0.0 };
  double ypos_total_weight[NSector] = { 0.0 };
  double total_energy[NSector]      = { 0.0 };

  const double minFrac = std::exp(-w0);

  for(int isec = 0; isec < NSector; ++isec){

    if (total_x_energy[isec] <= 0.0 || total_y_energy[isec] <= 0.0 ) continue;

    for (int ilay = 0; ilay < NLayer; ++ilay) {
      for(int imod =0; imod < NModule; imod++){

        if (ilay > 13) {
          for (int ibx = 0; ibx < (NBlockX / 2); ++ibx) {
            if (county[isec][ilay][imod][ibx] > 0.0) {
              const double frac = y_energy[isec][ilay][imod][ibx] / total_y_energy[isec];
              const double wi = (frac > minFrac) ? (w0 + std::log(frac)) : 0.0;
              ypos_dot_energy[isec]   += ypos[isec][ilay][imod][ibx] * wi;
              ypos_total_weight[isec] += wi;
            }
            if (countx[isec][ilay][imod][ibx] > 0.0) {
              const double frac = x_energy[isec][ilay][imod][ibx] / total_x_energy[isec];
              const double wi = (frac > minFrac) ? (w0 + std::log(frac)) : 0.0;
              xpos_dot_energy[isec]   += xpos[isec][ilay][imod][ibx] * wi;
              xpos_total_weight[isec] += wi;
            }
          }
        } else {
          for (int ibx = 0; ibx < NBlockX; ++ibx) {
            for (int iby = 0; iby < NBlockY; ++iby) {
              if (county2[isec][ilay][imod][ibx][iby] > 0.0) {
                const double frac = y_energy2[isec][ilay][imod][ibx][iby] / total_y_energy[isec];
                const double wi = (frac > minFrac) ? (w0 + std::log(frac)) : 0.0;
                ypos_dot_energy[isec]   += ypos2[isec][ilay][imod][ibx][iby] * wi;
                ypos_total_weight[isec] += wi;
              }
              if (countx2[isec][ilay][imod][ibx][iby] > 0.0) {
                const double frac = x_energy2[isec][ilay][imod][ibx][iby] / total_x_energy[isec];
                const double wi = (frac > minFrac) ? (w0 + std::log(frac)) : 0.0;
                xpos_dot_energy[isec]   += xpos2[isec][ilay][imod][ibx][iby] * wi;
                xpos_total_weight[isec] += wi;
              }
            } //block y close
          }//block x close
        } //layer >13 if-else close
      }//module close
    } //layer close

  total_energy[isec] += total_x_energy[isec] + total_y_energy[isec]; //sector energy for photon position

}//sector close

//centroid position inside the calorimeter
for (int isec = 0; isec < NSector; ++isec) {

  if (xpos_total_weight[isec] <= 0.0 || ypos_total_weight[isec] <= 0.0) continue;

  vx[isec] = xpos_dot_energy[isec] / xpos_total_weight[isec];
  vy[isec] = ypos_dot_energy[isec] / ypos_total_weight[isec];

  posRec[isec].SetXYZ(vx[isec], vy[isec], 0.0);
  energyRec[isec] = total_energy[isec];
}

  return;
}

//_________________________________________________________________________________________________________
void InitLumiSpecCALHistograms()
{
  // Event histograms
  hEvtHitE = new TH1D("hEvtHitE",
    "LumiSpecCAL Event Energy from GEANT4 Hits;E_{sum}^{hits} [GeV];Events",
    800, 0.0, 2.0);

  hEvtRecE_top = new TH1D("hEvtRecE_top",
    "LumiSpecCAL Event Energy from RecHits (top);E_{sum}^{rec} [GeV];Events",
    800, 0.0, 2.0);

  hEvtRecE_bot = new TH1D("hEvtRecE_bot",
    "LumiSpecCAL Event Energy from RecHits (bottom);E_{sum}^{rec} [GeV];Events",
    800, 0.0, 2.0);

  hEvtRecE_tot = new TH1D("hEvtRecE_tot",
    "LumiSpecCAL Event Energy from RecHits (total);E_{sum}^{rec} [GeV];Events",
    800, 0.0, 2.0);

  hEvtDiff_hits_minus_rec = new TH1D("hEvtDiff_hits_minus_rec",
    "Event Energy Difference;|E_{sum}^{hits} - E_{sum}^{rec}| [GeV];Events",
    800, 0.0, 2.0);

  hSamplingFraction_hits = new TH1D("hSamplingFraction_hits",
    "Sampling Fraction (Hits / E_{#gamma}^{gen});E_{hits}/E_{#gamma}^{gen};Events",
    500, 0.0, 1.5);

  hSamplingFraction_rec = new TH1D("hSamplingFraction_rec",
    "Sampling Fraction (RecHits / E_{#gamma}^{gen});E_{rec}/E_{#gamma}^{gen};Events",
    500, 0.0, 1.5);

  // Hit ID histograms
  hHitSystem = new TH1D("hHitSystem", "Hit system;system;Hits", 256, -0.5, 255.5);
  hHitSector = new TH1D("hHitSector", "Hit sector;sector;Hits", 2, -0.5, 1.5);
  hHitLayer  = new TH1D("hHitLayer",  "Hit layer;layer;Hits", 31, -0.5, 30.5);
  hHitModule = new TH1D("hHitModule", "Hit module;module;Hits", 3, -0.5, 2.5);
  hHitBlock  = new TH1D("hHitBlock",  "Hit block;block;Hits", 256, -0.5, 255.5);
  hHitFiberX = new TH1D("hHitFiberX", "Hit fiber x;fiber_x;Hits", 16, -0.5, 15.5);
  hHitFiberY = new TH1D("hHitFiberY", "Hit fiber y;fiber_y;Hits", 16, -0.5, 15.5);

  hRecHitSystem = new TH1D("hRecHitSystem", "RecHit system;system;RecHits", 256, -0.5, 255.5);
  hRecHitSector = new TH1D("hRecHitSector", "RecHit sector;sector;RecHits", 2, -0.5, 1.5);
  hRecHitLayer  = new TH1D("hRecHitLayer",  "RecHit layer;layer;RecHits", 31, -0.5, 30.5);
  hRecHitModule = new TH1D("hRecHitModule", "RecHit module;module;RecHits", 3, -0.5, 2.5);
  hRecHitBlock  = new TH1D("hRecHitBlock",  "RecHit block;block;RecHits", 256, -0.5, 255.5);
  hRecHitFiberX = new TH1D("hRecHitFiberX", "RecHit fiber x;fiber_x;RecHits", 16, -0.5, 15.5);
  hRecHitFiberY = new TH1D("hRecHitFiberY", "RecHit fiber y;fiber_y;RecHits", 16, -0.5, 15.5);

  // Acceptance histograms
  hGenPSCAL = new TH1D("hGenPSCAL", "Generated brems photon;E_{#gamma}^{gen} [GeV];Events",
    100, 0.25, 50.25);

  hAccPSCAL = new TH1D("hAccPSCAL", "Accepted brems photon;E_{#gamma}^{gen} [GeV];Events",
    100, 0.25, 50.25);

  hGenUpCAL = new TH1D("hGenUpCAL", "Generated brems photon;E_{#gamma}^{gen} [GeV];Events",
    100, 0.25, 50.25);

  hAccUpCAL = new TH1D("hAccUpCAL", "Accepted brems photon (top);E_{#gamma}^{gen} [GeV];Events",
    100, 0.25, 50.25);

  hGenDwCAL = new TH1D("hGenDwCAL", "Generated brems photon;E_{#gamma}^{gen} [GeV];Events",
    100, 0.25, 50.25);

  hAccDwCAL = new TH1D("hAccDwCAL", "Accepted brems photon (bottom);E_{#gamma}^{gen} [GeV];Events",
    100, 0.25, 50.25);

  // MC histograms
  hMCGenStatus3_Vx = new TH1D("hMCGenStatus3_Vx", "MC gen particles (status==3);vertex x [mm];Particles", 4000, -200.0, 200.0);
  hMCGenStatus3_Vy = new TH1D("hMCGenStatus3_Vy", "MC gen particles (status==3);vertex y [mm];Particles", 4000, -200.0, 200.0);
  hMCGenStatus3_Vz = new TH1D("hMCGenStatus3_Vz", "MC gen particles (status==3);vertex z [mm];Particles", 7100, -70000.0, 1000.0);

  hMCGenStatus3_Px = new TH1D("hMCGenStatus3_Px", "MC gen particles (status==3);p_{x} [GeV];Particles", 4000, -20.0, 20.0);
  hMCGenStatus3_Py = new TH1D("hMCGenStatus3_Py", "MC gen particles (status==3);p_{y} [GeV];Particles", 4000, -20.0, 20.0);
  hMCGenStatus3_Pz = new TH1D("hMCGenStatus3_Pz", "MC gen particles (status==3);p_{z} [GeV];Particles", 6000, -300.0, 300.0);
  hMCGenStatus3_E = new TH1D("hMCGenStatus3_E", "MC gen particles (status==3);E [GeV];Particles", 6000, 0.0, 300.0);

  hMCGenStatus1_Vx = new TH1D("hMCGenStatus1_Vx", "MC gen particles (status==1);vertex x [mm];Particles", 4000, -200.0, 200.0);
  hMCGenStatus1_Vy = new TH1D("hMCGenStatus1_Vy", "MC gen particles (status==1);vertex y [mm];Particles", 4000, -200.0, 200.0);
  hMCGenStatus1_Vz = new TH1D("hMCGenStatus1_Vz", "MC gen particles (status==1);vertex z [mm];Particles", 7100, -70000.0, 1000.0);

  hMCGenStatus1_Px = new TH1D("hMCGenStatus1_Px", "MC gen particles (status==1);p_{x} [GeV];Particles", 4000, -20.0, 20.0);
  hMCGenStatus1_Py = new TH1D("hMCGenStatus1_Py", "MC gen particles (status==1);p_{y} [GeV];Particles", 4000, -20.0, 20.0);
  hMCGenStatus1_Pz = new TH1D("hMCGenStatus1_Pz", "MC gen particles (status==1);p_{z} [GeV];Particles", 6000, -300.0, 300.0);
  hMCGenStatus1_E = new TH1D("hMCGenStatus1_E", "MC gen particles (status==1);E [GeV];Particles", 6000, 0.0, 300.0);

  hMCGenStatus0_Vx = new TH1D("hMCGenStatus0_Vx", "MC secondary particles (status==0);vertex x [mm];Particles", 4000, -200.0, 200.0);
  hMCGenStatus0_Vy = new TH1D("hMCGenStatus0_Vy", "MC secondary particles (status==0);vertex y [mm];Particles", 4000, -200.0, 200.0);
  hMCGenStatus0_Vz = new TH1D("hMCGenStatus0_Vz", "MC secondary particles (status==0);vertex z [mm];Particles", 7000, -70000.0, 0.0);

  hMCGenStatus0_Px = new TH1D("hMCGenStatus0_Px",
    "MC secondary particles (status==0);p_{x} [GeV];Particles", 4000, -20.0, 20.0);
  hMCGenStatus0_Py = new TH1D("hMCGenStatus0_Py",
    "MC secondary particles (status==0);p_{y} [GeV];Particles", 4000, -20.0, 20.0);
  hMCGenStatus0_Pz = new TH1D("hMCGenStatus0_Pz",
    "MC secondary particles (status==0);p_{z} [GeV];Particles", 6000, -300.0, 300.0);
  hMCGenStatus0_E = new TH1D("hMCGenStatus0_E",
    "MC secondary particles (status==0);E [GeV];Particles", 6000, 0.0, 300.0);

  hEvtPhotonE_Status3 = new TH1D("hEvtPhotonE_Status3",
    "Event generated photon energy (status==3);#Sigma E_{#gamma} [GeV];Events",
    6000, 0.0, 300.0);

  hEvtElectronE_Status3 = new TH1D("hEvtElectronE_Status3",
    "Event generated electron energy (status==3);#Sigma E_{e^{-}} [GeV];Events",
    6000, 0.0, 300.0);

  hEvtPositronE_Status3 = new TH1D("hEvtPositronE_Status3",
    "Event generated positron energy (status==3);#Sigma E_{e^{+}} [GeV];Events",
    6000, 0.0, 300.0);

  hEvtPhotonE_Status1 = new TH1D("hEvtPhotonE_Status1",
    "Event generated photon energy (status==1);#Sigma E_{#gamma} [GeV];Events",
    6000, 0.0, 300.0);

  hEvtElectronE_Status1 = new TH1D("hEvtElectronE_Status1",
    "Event generated electron energy (status==1);#Sigma E_{e^{-}} [GeV];Events",
    6000, 0.0, 300.0);

  hEvtPositronE_Status1 = new TH1D("hEvtPositronE_Status1",
    "Event generated positron energy (status==1);#Sigma E_{e^{+}} [GeV];Events",
    6000, 0.0, 300.0);

  hEvtPhotonE_Status0 = new TH1D("hEvtPhotonE_Status0",
    "Event secondary photon energy (status==0);#Sigma E_{#gamma} [GeV];Events",
    6000, 0.0, 300.0);

  hEvtElectronE_Status0 = new TH1D("hEvtElectronE_Status0",
    "Event secondary electron energy (status==0);#Sigma E_{e^{-}} [GeV];Events",
    6000, 0.0, 300.0);

  hEvtPositronE_Status0 = new TH1D("hEvtPositronE_Status0",
    "Event secondary positron energy (status==0);#Sigma E_{e^{+}} [GeV];Events",
    6000, 0.0, 300.0);

  // Photon endpoint maps
  hPhotonEndX_vs_Z_all = new TH2D("hPhotonEndX_vs_Z_all",
    "All generated photons endpoint;endpoint z [mm];endpoint x [mm]",
    7000, -70000.0, 0.0, 600, -300.0, 300.0);

  hPhotonEndY_vs_Z_all = new TH2D("hPhotonEndY_vs_Z_all",
    "All generated photons endpoint;endpoint z [mm];endpoint y [mm]",
    7000, -70000.0, 0.0, 600, -300.0, 300.0);

  hPhotonEndPx_vs_Z_all = new TH2D("hPhotonEndPx_vs_Z_all",
    "All generated photons momentum at endpoint;endpoint z [mm];p_{x}^{end} [GeV]",
    7000, -70000.0, 0.0, 600, -5.0, 5.0);

  hPhotonEndPy_vs_Z_all = new TH2D("hPhotonEndPy_vs_Z_all",
    "All generated photons momentum at endpoint;endpoint z [mm];p_{y}^{end} [GeV]",
    7000, -70000.0, 0.0, 600, -5.0, 5.0);

  hPhotonEndPz_vs_Z_all = new TH2D("hPhotonEndPz_vs_Z_all",
    "All generated photons momentum at endpoint;endpoint z [mm];p_{z}^{end} [GeV]",
    7000, -70000.0, 0.0, 800, -20.0, 20.0);

  hPhotonEndX_vs_Z_acc  = (TH2D*)hPhotonEndX_vs_Z_all->Clone("hPhotonEndX_vs_Z_acc");
  hPhotonEndY_vs_Z_acc  = (TH2D*)hPhotonEndY_vs_Z_all->Clone("hPhotonEndY_vs_Z_acc");
  hPhotonEndPx_vs_Z_acc = (TH2D*)hPhotonEndPx_vs_Z_all->Clone("hPhotonEndPx_vs_Z_acc");
  hPhotonEndPy_vs_Z_acc = (TH2D*)hPhotonEndPy_vs_Z_all->Clone("hPhotonEndPy_vs_Z_acc");
  hPhotonEndPz_vs_Z_acc = (TH2D*)hPhotonEndPz_vs_Z_all->Clone("hPhotonEndPz_vs_Z_acc");

  hPhotonEndX_vs_Z_lost  = (TH2D*)hPhotonEndX_vs_Z_all->Clone("hPhotonEndX_vs_Z_lost");
  hPhotonEndY_vs_Z_lost  = (TH2D*)hPhotonEndY_vs_Z_all->Clone("hPhotonEndY_vs_Z_lost");
  hPhotonEndPx_vs_Z_lost = (TH2D*)hPhotonEndPx_vs_Z_all->Clone("hPhotonEndPx_vs_Z_lost");
  hPhotonEndPy_vs_Z_lost = (TH2D*)hPhotonEndPy_vs_Z_all->Clone("hPhotonEndPy_vs_Z_lost");
  hPhotonEndPz_vs_Z_lost = (TH2D*)hPhotonEndPz_vs_Z_all->Clone("hPhotonEndPz_vs_Z_lost");

  hPhotonEndX_vs_Z_acc->SetTitle("Accepted photons endpoint;endpoint z [mm];endpoint x [mm]");
  hPhotonEndY_vs_Z_acc->SetTitle("Accepted photons endpoint;endpoint z [mm];endpoint y [mm]");
  hPhotonEndX_vs_Z_lost->SetTitle("Lost photons endpoint;endpoint z [mm];endpoint x [mm]");
  hPhotonEndY_vs_Z_lost->SetTitle("Lost photons endpoint;endpoint z [mm];endpoint y [mm]");

  hPhotonAcceptedZ = new TH1D("hPhotonAcceptedZ",
    "Accepted generated photons;endpoint z [mm];Photons",
    7000, -70000.0, 0.0);

  hPhotonLostZ = new TH1D("hPhotonLostZ",
    "Lost generated photons;endpoint z [mm];Photons",
    7000, -70000.0, 0.0);

  // Detector ranges
  detectorRanges = {
    {"Unknown1",         "Unknown1",        -4800.0,  -3400.0,  true,  false},
    {"Exit window",      "ExitWindow",      -18510.0, -18500.0, true,  true},
    {"Unknown2",         "Unknown2",        -25350.0, -25250.0, true,  false},
    {"Lumi wall",        "LumiWall",        -55100.0, -54900.0, true,  false},
    {"Sweeper magnet",   "SweeperMagnet",   -56766.5, -55233.5, true,  false},
    {"Analyser magnet",  "AnalyserMagnet",  -60766.5, -59233.5, true,  false},
    {"Photon exit cap",  "PhotonExitCap",   -60769.0, -60764.0, true,  false},
    {"Photon entry cap", "PhotonEntryCap",  -55238.5, -55228.5, true,  true},
    {"Thin film",        "ThinFilm",        -58001.0, -57999.0, true,  true},
    {"Air column",       "AirColumn",       -55228.5, -18510.0, true,  false}
  };

  hLostDet.clear();
  hAccDet.clear();

  for (const auto& det : detectorRanges) {
    if (det.saveLost)
      hLostDet[det.key] = MakePhotonDetectorHistSet("Lost", det, "Lost");

    if (det.saveAccepted)
      hAccDet[det.key] = MakePhotonDetectorHistSet("Accepted", det, "Accepted");
  }

  //==================================================
  // Photon position reconstruction: centroid method
  //==================================================

  TString name = "photon_pos_resol_centroid";

  hDiffX_Ctd = new TH1D(
      Form("hDiffX_%s", name.Data()),
      "Photon X resolution using centroid method;"
      "X_{rec}-X_{true} [cm];Events",
      200, -5.0, 5.0
  );

  hDiffY_Ctd = new TH1D(
      Form("hDiffY_%s", name.Data()),
      "Photon Y resolution using centroid method;"
      "Y_{rec}-Y_{true} [cm];Events",
      200, -5.0, 5.0
  );

  hRecPhotonAtThinFilmVx_Ctd = new TH1D(
      Form("hRecPhotonAtThinFilmVx_%s", name.Data()),
      "Reconstructed photon X position at thin film using centroid method;"
      "X_{rec} [cm];Events",
      400, -10.0, 10.0
  );

  hRecPhotonAtThinFilmVy_Ctd = new TH1D(
      Form("hRecPhotonAtThinFilmVy_%s", name.Data()),
      "Reconstructed photon Y position at thin film using centroid method;"
      "Y_{rec} [cm];Events",
      400, -10.0, 10.0
  );

  hGenPhotonAtThinFilmVx_Ctd = new TH1D(
      Form("hGenPhotonAtThinFilmVx_%s", name.Data()),
      "Generated photon X position at thin film for centroid method;"
      "X_{true} [cm];Events",
      400, -10.0, 10.0
  );

  hGenPhotonAtThinFilmVy_Ctd = new TH1D(
      Form("hGenPhotonAtThinFilmVy_%s", name.Data()),
      "Generated photon Y position at thin film for centroid method;"
      "Y_{true} [cm];Events",
      400, -10.0, 10.0
  );

  //===============================================
  // Photon position reconstruction: slope method
  //===============================================

  name = "photon_pos_resol_slope";

  hDiffX_Slope = new TH1D(
      Form("hDiffX_%s", name.Data()),
      "Photon X resolution using slope method;"
      "X_{rec}-X_{true} [cm];Events",
      200, -5.0, 5.0
  );

  hDiffY_Slope = new TH1D(
      Form("hDiffY_%s", name.Data()),
      "Photon Y resolution using slope method;"
      "Y_{rec}-Y_{true} [cm];Events",
      200, -5.0, 5.0
  );

  hRecPhotonAtThinFilmVx_Slope = new TH1D(
      Form("hRecPhotonAtThinFilmVx_%s", name.Data()),
      "Reconstructed photon X position at thin film using slope method;"
      "X_{rec} [cm];Events",
      400, -10.0, 10.0
  );

  hRecPhotonAtThinFilmVy_Slope = new TH1D(
      Form("hRecPhotonAtThinFilmVy_%s", name.Data()),
      "Reconstructed photon Y position at thin film using slope method;"
      "Y_{rec} [cm];Events",
      400, -10.0, 10.0
  );

  hGenPhotonAtThinFilmVx_Slope = new TH1D(
      Form("hGenPhotonAtThinFilmVx_%s", name.Data()),
      "Generated photon X position at thin film for slope method;"
      "X_{true} [cm];Events",
      400, -10.0, 10.0
  );

  hGenPhotonAtThinFilmVy_Slope = new TH1D(
      Form("hGenPhotonAtThinFilmVy_%s", name.Data()),
      "Generated photon Y position at thin film for slope method;"
      "Y_{true} [cm];Events",
      400, -10.0, 10.0
  );

  //===============================================
  // Photon position reconstruction: tracker method
  //===============================================

  name = "photon_pos_resol_tracker";

  hDiffX_Trk = new TH1D(
      Form("hDiffX_%s", name.Data()),
      "Photon X resolution using Trk method;"
      "X_{rec}-X_{true} [cm];Events",
      200, -5.0, 5.0
  );

  hDiffY_Trk = new TH1D(
      Form("hDiffY_%s", name.Data()),
      "Photon Y resolution using Trk method;"
      "Y_{rec}-Y_{true} [cm];Events",
      200, -5.0, 5.0
  );

  hRecPhotonAtThinFilmVx_Trk = new TH1D(
      Form("hRecPhotonAtThinFilmVx_%s", name.Data()),
      "Reconstructed photon X position at thin film using Trk method;"
      "X_{rec} [cm];Events",
      400, -10.0, 10.0
  );

  hRecPhotonAtThinFilmVy_Trk = new TH1D(
      Form("hRecPhotonAtThinFilmVy_%s", name.Data()),
    "Reconstructed photon Y position at thin film using Trk method;"
    "Y_{rec} [cm];Events",
    400, -10.0, 10.0
  );

  hGenPhotonAtThinFilmVx_Trk = new TH1D(
      Form("hGenPhotonAtThinFilmVx_%s", name.Data()),
      "Generated photon X position at thin film for Trk method;"
      "X_{true} [cm];Events",
      400, -10.0, 10.0
  );

  hGenPhotonAtThinFilmVy_Trk = new TH1D(
      Form("hGenPhotonAtThinFilmVy_%s", name.Data()),
      "Generated photon Y position at thin film for Trk method;"
      "Y_{true} [cm];Events",
      400, -10.0, 10.0
  );

  //=====================================
  // Electron/Positron Endpoint Check
  //=====================================

  hLeptonEndX = new TH1D(
    "hLeptonEndX",
    "e^{#pm} endpoint X;endpoint x [mm];Particles",
    600, -300.0, 300.0
  );

  hLeptonEndY = new TH1D(
      "hLeptonEndY",
      "e^{#pm} endpoint Y;endpoint y [mm];Particles",
      600, -300.0, 300.0
  );

  hLeptonEndZ = new TH1D(
      "hLeptonEndZ",
      "e^{#pm} endpoint Z;endpoint z [mm];Particles",
      1000, -70000.0, 0.0
  );

  hLeptonPDG_vs_EndY = new TH2D(
      "hLeptonPDG_vs_EndY",
      "e^{#pm} endpoint Y vs PDG;PDG;endpoint y [mm]",
      23, -11.5, 11.5,
      600, -300.0, 300.0
  );

  //==================================================
  // Electron position reconstruction: centroid method
  //==================================================

  name = "cal_pos_resol_centroid";

  for (int sec = 0; sec < NSector; ++sec) {

      TString sectorName = (sec == 0) ? "top" : "bot";

      hDiffElectronX_Ctd[sec] = new TH1D(
          Form("hDiffElectronX_%s_%s", name.Data(), sectorName.Data()),
          Form(
              "Electron X resolution using centroid method (%s);"
              "X_{rec}-X_{true} [cm];Events",
              sectorName.Data()
          ),
          200, -5.0, 5.0
      );

      hDiffElectronY_Ctd[sec] = new TH1D(
          Form("hDiffElectronY_%s_%s", name.Data(), sectorName.Data()),
          Form(
              "Electron Y resolution using centroid method (%s);"
              "Y_{rec}-Y_{true} [cm];Events",
              sectorName.Data()
          ),
          200, -5.0, 5.0
      );

      hRecElectronVx_Ctd[sec] = new TH1D(
          Form("hRecElectronVx_%s_%s", name.Data(), sectorName.Data()),
          Form(
              "Reconstructed electron X position using centroid method (%s);"
              "X_{rec} [cm];Events",
              sectorName.Data()
          ),
          800, -40.0, 40.0
      );

      hRecElectronVy_Ctd[sec] = new TH1D(
          Form("hRecElectronVy_%s_%s", name.Data(), sectorName.Data()),
          Form(
              "Reconstructed electron Y position using centroid method (%s);"
              "Y_{rec} [cm];Events",
              sectorName.Data()
          ),
          800, -40.0, 40.0
      );

      hGenElectronVx_Ctd[sec] = new TH1D(
          Form("hGenElectronVx_%s_%s", name.Data(), sectorName.Data()),
          Form(
              "Generated electron X position for centroid method (%s);"
              "X_{true} [cm];Events",
              sectorName.Data()
          ),
          800, -40.0, 40.0
      );

      hGenElectronVy_Ctd[sec] = new TH1D(
          Form("hGenElectronVy_%s_%s", name.Data(), sectorName.Data()),
          Form(
              "Generated electron Y position for centroid method (%s);"
              "Y_{true} [cm];Events",
              sectorName.Data()
          ),
          800, -40.0, 40.0
      );
  }


  //===============================================
  // Electron position reconstruction: slope method
  //===============================================

  name = "cal_pos_resol_slope";

  for (int sec = 0; sec < NSector; ++sec) {

      TString sectorName = (sec == 0) ? "top" : "bot";

      hDiffElectronX_Slope[sec] = new TH1D(
          Form("hDiffElectronX_%s_%s", name.Data(), sectorName.Data()),
          Form(
              "Electron X resolution using slope method (%s);"
              "X_{rec}-X_{tracker} [cm];Events",
              sectorName.Data()
          ),
          200, -5.0, 5.0
      );

      hDiffElectronY_Slope[sec] = new TH1D(
          Form("hDiffElectronY_%s_%s", name.Data(), sectorName.Data()),
          Form(
              "Electron Y resolution using slope method (%s);"
              "Y_{rec}-Y_{tracker} [cm];Events",
              sectorName.Data()
          ),
          200, -5.0, 5.0
      );

      hRecElectronVx_Slope[sec] = new TH1D(
          Form("hRecElectronVx_%s_%s", name.Data(), sectorName.Data()),
          Form(
              "Reconstructed electron X position using slope method (%s);"
              "X_{rec} [cm];Events",
              sectorName.Data()
          ),
          800, -40.0, 40.0
      );

      hRecElectronVy_Slope[sec] = new TH1D(
          Form("hRecElectronVy_%s_%s", name.Data(), sectorName.Data()),
          Form(
              "Reconstructed electron Y position using slope method (%s);"
              "Y_{rec} [cm];Events",
              sectorName.Data()
          ),
          800, -40.0, 40.0
      );

      hGenElectronVx_Slope[sec] = new TH1D(
          Form("hGenElectronVx_%s_%s", name.Data(), sectorName.Data()),
          Form(
              "Tracker electron X position for slope method (%s);"
              "X_{tracker} [cm];Events",
              sectorName.Data()
          ),
          800, -40.0, 40.0
      );

      hGenElectronVy_Slope[sec] = new TH1D(
          Form("hGenElectronVy_%s_%s", name.Data(), sectorName.Data()),
          Form(
              "Tracker electron Y position for slope method (%s);"
              "Y_{tracker} [cm];Events",
              sectorName.Data()
          ),
          800, -40.0, 40.0
      );
  }

  //==================================================
  // 2D Electron position reconstruction: centroid
  //==================================================

  name = "cal_pos_energy_corr_centroid";

  for (int sec = 0; sec < NSector; ++sec) {

      TString sectorName = (sec == top) ? "top" : "bot";

      hDiffElectronX_vs_gE_Ctd[sec] = new TH2D(
          Form("hDiffElectronX_vs_gE_%s_%s",
              name.Data(), sectorName.Data()),
          Form(
              "Electron X residual vs generated energy, centroid method (%s);"
              "E_{gen} [GeV];X_{rec}-X_{trk} [cm]",
              sectorName.Data()
          ),
          4000, 0.0, 40.0,
          200, -5.0, 5.0
      );

      hDiffElectronY_vs_gE_Ctd[sec] = new TH2D(
          Form("hDiffElectronY_vs_gE_%s_%s",
              name.Data(), sectorName.Data()),
          Form(
              "Electron Y residual vs generated energy, centroid method (%s);"
              "E_{gen} [GeV];Y_{rec}-Y_{trk} [cm]",
              sectorName.Data()
          ),
          4000, 0.0, 40.0,
          200, -5.0, 5.0
      );

      hRecElectronVx_vs_gE_Ctd[sec] = new TH2D(
          Form("hRecElectronVx_vs_gE_%s_%s",
              name.Data(), sectorName.Data()),
          Form(
              "Reconstructed electron X vs generated energy, centroid method (%s);"
              "E_{gen} [GeV];X_{rec} [cm]",
              sectorName.Data()
          ),
          4000, 0.0, 40.0,
          800, -40.0, 40.0
      );

      hRecElectronVy_vs_gE_Ctd[sec] = new TH2D(
          Form("hRecElectronVy_vs_gE_%s_%s",
              name.Data(), sectorName.Data()),
          Form(
              "Reconstructed electron Y vs generated energy, centroid method (%s);"
              "E_{gen} [GeV];Y_{rec} [cm]",
              sectorName.Data()
          ),
          4000, 0.0, 40.0,
          800, -40.0, 40.0
      );

      hTrkHitx_vs_gE_Ctd[sec] = new TH2D(
          Form("hTrkHitx_vs_gE_%s_%s",
              name.Data(), sectorName.Data()),
          Form(
              "Tracker hit X vs generated energy, centroid method (%s);"
              "E_{gen} [GeV];X_{trk} [cm]",
              sectorName.Data()
          ),
          4000, 0.0, 40.0,
          800, -40.0, 40.0
      );

      hTrkHity_vs_gE_Ctd[sec] = new TH2D(
          Form("hTrkHity_vs_gE_%s_%s",
              name.Data(), sectorName.Data()),
          Form(
              "Tracker hit Y vs generated energy, centroid method (%s);"
              "E_{gen} [GeV];Y_{trk} [cm]",
              sectorName.Data()
          ),
          4000, 0.0, 40.0,
          800, -40.0, 40.0
      );
  }


  //===============================================
  // 2D Electron position reconstruction: slope
  //===============================================

  name = "cal_pos_energy_corr_slope";

  for (int sec = 0; sec < NSector; ++sec) {

      TString sectorName = (sec == top) ? "top" : "bot";

      hDiffElectronX_vs_gE_Slope[sec] = new TH2D(
          Form("hDiffElectronX_vs_gE_%s_%s",
              name.Data(), sectorName.Data()),
          Form(
              "Electron X residual vs generated energy, slope method (%s);"
              "E_{gen} [GeV];X_{rec}-X_{trk} [cm]",
              sectorName.Data()
          ),
          4000, 0.0, 40.0,
          200, -5.0, 5.0
      );

      hDiffElectronY_vs_gE_Slope[sec] = new TH2D(
          Form("hDiffElectronY_vs_gE_%s_%s",
              name.Data(), sectorName.Data()),
          Form(
              "Electron Y residual vs generated energy, slope method (%s);"
              "E_{gen} [GeV];Y_{rec}-Y_{trk} [cm]",
              sectorName.Data()
          ),
          4000, 0.0, 40.0,
          200, -5.0, 5.0
      );

      hRecElectronVx_vs_gE_Slope[sec] = new TH2D(
          Form("hRecElectronVx_vs_gE_%s_%s",
              name.Data(), sectorName.Data()),
          Form(
              "Reconstructed electron X vs generated energy, slope method (%s);"
              "E_{gen} [GeV];X_{rec} [cm]",
              sectorName.Data()
          ),
          4000, 0.0, 40.0,
          800, -40.0, 40.0
      );

      hRecElectronVy_vs_gE_Slope[sec] = new TH2D(
          Form("hRecElectronVy_vs_gE_%s_%s",
              name.Data(), sectorName.Data()),
          Form(
              "Reconstructed electron Y vs generated energy, slope method (%s);"
              "E_{gen} [GeV];Y_{rec} [cm]",
              sectorName.Data()
          ),
          4000, 0.0, 40.0,
          800, -40.0, 40.0
      );

      hTrkHitx_vs_gE_Slope[sec] = new TH2D(
          Form("hTrkHitx_vs_gE_%s_%s",
              name.Data(), sectorName.Data()),
          Form(
              "Tracker hit X vs generated energy, slope method (%s);"
              "E_{gen} [GeV];X_{trk} [cm]",
              sectorName.Data()
          ),
          4000, 0.0, 40.0,
          800, -40.0, 40.0
      );

      hTrkHity_vs_gE_Slope[sec] = new TH2D(
          Form("hTrkHity_vs_gE_%s_%s",
              name.Data(), sectorName.Data()),
          Form(
              "Tracker hit Y vs generated energy, slope method (%s);"
              "E_{gen} [GeV];Y_{trk} [cm]",
              sectorName.Data()
          ),
          4000, 0.0, 40.0,
          800, -40.0, 40.0
      );
  }

  //===============================================
  // 2D Electron energy reconstruction: centroid
  //===============================================
  name = "cal_energy_resol_centroid";
  for (int sec = 0; sec < NSector; ++sec) {

      TString sectorName = (sec == top) ? "top" : "bot";

      hRecEnergy_vs_gE_Ctd[sec] = new TH2D(
          Form(
              "hRecEnergy_vs_gE_%s_%s",
              name.Data(),
              sectorName.Data()
          ),
          Form(
              "Reconstructed vs generated electron energy, centroid method (%s);"
              "E_{gen} [GeV];E_{rec} [GeV]",
              sectorName.Data()
          ),
          200, 0.0, 20.0,
          200, 0.0, 20.0
      );
  }

  //===============================================
  // 2D Electron energy reconstruction: slope
  //===============================================
  name = "cal_energy_resol_slope";
  for (int sec = 0; sec < NSector; ++sec) {

      TString sectorName = (sec == top) ? "top" : "bot";

      hRecEnergy_vs_gE_Slope[sec] = new TH2D(
          Form(
              "hRecEnergy_vs_gE_%s_%s",
              name.Data(),
              sectorName.Data()
          ),
          Form(
              "Reconstructed vs generated electron energy, slope method (%s);"
              "E_{gen} [GeV];E_{rec} [GeV]",
              sectorName.Data()
          ),
          200, 0.0, 20.0,
          200, 0.0, 20.0
      );
  }

  //==========================
  // Reconstructed Tracker Hits
  //==========================

  hTrkRecHitEnergy = new TH1D(
      "hTrkRecHitEnergy",
      "Tracker Hit Energy Deposition;E_{dep} (GeV);Entries",
      200, 0.0, 0.02);

  hTrkRecHitTime = new TH1D(
      "hTrkRecHitTime",
      "Tracker Hit Time;Time (ns);Entries",
      500, 0.0, 500.0);

  hTrkRecHitPrimary = new TH1D(
      "hTrkRecHitPrimary",
      "Tracker Hit Origin;Primary Hit;Entries",
      2, -0.5, 1.5);

  hTrkRecHitSystem = new TH1D(
      "hTrkRecHitSystem",
      "Tracker Hit System ID;System ID;Entries",
      10, -0.5, 9.5);

  hTrkRecHitSector = new TH1D(
      "hTrkRecHitSector",
      "Tracker Hit Sector ID;Sector ID;Entries",
      10, -0.5, 9.5);

  hTrkRecHitModule = new TH1D(
      "hTrkRecHitModule",
      "Tracker Hit Module ID;Module ID;Entries",
      20, -0.5, 19.5);

  hTrkRecHitXId = new TH1D(
      "hTrkRecHitXId",
      "Tracker Hit X Pixel ID;X Pixel Index;Entries",
      400, -200.5, 199.5);

  hTrkRecHitYId = new TH1D(
      "hTrkRecHitYId",
      "Tracker Hit Y Pixel ID;Y Pixel Index;Entries",
      400, -200.5, 199.5);

  hTrkRecHitGlobalX = new TH1D(
      "hTrkRecHitGlobalX",
      "Tracker Hit Global X;Global X (mm);Entries",
      500, -250.0, 250.0);

  hTrkRecHitGlobalY = new TH1D(
      "hTrkRecHitGlobalY",
      "Tracker Hit Global Y;Global Y (mm);Entries",
      500, -250.0, 250.0);

  hTrkRecHitGlobalZ = new TH1D(
      "hTrkRecHitGlobalZ",
      "Tracker Hit Global Z;Global Z (mm);Entries",
      1000, -1000.0, 1000.0);

  hTrkRecHitGlobalXY = new TH2D(
      "hTrkRecHitGlobalXY",
      "Tracker Hit Global XY;Global X (mm);Global Y (mm)",
      500, -250.0, 250.0,
      500, -250.0, 250.0);

  hTrkRecHitGlobalXZ = new TH2D(
      "hTrkRecHitGlobalXZ",
      "Tracker Hit Global XZ;Global X (mm);Global Z (mm)",
      500, -250.0, 250.0,
      1000, -1000.0, 1000.0);

}

//_________________________________________________________________________________________________________
void SaveInFile(TFile* fout){

  fout->cd();

  TDirectory* dirEvt      = fout->mkdir("Event");
  TDirectory* dirMC       = fout->mkdir("MCParticles");
  TDirectory* dirHit      = fout->mkdir("GeantHits");
  TDirectory* dirRec      = fout->mkdir("RecHits");
  TDirectory* dirAccept   = fout->mkdir("Acceptance");
  TDirectory* dirLoss     = fout->mkdir("PhotonLossMap");
  TDirectory* dirPosResol = fout->mkdir("PositionResolution");
  TDirectory* electronCentroidDir = dirPosResol->mkdir("ElectronPosResol_Centroid");
  TDirectory* electronSlopeDir    = dirPosResol->mkdir("ElectronPosResol_Slope");
  TDirectory* photonCentroidDir   = dirPosResol->mkdir("PhotonPosResol_Centroid");
  TDirectory* photonSlopeDir      = dirPosResol->mkdir("PhotonPosResol_Slope");
  TDirectory* photonTrkDir      = dirPosResol->mkdir("PhotonPosResol_Tracker");
  TDirectory* dirTree     = fout->mkdir("Tree");

  dirEvt->cd();
  hEvtHitE->Write();
  hEvtRecE_top->Write();
  hEvtRecE_bot->Write();
  hEvtRecE_tot->Write();
  hEvtDiff_hits_minus_rec->Write();
  hSamplingFraction_hits->Write();
  hSamplingFraction_rec->Write();

  dirMC->cd();
  hMCGenStatus3_Vx->Write();
  hMCGenStatus3_Vy->Write();
  hMCGenStatus3_Vz->Write();
  hMCGenStatus3_Px->Write();
  hMCGenStatus3_Py->Write();
  hMCGenStatus3_Pz->Write();
  hMCGenStatus3_E->Write();

  hMCGenStatus1_Vx->Write();
  hMCGenStatus1_Vy->Write();
  hMCGenStatus1_Vz->Write();
  hMCGenStatus1_Px->Write();
  hMCGenStatus1_Py->Write();
  hMCGenStatus1_Pz->Write();
  hMCGenStatus1_E->Write();

  hMCGenStatus0_Vx->Write();
  hMCGenStatus0_Vy->Write();
  hMCGenStatus0_Vz->Write();
  hMCGenStatus0_Px->Write();
  hMCGenStatus0_Py->Write();
  hMCGenStatus0_Pz->Write();
  hMCGenStatus0_E->Write();

  hEvtPhotonE_Status3->Write();
  hEvtElectronE_Status3->Write();
  hEvtPositronE_Status3->Write();

  hEvtPhotonE_Status1->Write();
  hEvtElectronE_Status1->Write();
  hEvtPositronE_Status1->Write();

  hEvtPhotonE_Status0->Write();
  hEvtElectronE_Status0->Write();
  hEvtPositronE_Status0->Write();

  dirHit->cd();
  hHitSystem->Write();
  hHitSector->Write();
  hHitLayer->Write();
  hHitModule->Write();
  hHitBlock->Write();
  hHitFiberX->Write();
  hHitFiberY->Write();

  dirRec->cd();
  hRecHitSystem->Write();
  hRecHitSector->Write();
  hRecHitLayer->Write();
  hRecHitModule->Write();
  hRecHitBlock->Write();
  hRecHitFiberX->Write();
  hRecHitFiberY->Write();

  dirAccept->cd();
  hGenPSCAL->Write();
  hGenUpCAL->Write();
  hGenDwCAL->Write();
  hAccPSCAL->Write();
  hAccUpCAL->Write();
  hAccDwCAL->Write();
  hAcceptancePSCAL->Write();
  hAcceptanceUpCAL->Write();
  hAcceptanceDwCAL->Write();

  dirPosResol->cd();

  hLeptonEndX->Write();
  hLeptonEndY->Write();
  hLeptonEndZ->Write();
  hLeptonPDG_vs_EndY->Write();

  //==============================
  // Centroid method
  //==============================
  electronCentroidDir->cd();
  for (int sec = 0; sec < NSector; ++sec) {
      // 1D
      if (hDiffElectronX_Ctd[sec])
          hDiffElectronX_Ctd[sec]->Write();
      if (hDiffElectronY_Ctd[sec])
          hDiffElectronY_Ctd[sec]->Write();
      if (hRecElectronVx_Ctd[sec])
          hRecElectronVx_Ctd[sec]->Write();
      if (hRecElectronVy_Ctd[sec])
          hRecElectronVy_Ctd[sec]->Write();
      if (hGenElectronVx_Ctd[sec])
          hGenElectronVx_Ctd[sec]->Write();
      if (hGenElectronVy_Ctd[sec])
          hGenElectronVy_Ctd[sec]->Write();
      // 2D
      if (hDiffElectronX_vs_gE_Ctd[sec])
          hDiffElectronX_vs_gE_Ctd[sec]->Write();
      if (hDiffElectronY_vs_gE_Ctd[sec])
          hDiffElectronY_vs_gE_Ctd[sec]->Write();
      if (hRecElectronVx_vs_gE_Ctd[sec])
          hRecElectronVx_vs_gE_Ctd[sec]->Write();
      if (hRecElectronVy_vs_gE_Ctd[sec])
          hRecElectronVy_vs_gE_Ctd[sec]->Write();
      if (hTrkHitx_vs_gE_Ctd[sec])
          hTrkHitx_vs_gE_Ctd[sec]->Write();
      if (hTrkHity_vs_gE_Ctd[sec])
          hTrkHity_vs_gE_Ctd[sec]->Write();
      if(hRecEnergy_vs_gE_Ctd[sec])
          hRecEnergy_vs_gE_Ctd[sec]->Write();
  }

  //==============================
  // Slope method
  //==============================
  electronSlopeDir->cd();
  for (int sec = 0; sec < NSector; ++sec) {
      // 1D
      if (hDiffElectronX_Slope[sec])
          hDiffElectronX_Slope[sec]->Write();
      if (hDiffElectronY_Slope[sec])
          hDiffElectronY_Slope[sec]->Write();
      if (hRecElectronVx_Slope[sec])
          hRecElectronVx_Slope[sec]->Write();
      if (hRecElectronVy_Slope[sec])
          hRecElectronVy_Slope[sec]->Write();
      if (hGenElectronVx_Slope[sec])
          hGenElectronVx_Slope[sec]->Write();
      if (hGenElectronVy_Slope[sec])
          hGenElectronVy_Slope[sec]->Write();
      // 2D
      if (hDiffElectronX_vs_gE_Slope[sec])
          hDiffElectronX_vs_gE_Slope[sec]->Write();
      if (hDiffElectronY_vs_gE_Slope[sec])
          hDiffElectronY_vs_gE_Slope[sec]->Write();
      if (hRecElectronVx_vs_gE_Slope[sec])
          hRecElectronVx_vs_gE_Slope[sec]->Write();
      if (hRecElectronVy_vs_gE_Slope[sec])
          hRecElectronVy_vs_gE_Slope[sec]->Write();
      if (hTrkHitx_vs_gE_Slope[sec])
          hTrkHitx_vs_gE_Slope[sec]->Write();
      if (hTrkHity_vs_gE_Slope[sec])
          hTrkHity_vs_gE_Slope[sec]->Write();
      if(hRecEnergy_vs_gE_Slope[sec])
          hRecEnergy_vs_gE_Slope[sec]->Write();
  }

  photonCentroidDir->cd();
  hDiffX_Ctd->Write();
  hDiffY_Ctd->Write();
  hRecPhotonAtThinFilmVx_Ctd->Write();
  hRecPhotonAtThinFilmVy_Ctd->Write();
  hGenPhotonAtThinFilmVx_Ctd->Write();
  hGenPhotonAtThinFilmVy_Ctd->Write();

  photonSlopeDir->cd();
  hDiffX_Slope->Write();
  hDiffY_Slope->Write();
  hRecPhotonAtThinFilmVx_Slope->Write();
  hRecPhotonAtThinFilmVy_Slope->Write();
  hGenPhotonAtThinFilmVx_Slope->Write();
  hGenPhotonAtThinFilmVy_Slope->Write();

  photonTrkDir->cd();
  hDiffX_Trk->Write();
  hDiffY_Trk->Write();
  hRecPhotonAtThinFilmVx_Trk->Write();
  hRecPhotonAtThinFilmVy_Trk->Write();
  hGenPhotonAtThinFilmVx_Trk->Write();
  hGenPhotonAtThinFilmVy_Trk->Write();

  dirLoss->cd();
  hPhotonEndX_vs_Z_all->Write();
  hPhotonEndY_vs_Z_all->Write();
  hPhotonEndPx_vs_Z_all->Write();
  hPhotonEndPy_vs_Z_all->Write();
  hPhotonEndPz_vs_Z_all->Write();

  hPhotonEndX_vs_Z_acc->Write();
  hPhotonEndY_vs_Z_acc->Write();
  hPhotonEndPx_vs_Z_acc->Write();
  hPhotonEndPy_vs_Z_acc->Write();
  hPhotonEndPz_vs_Z_acc->Write();

  hPhotonEndX_vs_Z_lost->Write();
  hPhotonEndY_vs_Z_lost->Write();
  hPhotonEndPx_vs_Z_lost->Write();
  hPhotonEndPy_vs_Z_lost->Write();
  hPhotonEndPz_vs_Z_lost->Write();

  hPhotonAcceptedZ->Write();
  hPhotonLostZ->Write();

  TDirectory* dirLostDet = dirLoss->mkdir("DetectorWiseLostPhotons");
  dirLostDet->cd();

  for (const auto& item : hLostDet) {
    TDirectory* d = dirLostDet->mkdir(item.first);
    d->cd();

    item.second.hEndY_vs_EndX->Write();
    item.second.hGenY_vs_GenX->Write();
    item.second.hGenZ->Write();
    item.second.hGenE->Write();

    dirLostDet->cd();
  }

  TDirectory* dirAccDet = dirLoss->mkdir("DetectorWiseAcceptedPhotons");
  dirAccDet->cd();

  for (const auto& item : hAccDet) {
    TDirectory* d = dirAccDet->mkdir(item.first);
    d->cd();

    item.second.hEndY_vs_EndX->Write();
    item.second.hGenY_vs_GenX->Write();
    item.second.hGenZ->Write();
    item.second.hGenE->Write();

    dirAccDet->cd();
  }

  // dirTree->cd();
  // outTree->Write(); // ***** NO NEED TO SAVE TREE LEVEL INFO FOR BIG DATA *****

  fout->Close();
}