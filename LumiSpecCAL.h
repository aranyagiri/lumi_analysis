# Author: Aranya Giri
# Description: Parameters used in LumiSpecCAL_GETaLMAnalysis.C code.

#ifndef LUMISPECCAL_H
#define LUMISPECCAL_H

#include <map>
#include <vector>
#include <string>
#include <cstdint>
#include <iostream>
#include <fstream>
#include <cmath>

#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "TFile.h"
#include "TVector3.h"

enum sector_t { top = 0, bot = 1};
const int NSector = 2;
const int NLayer  = 20;
const int NModule = 3;
const int NBlockX = 14;
const int NBlockY = 2;
const int NSlice  = NLayer*NBlockY;
const int NTracker = 2; // two tracking walls, can be 1 or 3. 

// Geometry
// root -l -q 'PrintLumiConstants.C("/opt/detector/epic-26.06.0/share/epic/compact/far_backward/definitions.xml")'
const double block_xy_size   = 0.43;
//const double mod_coat_size   = 0.02; It is only present, in my newest aranya branch.
const double mod_coat_size   = 0.0;
const double lay_coat_size   = 0.25;

const double mod_x_size      = 14*block_xy_size + 2*mod_coat_size;
const double mod_y_size      = 2*block_xy_size + 2*mod_coat_size;

const double lay_x_size      = (3.0 * mod_x_size) + lay_coat_size * 2.0;
const double lay_y_size      = mod_y_size + lay_coat_size * 2.0;

const double cal_dxy         = lay_x_size;
const double cal_dz          = 20*lay_y_size;

const double cal_x_low_pos          = - (cal_dxy /2.0);
const double fivesigma_gap_atCAL    =  6.766348; //from y=0
const double cal_y_low_pos[NSector] = {fivesigma_gap_atCAL, - fivesigma_gap_atCAL - cal_dxy}; 
const double cal_z_mid_pos          = -6414.0; 
const double cal_z_front_pos        = -6414.0 + (cal_dz /2.0); 

//variable initialization
Long64_t event_id = 0;

std::vector<int> gen_pdg_vec, gen_status_vec;
std::vector<float> gen_charge_vec;
std::vector<float> gen_px_vec, gen_py_vec, gen_pz_vec;
std::vector<float> gen_mass_vec, gen_energy_vec;
std::vector<float> gen_vx_vec, gen_vy_vec, gen_vz_vec;

// CAL RecHits
std::vector<float> cal_rec_energy_vec;
std::vector<unsigned long long> cal_rec_cellID_vec;

std::vector<int> cal_rec_system_vec;
std::vector<int> cal_rec_sector_vec;
std::vector<int> cal_rec_layer_vec;
std::vector<int> cal_rec_module_vec;
std::vector<int> cal_rec_block_vec;
std::vector<int> cal_rec_fiberx_vec;
std::vector<int> cal_rec_fibery_vec;

// Tracker RecHits
/*
std::vector<unsigned long long> trk_rec_cellID_vec;

std::vector<int> trk_rec_system_vec;
std::vector<int> trk_rec_sector_vec;
std::vector<int> trk_rec_module_vec;
std::vector<int> trk_rec_x_id_vec;
std::vector<int> trk_rec_y_id_vec;

std::vector<double> trk_rec_glob_x_vec;
std::vector<double> trk_rec_glob_y_vec;
std::vector<double> trk_rec_glob_z_vec;

std::vector<float> trk_rec_energy_vec;
std::vector<int> trk_rec_quality_vec;
std::vector<float> trk_rec_time_vec;
*/

std::vector<unsigned long long> trk_hits_cellID_vec;

std::vector<int> trk_hits_system_vec;
std::vector<int> trk_hits_sector_vec;
std::vector<int> trk_hits_module_vec;
std::vector<int> trk_hits_x_id_vec;
std::vector<int> trk_hits_y_id_vec;

std::vector<double> trk_hits_glob_x_vec;
std::vector<double> trk_hits_glob_y_vec;
std::vector<double> trk_hits_glob_z_vec;

std::vector<float> trk_hits_energy_vec;
std::vector<int> trk_hits_quality_vec;
std::vector<float> trk_hits_time_vec;

double total_gen_energy_all_status1 = 0.0;
double total_gen_energy_all_status3 = 0.0;

double total_hit_energy = 0.0;
double total_rec_energy = 0.0;
double total_rec_energy_top = 0.0;
double total_rec_energy_bot = 0.0;

double evt_genE_photon_status3 = 0.0;
double evt_genE_electron_status3 = 0.0;
double evt_genE_positron_status3 = 0.0;

double evt_genE_photon_status1 = 0.0;
double evt_genE_electron_status1 = 0.0;
double evt_genE_positron_status1 = 0.0;

double evt_genE_photon_status0 = 0.0;
double evt_genE_electron_status0 = 0.0;
double evt_genE_positron_status0 = 0.0;

// Event histograms
TH1D* hEvtHitE = nullptr;
TH1D* hEvtRecE_top = nullptr;
TH1D* hEvtRecE_bot = nullptr;
TH1D* hEvtRecE_tot = nullptr;
TH1D* hEvtDiff_hits_minus_rec = nullptr;
TH1D* hSamplingFraction_hits = nullptr;
TH1D* hSamplingFraction_rec = nullptr;

// Hit ID histograms
TH1D* hHitSystem = nullptr;
TH1D* hHitSector = nullptr;
TH1D* hHitLayer = nullptr;
TH1D* hHitModule = nullptr;
TH1D* hHitBlock = nullptr;
TH1D* hHitFiberX = nullptr;
TH1D* hHitFiberY = nullptr;

TH1D* hRecHitSystem = nullptr;
TH1D* hRecHitSector = nullptr;
TH1D* hRecHitLayer = nullptr;
TH1D* hRecHitModule = nullptr;
TH1D* hRecHitBlock = nullptr;
TH1D* hRecHitFiberX = nullptr;
TH1D* hRecHitFiberY = nullptr;

// Acceptance
TH1D* hGenPSCAL = nullptr;
TH1D* hAccPSCAL = nullptr;
TH1D* hGenUpCAL = nullptr;
TH1D* hAccUpCAL = nullptr;
TH1D* hGenDwCAL = nullptr;
TH1D* hAccDwCAL = nullptr;

// MC particle histograms
TH1D* hMCGenStatus3_Vx = nullptr;
TH1D* hMCGenStatus3_Vy = nullptr;
TH1D* hMCGenStatus3_Vz = nullptr;
TH1D* hMCGenStatus3_Px = nullptr;
TH1D* hMCGenStatus3_Py = nullptr;
TH1D* hMCGenStatus3_Pz = nullptr;
TH1D* hMCGenStatus3_E = nullptr;

TH1D* hMCGenStatus1_Vx = nullptr;
TH1D* hMCGenStatus1_Vy = nullptr;
TH1D* hMCGenStatus1_Vz = nullptr;
TH1D* hMCGenStatus1_Px = nullptr;
TH1D* hMCGenStatus1_Py = nullptr;
TH1D* hMCGenStatus1_Pz = nullptr;
TH1D* hMCGenStatus1_E = nullptr;

TH1D* hMCGenStatus0_Vx = nullptr;
TH1D* hMCGenStatus0_Vy = nullptr;
TH1D* hMCGenStatus0_Vz = nullptr;
TH1D* hMCGenStatus0_Px = nullptr;
TH1D* hMCGenStatus0_Py = nullptr;
TH1D* hMCGenStatus0_Pz = nullptr;
TH1D* hMCGenStatus0_E = nullptr;

TH1D* hEvtPhotonE_Status3 = nullptr;
TH1D* hEvtElectronE_Status3 = nullptr;
TH1D* hEvtPositronE_Status3 = nullptr;

TH1D* hEvtPhotonE_Status1 = nullptr;
TH1D* hEvtElectronE_Status1 = nullptr;
TH1D* hEvtPositronE_Status1 = nullptr;

TH1D* hEvtPhotonE_Status0 = nullptr;
TH1D* hEvtElectronE_Status0 = nullptr;
TH1D* hEvtPositronE_Status0 = nullptr;

// Photon loss maps
TH2D* hPhotonEndX_vs_Z_all = nullptr;
TH2D* hPhotonEndY_vs_Z_all = nullptr;
TH2D* hPhotonEndPx_vs_Z_all = nullptr;
TH2D* hPhotonEndPy_vs_Z_all = nullptr;
TH2D* hPhotonEndPz_vs_Z_all = nullptr;

TH2D* hPhotonEndX_vs_Z_acc = nullptr;
TH2D* hPhotonEndY_vs_Z_acc = nullptr;
TH2D* hPhotonEndPx_vs_Z_acc = nullptr;
TH2D* hPhotonEndPy_vs_Z_acc = nullptr;
TH2D* hPhotonEndPz_vs_Z_acc = nullptr;

TH2D* hPhotonEndX_vs_Z_lost = nullptr;
TH2D* hPhotonEndY_vs_Z_lost = nullptr;
TH2D* hPhotonEndPx_vs_Z_lost = nullptr;
TH2D* hPhotonEndPy_vs_Z_lost = nullptr;
TH2D* hPhotonEndPz_vs_Z_lost = nullptr;

TH1D* hPhotonAcceptedZ = nullptr;
TH1D* hPhotonLostZ = nullptr;

TH1D* hAcceptancePSCAL = nullptr;
TH1D* hAcceptanceUpCAL = nullptr;
TH1D* hAcceptanceDwCAL = nullptr;

//=====================================
// Electron/Positron Endpoint Check
//=====================================
TH1D* hLeptonEndX = nullptr;
TH1D* hLeptonEndY = nullptr;
TH1D* hLeptonEndZ = nullptr;

TH2D* hLeptonPDG_vs_EndY = nullptr;

//=====================================
// Electron position: centroid method
//=====================================

TH1D* hDiffElectronX_Ctd[NSector] = {nullptr};
TH1D* hDiffElectronY_Ctd[NSector] = {nullptr};

TH1D* hRecElectronVx_Ctd[NSector] = {nullptr};
TH1D* hRecElectronVy_Ctd[NSector] = {nullptr};

TH1D* hGenElectronVx_Ctd[NSector] = {nullptr};
TH1D* hGenElectronVy_Ctd[NSector] = {nullptr};

//=====================================
// Electron position: slope method
//=====================================

TH1D* hDiffElectronX_Slope[NSector] = {nullptr};
TH1D* hDiffElectronY_Slope[NSector] = {nullptr};

TH1D* hRecElectronVx_Slope[NSector] = {nullptr};
TH1D* hRecElectronVy_Slope[NSector] = {nullptr};

TH1D* hGenElectronVx_Slope[NSector] = {nullptr};
TH1D* hGenElectronVy_Slope[NSector] = {nullptr};

//=====================================
// Photon position: centroid method
//=====================================

TH1D* hDiffX_Ctd = nullptr;
TH1D* hDiffY_Ctd = nullptr;

TH1D* hRecPhotonAtThinFilmVx_Ctd = nullptr;
TH1D* hRecPhotonAtThinFilmVy_Ctd = nullptr;

TH1D* hGenPhotonAtThinFilmVx_Ctd = nullptr;
TH1D* hGenPhotonAtThinFilmVy_Ctd = nullptr;

//=====================================
// Photon position: slope method
//=====================================

TH1D* hDiffX_Slope = nullptr;
TH1D* hDiffY_Slope = nullptr;

TH1D* hRecPhotonAtThinFilmVx_Slope = nullptr;
TH1D* hRecPhotonAtThinFilmVy_Slope = nullptr;

TH1D* hGenPhotonAtThinFilmVx_Slope = nullptr;
TH1D* hGenPhotonAtThinFilmVy_Slope = nullptr;

//=====================================
// Photon position: Trk method
//=====================================

TH1D* hDiffX_Trk = nullptr;
TH1D* hDiffY_Trk = nullptr;

TH1D* hRecPhotonAtThinFilmVx_Trk = nullptr;
TH1D* hRecPhotonAtThinFilmVy_Trk = nullptr;

TH1D* hGenPhotonAtThinFilmVx_Trk = nullptr;
TH1D* hGenPhotonAtThinFilmVy_Trk = nullptr;

// Reconstructed tracker hit histograms
TH1D* hTrkRecHitEnergy  = nullptr;
TH1D* hTrkRecHitTime    = nullptr;
TH1D* hTrkRecHitPrimary = nullptr;

TH1D* hTrkRecHitSystem  = nullptr;
TH1D* hTrkRecHitSector  = nullptr;
TH1D* hTrkRecHitModule  = nullptr;
TH1D* hTrkRecHitXId     = nullptr;
TH1D* hTrkRecHitYId     = nullptr;

TH1D* hTrkRecHitGlobalX = nullptr;
TH1D* hTrkRecHitGlobalY = nullptr;
TH1D* hTrkRecHitGlobalZ = nullptr;

TH2D* hTrkRecHitGlobalXY = nullptr;
TH2D* hTrkRecHitGlobalXZ = nullptr;

//============================================================
// 2D histograms: generated energy vs position quantities
//============================================================

// Centroid method
TH2D* hDiffElectronX_vs_gE_Ctd[NSector] = {nullptr};
TH2D* hDiffElectronY_vs_gE_Ctd[NSector] = {nullptr};
TH2D* hRecElectronVx_vs_gE_Ctd[NSector] = {nullptr};
TH2D* hRecElectronVy_vs_gE_Ctd[NSector] = {nullptr};
TH2D* hTrkHitx_vs_gE_Ctd[NSector]       = {nullptr};
TH2D* hTrkHity_vs_gE_Ctd[NSector]       = {nullptr};

// Slope method
TH2D* hDiffElectronX_vs_gE_Slope[NSector] = {nullptr};
TH2D* hDiffElectronY_vs_gE_Slope[NSector] = {nullptr};
TH2D* hRecElectronVx_vs_gE_Slope[NSector] = {nullptr};
TH2D* hRecElectronVy_vs_gE_Slope[NSector] = {nullptr};
TH2D* hTrkHitx_vs_gE_Slope[NSector]       = {nullptr};
TH2D* hTrkHity_vs_gE_Slope[NSector]       = {nullptr};

// Generated energy vs reconstructed energy
TH2D* hRecEnergy_vs_gE_Ctd[NSector]   = {nullptr};
TH2D* hRecEnergy_vs_gE_Slope[NSector] = {nullptr};

// Detector-wise photon maps
struct DetectorZRange {
  TString label;
  TString key;
  double zmin;
  double zmax;
  bool saveLost;
  bool saveAccepted;
};

struct PhotonDetectorHistSet {
  TH2D* hEndY_vs_EndX = nullptr;
  TH2D* hGenY_vs_GenX = nullptr;
  TH1D* hGenZ = nullptr;
  TH1D* hGenE = nullptr;
};

// Function declarations
void InitLumiSpecCALHistograms();
void SaveInFile(TFile* fout);

PhotonDetectorHistSet MakePhotonDetectorHistSet(
  const TString& prefix,
  const DetectorZRange& det,
  const TString& sample
);

std::vector<DetectorZRange> detectorRanges;
std::map<TString, PhotonDetectorHistSet> hLostDet;
std::map<TString, PhotonDetectorHistSet> hAccDet;

void AccumulateWeightedCentroid(const std::vector<float>* rec_energy,
                                       const std::vector<int>* rec_sector,
                                       const std::vector<int>* rec_layer,
                                       const std::vector<int>* rec_module,
                                       const std::vector<int>* rec_block,
                                       double w0,
                                       TVector3 posRec[NSector], 
                                       double energyRec[NSector]);

void AccumulateWeightedSlope(const std::vector<float>* rec_energy,
                                       const std::vector<int>* rec_sector,
                                       const std::vector<int>* rec_layer,
                                       const std::vector<int>* rec_module,
                                       const std::vector<int>* rec_block,
                                       const std::vector<TVector3> tracker_pos[NSector][NTracker],
                                       TVector3 posRec[NSector],
                                       double energyRec[NSector] );

#endif
