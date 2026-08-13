// lumi_particles_Resolution.cxx
// Build (example):
// g++ -O2 -g lumi_particles_Resolution.cxx -o lumi_particles_Resolution \
//   $(root-config --cflags --libs) -lHepMC3

#include <HepMC3/GenEvent.h>
#include <HepMC3/GenParticle.h>
#include <HepMC3/GenVertex.h>
#include <HepMC3/WriterAscii.h>
#include <HepMC3/Print.h>
#include <HepMC3/Units.h>

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TF2.h>
#include <TRandom3.h>
#include <TMath.h>
#include <TLorentzVector.h>

#include <cmath>
#include <iostream>
#include <tuple>
#include <string>
#include <memory>

using std::string;

// -----------------------------------------------------------------------------
// Globals (keep these consistent with your setup)
// -----------------------------------------------------------------------------

static constexpr double electronMass  = 0.00051099895; // GeV
static constexpr double protonMass    = 0.93827208816;  // GeV
static constexpr double photonMass    = 0.0;            // GeV
static constexpr double muonMass      = 0.1056583755;   // GeV
static constexpr double pionMass      = 0.13957039;     // GeV
static constexpr double pionZeroMass  = 0.1349768;      // GeV

// Beam momenta along z (set these correctly for your beam setup)
static double electronPz = -10.0; // GeV
static double hadronPz   = +100.0; // GeV

// Parameters used in InitializeFunctions (set to your physics model needs)
static double Z = 1.0;
static double prefactor = 1.0;

// Vectors used in InitializeFunctions
static TLorentzVector electron, hadron, electron_trf;

// Functions used in BH models
static TF1* BH_E       = nullptr;
static TF1* BH_E_trf   = nullptr;
static TF2* QED_BH     = nullptr;
static TF1* BH_Theta   = nullptr;
static TF1* BH_Phi     = nullptr;
static TF1* PDGsplitting = nullptr;

// -----------------------------------------------------------------------------
// Helpers
// -----------------------------------------------------------------------------

std::tuple<int, double> extract_particle_parameters(const std::string& particle_name) {
  if (particle_name == "electron") return std::make_tuple(11,   electronMass);
  if (particle_name == "photon")   return std::make_tuple(22,   photonMass);
  if (particle_name == "positron") return std::make_tuple(-11,  electronMass);
  if (particle_name == "proton")   return std::make_tuple(2212, protonMass);
  if (particle_name == "muon")     return std::make_tuple(13,   muonMass);
  if (particle_name == "pi0")      return std::make_tuple(111,  pionZeroMass);
  if (particle_name == "piplus")   return std::make_tuple(211,  pionMass);
  if (particle_name == "piminus")  return std::make_tuple(-211, pionMass);

  std::cerr << "Error: wrong particle name: " << particle_name << std::endl;
  std::abort();
}

void InitializeFunctions() {
  electron.SetXYZM(0.0, 0.0, electronPz, electronMass);
  hadron.SetXYZM(0.0, 0.0, hadronPz, protonMass);

  electron_trf = electron;
  electron_trf.Boost(0.0, 0.0, -hadron.Beta());

  // Bethe-Heitler photon energy in the Lab Frame
  if (BH_E) delete BH_E;
  BH_E = new TF1(
    "BH_E",
    "[0] * ([1] - x)/(x*[1]) * ([1]/([1] - x) + ([1] - x)/[1] - 2/3.) * "
    "(log(4*[2]*[1]*([1] - x)/([3]*[4]*x)) - 0.5)",
    0.1, 5.0
  );
  BH_E->SetParameter(0, Z*Z*prefactor);
  BH_E->SetParameter(1, std::fabs(electron.E()));
  BH_E->SetParameter(2, std::fabs(hadron.E()));
  BH_E->SetParameter(3, protonMass);
  BH_E->SetParameter(4, electronMass);
  BH_E->SetNpx(10000);

  // Bethe-Heitler photon energy in the hadron rest frame
  if (BH_E_trf) delete BH_E_trf;
  BH_E_trf = new TF1(
    "BH_E_trf",
    "[0] * ([1] - x)/(x*[1]) * ([1]/([1] - x) + ([1] - x)/[1] - 2/3.) * "
    "(log(2*[1]*([1] - x)/([2]*x)) - 0.5)",
    0.1, 5.0
  );
  BH_E_trf->SetParameter(0, Z*Z*prefactor);
  BH_E_trf->SetParameter(1, std::fabs(electron_trf.E()));
  BH_E_trf->SetParameter(2, electronMass);
  BH_E_trf->SetNpx(10000);

  // Lifshitz doubly-differential Bethe-Heitler in hadron rest frame
  if (QED_BH) delete QED_BH;
  QED_BH = new TF2(
    "QED_BH",
    "[0] * 1/x * ([1] - x)/[1] * y/pow(1 + y*y, 2) * "
    "(( [1]/([1] - x) + ([1] - x)/[1] - 4*y*y/pow(1 + y*y, 2) ) * "
    "log(2*[1]*([1] - x)/([2]*x)) - "
    "0.5*( [1]/([1] - x) + ([1] - x)/[1] + 2 - 16*y*y/pow(1 + y*y, 2)) )",
    0.1, 5.0, 0.0, 10.0
  );
  QED_BH->SetParameter(0, 2*Z*Z*prefactor);
  QED_BH->SetParameter(1, std::fabs(electron_trf.E()));
  QED_BH->SetParameter(2, electronMass);
  QED_BH->SetNpx(10000);
  QED_BH->SetNpy(100);

  // Bethe-Heitler photon angle wrt beam electron in the Lab Frame
  if (BH_Theta) delete BH_Theta;
  BH_Theta = new TF1("BH_Theta", "x / pow( [0]*[0] + x*x, 2)", 0.0, TMath::Pi());
  BH_Theta->SetParameter(0, electronMass / std::sqrt(electronMass*electronMass + electronPz*electronPz));
  BH_Theta->SetNpx(100000);

  // Flat azimuthal angle
  if (BH_Phi) delete BH_Phi;
  BH_Phi = new TF1("BH_Phi", "1", 0.0, 2.0*TMath::Pi());

  // Photon splitting function
  if (PDGsplitting) delete PDGsplitting;
  PDGsplitting = new TF1("PDGsplitting", "1 - 4/3.*x*(1-x)", 0.0, 1.0);
}

// -----------------------------------------------------------------------------
// Main generator
// -----------------------------------------------------------------------------

void lumi_particles_Resolution(long long n_events = 1000000,
                               bool flat = true,
                               bool convert = true,
                               double Egamma_start = 10.0,
                               double Egamma_end   = 10.0,
                               string out_fname = "genParticles.hepmc")
{
  (void)flat;
  (void)convert;
  (void)Egamma_end;

  // Output ROOT diagnostics
  TFile fout("genEventsDiagnostics.root", "RECREATE");
  TH1D BH_h1("BH_h1", "E", 100, 0, 10);
  TH1D BH_trf_h1("BH_trf_h1", "E", 100, 0, 10);
  TH2D BH_h2("BH_h2", "E vs theta", 100, 0, 10, 100, 0, 0.0035);
  TH2D QED_BH_h2("QED_BH_h2", "E vs theta", 100, 0, 10, 100, 0, 0.0035);
  TH1D QED_BH_h1("QED_BH_h1", "E", 100, 0, 10);

  // HepMC output
  HepMC3::WriterAscii hepmc_output(out_fname);

  // RNGs (stack objects, no manual delete)
  TRandom3 rtheta(0);
  TRandom3 rphi(0);
  TRandom3 rx(0);
  TRandom3 ry(0);

  InitializeFunctions();

  // Vertex origin (will be overwritten per event)
  double Vx = 0.0;
  double Vy = 0.0;
  double Vz = -62500.0; // mm

  for (long long events_parsed = 0; events_parsed < n_events; ++events_parsed) {

    // Create a fresh GenEvent every iteration (avoids internal reuse issues)
    HepMC3::GenEvent evt(HepMC3::Units::GEV, HepMC3::Units::MM);

    // One vertex
    HepMC3::GenVertexPtr v1 = std::make_shared<HepMC3::GenVertex>();

    // Beam particles (status 4)
    const double e1 = std::sqrt(electronPz*electronPz + electronMass*electronMass);
    const double e2 = std::sqrt(hadronPz*hadronPz + protonMass*protonMass);

    HepMC3::GenParticlePtr p1 = std::make_shared<HepMC3::GenParticle>(
      HepMC3::FourVector(0.0, 0.0, electronPz, e1), 11, 4
    );
    HepMC3::GenParticlePtr p2 = std::make_shared<HepMC3::GenParticle>(
      HepMC3::FourVector(0.0, 0.0, hadronPz, e2), 2212, 4
    );
    v1->add_particle_in(p1);
    v1->add_particle_in(p2);

    // Positron kinematics
    const double E_electron = Egamma_start;
    if (E_electron <= electronMass) {
      std::cerr << "Error: E_electron <= electronMass, cannot compute momentum\n";
      break;
    }
    const double p_electron = std::sqrt(E_electron*E_electron - electronMass*electronMass);

    // Your current settings (almost collinear in -z)
    double theta = 0.0;             // radians
    theta = TMath::Pi() - theta;    // convert to angle from +z axis
    double phi = 0.0;               // radians

    // If you want small angular smearing, uncomment:
    theta = TMath::Pi() - (TMath::Pi()/180.0) * rtheta.Uniform(-2.5, 2.5); //the random generation is in deg.
    phi   = (TMath::Pi()/180.0) * rphi.Uniform(-180, 180);

    //theta = TMath::Pi() - (TMath::Pi()/180.0) * 4.5; //constant high angle 
    //phi = TMath::Pi()/2.0; 90 deg along y-z plane , 0 deg along x-z plane
    //phi = TMath::Pi()/2.0 + (TMath::Pi()/180.0) * rphi.Uniform(-4.5, 4.5);

    // Vertex position smearing
    Vx = rx.Uniform(-15.0, 15.0);
    //Vy = ry.Uniform(145.0, 175.0); //top CAL
    Vy = ry.Uniform(-175.0, -145.0); //bottom CAL
    Vz = -62500.0;

    auto [id, mass] = extract_particle_parameters("positron");
    (void)mass;

    const double px = p_electron * TMath::Sin(theta) * TMath::Cos(phi);
    const double py = p_electron * TMath::Sin(theta) * TMath::Sin(phi);
    const double pz = p_electron * TMath::Cos(theta);

    HepMC3::GenParticlePtr p3 = std::make_shared<HepMC3::GenParticle>(
      HepMC3::FourVector(px, py, pz, E_electron), id, 1
    );

    v1->add_particle_out(p3);
    evt.add_vertex(v1);

    // Shift event position
    evt.shift_position_to(HepMC3::FourVector(Vx, Vy, Vz, 0.0));

    if (events_parsed == 0) {
      std::cout << "First event:\n";
      HepMC3::Print::listing(evt);
    }

    hepmc_output.write_event(evt);

    if (events_parsed % 10000 == 0) {
      std::cout << "Event: " << events_parsed << "\n";
    }

    // Fill some diagnostics if you want (examples, keep or remove)
    // BH_h1.Fill(E_electron);
    // BH_h2.Fill(E_electron, std::fabs(TMath::Pi() - theta));
  }

  hepmc_output.close();
  std::cout << "Done. Events written: " << n_events << std::endl;

  fout.cd();
  BH_h1.Write();
  BH_trf_h1.Write();
  QED_BH_h1.Write();
  BH_h2.Write();
  QED_BH_h2.Write();
  fout.Close();
}
