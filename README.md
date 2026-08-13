Both are different method, first one is used for personalised EICrecon analysis code, second one uses general EICrecon output.

1. To use analyzeLumiHits

cmake -S analyzeLumiHits -B analyzeLumiHits/build

cmake --build analyzeLumiHits/build --target install 

NOTE : The second command should be executed after each modifications, the first after significant modifications or some errors.

2. To perform lumi analysis in BNL-SDCC

a. lumi_particles_Acceptance.cxx - to generate particle's HepMC files for full simulation.

b. run_UserDef_reco.sh performs the DD4hep simulation. NOTE : Make sure to edit the source file, detector path, Acceptance Hepmc generator/Resolution Hepmc generator accordingly.

c. LumiSpecCAL_GETaLMAnalysis.C and LumiSpecCAL.h - analysis code 

  1. CAL energy and position resolution, sampling fraction.
  2. photon's energy and position resolution.
  3. Acceptance of PSCAL's
  4. Tracker's hits are used as of now, need a general EICrecon update to get rec. hits.
     
d. condor_submit run_UserDef.job . Note : Make sure to edit the correct ____Acceptance.cxx / _____Resolution.cxx line.

e. lumi_particles_Resolution.cxx - to generate particle's HepMC files for CAL standalone simulation.


