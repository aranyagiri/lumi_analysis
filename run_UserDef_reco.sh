#!/bin/bash

#source /opt/detector/epic-main/bin/thisepic.sh #nightly built epic main branch
source /gpfs/mnt/gpfs02/eic/aranya/epic/install/bin/thisepic.sh #andre's magnet branch

FILENUM=$1 #file number is process id, first argument

#Set number of events in each job
NEVENTS=50000

# -------- reconstruction paths --------
#subPath="brems_egen_9_to_18_gev_step_0p5_gev"
subPath="photon_pos_resol"
#subPath="cal_pos_resol"

ComDir="/eic/u/aranya/eic/data_dir/data_UserDef"
genPath="${ComDir}/genFiles/${subPath}"
recPath="${ComDir}/recFiles/${subPath}"
anaPath="${ComDir}/anaFiles/${subPath}"
#decPath="/opt/detector/epic-26.06.0/share/epic/epic_ip6.xml" #nightly built main epic branch
decPath="/gpfs/mnt/gpfs02/eic/aranya/epic/install/share/epic/epic_ip6.xml" #andre's new magnet branch
tempPath="${ComDir}/tempFiles"

mkdir -p "$genPath"
mkdir -p "$recPath"
mkdir -p "$anaPath"
mkdir -p "$tempPath"

# -------- filenames (IMPORTANT) --------

hepmcFile="${genPath}/input_${FILENUM}.hepmc" #generate events and then use.
simOut="${tempPath}/output_${FILENUM}.edm4hep.root"
recOut="${recPath}/eicrecon_${FILENUM}.root"
anaOut="${anaPath}/anaout_${FILENUM}.root"

# -------- HEPMC FILE MAKER -----
# for 0p5 step : /2
#
E=$(echo "10 + $FILENUM/1" | bc -l)

root -l -b -q "lumi_particles_Acceptance.cxx(${NEVENTS},true,true,true,${E},${E},\"${hepmcFile}\")"

#root -l -b -q "lumi_particles_Resolution.cxx(${NEVENTS},true,true,${E},${E},\"${hepmcFile}\")"

if [ ! -f "$hepmcFile" ]; then
    echo "HEPMC generation failed"
    exit 1
fi

# -------- DDSIM --------
echo "Running ddsim..."
npsim --inputFiles "$hepmcFile" \
      --outputFile "$simOut" \
      --compactFile "$decPath" \
      --numberOfEvents ${NEVENTS} \
      
if [ $? -ne 0 ]; then
    echo "ddsim failed for fileNum=$FILENUM"
    exit 1
fi

# ------- cleanup hemp file ----
echo " cleaning HEPMC File ... "
rm -f "$hepmcFile"

# -------- EICRECON --------
echo "Running eicrecon..."
eicrecon -Ppodio:output_file="$recOut" \
         -Ppodio:output_collections=MCParticles,EcalLumiSpecRawHits,EcalLumiSpecRecHits,EcalLumiSpecHits,LumiSpecTrackerRawHits,LumiSpecTrackerRecHits,LumiSpecTrackerHits \
         -Pjana:nevents=${NEVENTS} \
         -Pdd4hep:xml_files="$decPath" \
         "$simOut"

if [ $? -ne 0 ]; then
    echo "eicrecon failed for fileNum=$FILENUM"
    exit 1
fi

# ------- cleanup ddsim file ----
echo " cleaning DDSIM-EDM4HEP File ... "
rm -f "$simOut"

# --------- Analysis Code ----
root -l -q -b "LumiSpecCAL_GETaLMAnalysis.C(\"$recOut\",\"$anaOut\", false)"

if [ $? -ne 0 ]; then
    echo "Analysis failed for fileNum=$FILENUM"
    exit 1
fi

# -------- cleanup eicrecon file  --------
echo " cleaning EICRECON File ... "
rm -f "$recOut"

echo "Done for fileNum=$FILENUM"
