# EIC_Root_F4A_Merrick_Lavinsky

# This document is for the EIC R&D efforts to determine the best setup for the new detector. These simulations are of the EIC with combined silicon and gaseous detectors. The inal result will be a combination of Alexander Kiselev's repository (https://github.com/eic/EicToyModel.git) and Hakann Wennlof's repository (https://gitlab.com/hwennlof/fun4allgdmlimport.git). Shown below are the commands I use to start each simulation in Fun4All.


# 1) Using Singularity and cvmfs, I visualized Alexanders Fun4All_G4_Tracking.C macro via the following commands:

#singularity shell -B /cvmfs:/cvmfs -B /home/merrick/EicToyModel:/scratch /cvmfs/eic.opensciencegrid.org/singularity/rhic_sl7_ext.simg

#source /cvmfs/eic.opensciencegrid.org/x8664_sl7/opt/fun4all/core/bin/eic_setup.sh -n

#cd /scratch && mkdir -p /scratch/build

#cd /scratch/build

#cmake -DCMAKE_INSTALL_PREFIX=. -DGEANT=YES -Wno-dev -DVGM=${OPT_SPHENIX}/vgm ..

#make -j6 install

#export LD_LIBRARY_PATH=/scratch/build/lib:${OPT_SPHENIX}/vgm/lib64:${LD_LIBRARY_PATH}

#root -l ../scripts/eicroot.C
#.q

#cd /scratch/fun4all_with_eicroot && mkdir -p build && cd build

#../sandbox/autogen.sh --prefix=/scratch/fun4all_with_eicroot

#make -j6 install

#source /cvmfs/eic.opensciencegrid.org/x8664_sl7/opt/fun4all/core/bin/setup_local.sh /scratch/fun4all_with_eicroot

#export ROOT_INCLUDE_PATH=/scratch/build/include/etm:${ROOT_INCLUDE_PATH}

#cd ../macro

#root -l Fun4All_G4_Tracking.C
#.q

