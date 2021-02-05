#!/bin/bash
#$ -cwd
#$ -j y
#$ -N twoD_singleP_1
#$ -o $JOB_NAME.$JOB_ID.output
#$ -M amir.golparvar@ufz.de
#$ -m beas
#$ -l h_rt=12:0:00
#$ -l h_vmem=3G
#$ -pe openmpi-orte 6
##$ -l cpu_model=E5-2690v4
module load foss/2018b
#module load OpenFOAM/6
module load foss/2018a OpenFOAM-Extend/4.0
echo "Executing job commands, current working directory is `pwd`"
source $FOAM_BASH
source $HOME/works/GeoChemFoam-4.1/etc/bashrc
#source /home/golparva/works/poreFoam-singlePhase/src/bashrc
#export PATH=$PATH:$HOME/works/poreFoam-singlePhase/src/scripts/singlePhase
#export nProcX=1
#export nProcY=1
#export nProcZ=1

#AllRunImagePar HalleOrg.mhd
#cp system/decomposeParDict.16 system/decomposeParDict
#decomposePar
mpirun -np 6 bioGrowthDOM -parallel
#mpiexec -np 4 interTransportFoam -parallel
#mpiexec -np 16 interFoam -parallel
