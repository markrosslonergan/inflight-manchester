run(){
  COLOR='\033[1;33m'
  DEFAULT='\033[0m'
  echo -e "${COLOR}-> ${1}${DEFAULT}";
  eval ${1};
}

#  Settings
MASS="0.25" #in GeV
NUMBER="20000"
CHANNEL="2"
FLUXFILE="FluxFiles/flux_numu_cut${MASS}.dat"
OUTFILE="Output/sterileEvents_mass${MASS}_channel${CHANNEL}_n${NUMBER}.dat"

IF_CUTS="--using-no-detector"
IF_HEPEVT="--hepevt"

# Commands
run "python FluxFiles/generateSterileFlux.py ${MASS}"
run "./inflight --mass ${MASS} --number ${NUMBER} --channel ${CHANNEL} --flux-file ${FLUXFILE} ${IF_CUTS} ${IF_HEPEVT} > ${OUTFILE}"
# run "./inflight --mass ${MASS} --number ${NUMBER} --channel ${CHANNEL} --flux-file ${FLUXFILE} ${IF_CUTS} ${IF_HEPEVT}"


# ******************************************
# Allowed arguments:
# 	-m, --mass			Sets the parent sterile mass. [default = 0.1500]
#
# 	-n, --number			How many events will we generate? [Default 100]
#
# 	-C, --channel			sets the decay channel [default = 1]
#
# 						0: 3-body (nu e e)
#
# 						1: Isotropic 2-body e pi
#
# 						2: Isotropic 2-body mu pi
#
# 						3: Isotropic 2-body nu Pi0
#
# 						4: Isotropic 2-body nu gamma
#
# 						5: 3-body (nu mu mu)
#
# 						6: 3-body (nu mu e)
#
# 	--using-no-detector		runs with no detector cuts. [DEFAULT]
#
# 	--using-muboone			runs with muBooNE detector cuts.
#
# 	-f, --flux-file			The name of the file containing flux to sample from, no normalisation needed. 1st column energy in GeV, 2nd column flux.
#
# 	-v, --hepevt			Outputs to hepevt format.
#
# 	-h, --help			Displays this message
#
# ******************************************
# Example usage, if file flux_150.dat exists in the directory then to generate 2000 N->mu-pi+ decays in the hepevt format either
#
#   ./inflight -m 0.150 -n 2000 -C 2 -v -f flux_150.dat
#
#   ./inflight --mass 0.150 --number 2000 --channel 2 --hepevt --flux-file flux_150.dat
