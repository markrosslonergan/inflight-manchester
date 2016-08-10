# inflight-manchester

Inflight, an event generator for heavy sterile decay spectrums at SBL facilities
If you have any questions please contact the authors mark.ross-lonergan@durham.ac.uk or peter.ballett@durham.ac.uk
The authors make no guarrentee of this codes stability or bug-free-ness. Use is at own risk.
******************************************
Allowed arguments:
	-m, --mass			Sets the parent sterile mass. [default = 0.1500]

	-n, --number			How many events will we generate? [Default 100]

	-C, --channel			sets the decay channel [default = 1] 

						0: 3-body (nu e e).

						1: Isotroic 2-body e pi

						2: Isotropic 2-body mu pi

						3: Isotropic 2-body nu Pi0

						4: Isotropic 2-body nu gamma

	--using-no-detector		runs with no detector cuts. [DEFAULT] 

	--using-muboone			runs with muBooNE detector cuts.

	-f, --flux-file			The name of the file containing flux to sample from, no normalisation needed. 1st column energy in GeV, 2nd column flux.

	-v, --hepevt			Outputs to hepevt format.

	-h, --help			Displays this message

******************************************
Example usage, if file flux_150.dat exists in the directory then to generate 2000 N->mu-pi+ decays in the hepevt format either

  ./inflight -m 0.150 -n 2000 -C 2 -v -f flux_150.dat 

  ./inflight --mass 0.150 --number 2000 --channel 2 --hepevt --flux-file flux_150.dat


