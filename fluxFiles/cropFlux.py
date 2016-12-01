import numpy as np
import os, sys, argparse

def getScriptPath():
    return os.path.dirname(os.path.realpath(sys.argv[0]))+'/'

def main():
    # Parser snippet
    parser = argparse.ArgumentParser()
    parser.add_argument("energy", help="Energy (in GeV) value at which crop is applied",type=float)
    args = parser.parse_args()

    # File path
    inFile = getScriptPath()+'OriginalFlux/flux_numu.dat'
    outFile = getScriptPath()+'flux_numu_cut%.2f.dat' %(args.energy)
    inData = np.genfromtxt(inFile,delimiter=' ',skip_header=1,names=True,dtype=None)

    # Crop data
    for i,data in enumerate(inData):
        inData['Energy'][i] = inData['Energy'][i]/1000.
        if data['Energy']<=args.energy:
            inData['Flux'][i] = 0

    # Save output
    np.savetxt(outFile,inData,delimiter=' ',fmt='%s')

if __name__ == "__main__":
    main()
