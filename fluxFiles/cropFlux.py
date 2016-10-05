import numpy as np
import os, sys, argparse

def getScriptPath():
    return os.path.dirname(os.path.realpath(sys.argv[0]))

def main():
    # Parser snippet
    parser = argparse.ArgumentParser()
    parser.add_argument("energy", help="Energy (in GeV) value at which crop is applied",type=float)
    args = parser.parse_args()

    # File path
    fluxDir = '/Users/sdporzio/MyPythonPackages/SterileNeutrinosUtils/BnbNuFlux/'
    fluxFile = 'flux_numu.dat'
    outFile = getScriptPath()+'/'+fluxFile.replace('.dat','_%.2f.dat' %(args.energy))
    print outFile

    pathIn = fluxDir + fluxFile
    dataIn = np.genfromtxt(pathIn,delimiter=' ',skip_header=1,names=True,dtype=None)

    # Crop data
    for i,data in enumerate(dataIn):
        dataIn['Energy'][i] = dataIn['Energy'][i]/1000.
        if data['Energy']<=args.energy:
            dataIn['Flux'][i] = 0

    # Save output
    np.savetxt(outFile,dataIn,delimiter=' ',fmt='%s')

if __name__ == "__main__":
    main()
