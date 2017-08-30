import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from collections import OrderedDict as OD
plt.rcParams['font.family']='serif'
plt.rcParams['font.weight']='light'
plt.rcParams['font.size']=14
figsize = (12,6)

def Hepevt2Pandas(filePath):
    line_l = [line[:-1].split(' ') for line in open(filePath,'r').readlines()]

    data = OD([
            ('Event',[]),
            ('nPart',[]),
            ('StatusCode',[]),
            ('Pdg',[]),
            ('FirstMother',[]),
            ('SecondMother',[]),
            ('FirstDaughter',[]),
            ('SecondDaughter',[]),
            ('pX',[]),
            ('pY',[]),
            ('pZ',[]),
            ('Energy',[]),
            ('Mass',[]),
            ('x',[]),
            ('y',[]),
            ('z',[]),
            ('Time',[]),
    ])
    event = 0
    eventLine = 0
    nParticles = 0
    Pdg = []
    FirstMother = []
    SecondMother = []
    FirstDaughter = []
    SecondDaughter = []
    pX = []
    pY = []
    pZ = []
    Energy = []
    Mass = []
    x = []
    y = []
    z = []
    Time = []

    for i,line in enumerate(line_l):
        if len(line)==2:
            eventLine = i
            event = int(line[0])
            nParticles = int(line[1])
            loopPart = 0
        if len(line)!=2:
            loopPart += 1
            Pdg.append(int(line[1]))
            FirstMother.append(int(line[2]))
            SecondMother.append(int(line[3]))
            FirstDaughter.append(int(line[4]))
            SecondDaughter.append(int(line[5]))
            pX.append(float(line[6]))
            pY.append(float(line[7]))
            pZ.append(float(line[8]))
            Energy.append(float(line[9]))
            Mass.append(float(line[10]))
            x.append(float(line[11]))
            y.append(float(line[12]))
            z.append(float(line[13]))
            Time.append(float(line[14]))
            if loopPart == nParticles:
                data['Event'].append(int(event))
                data['nPart'].append(int(nParticles))
                data['StatusCode'].append(int(line[0]))
                data['Pdg'].append(Pdg)
                data['FirstMother'].append(FirstMother)
                data['SecondMother'].append(SecondMother)
                data['FirstDaughter'].append(FirstDaughter)
                data['SecondDaughter'].append(SecondDaughter)
                data['pX'].append(pX)
                data['pY'].append(pY)
                data['pZ'].append(pZ)
                data['Energy'].append(Energy)
                data['Mass'].append(Mass)
                data['x'].append(x)
                data['y'].append(y)
                data['z'].append(z)
                data['Time'].append(Time)
                Pdg = []
                FirstMother = []
                SecondMother = []
                FirstDaughter = []
                SecondDaughter = []
                pX = []
                pY = []
                pZ = []
                Energy = []
                Mass = []
                x = []
                y = []
                z = []
                Time = []
    df = pd.DataFrame(data)
    return df

def PlotDistributions(filePath_l,var,pdg):
    for filePath in filePath_l:
        df = Hepevt2Pandas(filePath)
        if pdg!=0:
            df = df[df['Pdg']==pdg]
        fig = plt.figure(figsize=figsize)
        if var=='Time':
            gOffset = 3125
            rOffset = gOffset+1600
            n,hist,_ = plt.hist(df['Time'],bins=np.linspace(2000,10000,100))
            plt.axvline(gOffset,color='black',ls='--',lw=5)
            plt.axvline(rOffset,color='black',ls='--',lw=5)
        else:
            n,hist,_ = plt.hist(df[var],bins=100)
        plt.title('%s | %s' %(filePath.split('/')[-1],var))
        plt.xlabel(var)
        plt.grid()
        plt.show()

def Plot2dDistributions(filePath_l,var1,var2,pdg):
    for filePath in filePath_l:
        df = Hepevt2Pandas(filePath)
        if pdg!=0:
            df = df[df['Pdg']==pdg]
        fig = plt.figure(figsize=figsize)
        plt.plot(df[var1],df[var2],ls='',marker='.')
        plt.title('%s | %s vs. %s' %(filePath.split('/')[-1],var1,var2))
        plt.xlabel(var1)
        plt.ylabel(var2)
        plt.grid()
        plt.show()
