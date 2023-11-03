import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pymatgen


def plot(path, mat):

    #plot bandlines
    bands = pd.read_csv(path+"BAND.OUT", header=None, delim_whitespace=True,names=["x","y"])
    xvals = np.array(list(bands['x']))
    yvals = np.array(list(bands['y']))*27.2138
    length = np.where(xvals==0.000000)[0][1]
    for i in range(0, len(xvals), length):
        plt.plot(xvals[i:i+length], yvals[i:i+length],label='_nolegned_')
    plt.plot([0, np.max(xvals)], [0, 0], linestyle='dashed',color="black")#,label="Fermi Level")
#    if float("NaN") not in yvals:
#        plt.ylim(np.min(yvals), np.max(yvals))

    #plot significant points
    sigpoints = pd.read_csv(path+"BANDLINES.OUT", header=None, delim_whitespace=True, names=["x", "y"])
    xvals = np.array(list(sigpoints['x']))
    yvals = np.array(list(sigpoints['y']))*27.2138
    length = 2
    ticks=[]
    for i in range(0, len(xvals), length):
        plt.plot(xvals[i:i+length], yvals[i:i+length],label='_nolegned_', color='black',linestyle='dashed')
        ticks.append(xvals[i])

    sigpointlabels = []
    for point in mat["kpoints"]:
        if point[0] == "\\Gamma":
            sigpointlabels.append("\u0393")
        else:
            sigpointlabels.append(point[0])
    plt.xticks(ticks, sigpointlabels)

    plt.ylim(-5, 8)
    plt.legend()
    plt.ylabel("E-Ef (eV)")
    plt.xlabel("Wave Vector")
    plt.title(str(mat["formula_pretty"])+" Band Structure")
    plt.show()
    # plt.savefig(path+"BandPlot.png")


def plotDOS(path, mat):
    # plot bandlines
    bands = pd.read_csv(path + "PDOS_S01_A0001.OUT", header=None, delim_whitespace=True, names=["y", "x"])
    xvals = np.array(list(bands['x']))
    yvals = np.array(list(bands['y'])) * 27.2138
    length = np.where(xvals == 0.000000)[0][1]
    for i in range(0, len(xvals), length):
        plt.plot(xvals[i:i + length], yvals[i:i + length], label='_nolegned_',
                 linewidth=0.2)
    plt.plot([0, np.max(xvals)], [0, 0], linestyle='dashed', color="black")  # ,label="Fermi Level")


    plt.legend()
    plt.ylabel("E-Ef (eV)")
    plt.xlabel("Wave Vector")
    plt.title(str(mat["formula_pretty"]) + " Band Structure")
    plt.show()
    plt.savefig(path + "DOSPlot.png")