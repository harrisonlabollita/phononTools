#!/usr/bin/python3
################################################################################
#                                   PhononTools                                #
# Program: band_dos                                                            #
# Descripton: plot phonon dispersion and density of states                     #
# Author: Harrison LaBollita                                                   #
################################################################################
import matplotlib.pyplot as plt
import yaml
import sys
import warnings
import numpy as np
from phononTools.base.units import convert2cm1, convert2meV


def unique(atom_labels):
    uniq, index, counts = np.unique(atom_labels, return_index=True, return_counts=True)
    return uniq[index.argsort()], counts[index.argsort()]
            

def band_dos(args):
    myRed = '#FF503D'
    myGreen = '#90AF3D'
    myOrange = '#DF9B34'
    myBlue   = '#5F82AF'

    colors = [myOrange, myBlue, myGreen]
    fig, ax = plt.subplots(1, 2, sharey=True, figsize = (7,5))

    bands = yaml.load(open(args.band_file), Loader = yaml.FullLoader)

    qpoints = bands["nqpoint"]
    segments = bands["segment_nqpoint"]
    atom_labels = [bands["points"][i]["symbol"] for i in range(len(bands["points"]))]
    labels, count = unique(atom_labels)

    print("atoms :", labels)
    print("mult. :", count)

    frequencies = []
    distance = []

    distance = np.array([bands["phonon"][i]["distance"] for i in range(len(bands["phonon"]))])

    qpoint_ticks = [0, ]
    for q in segments:
        qpoint_ticks.append(q + qpoint_ticks[-1])

    qpoint_ticks[-1] -= 1
    qpoint_ticks = [distance[qpoint_ticks[q]] for q in range(len(qpoint_ticks))]

    for i in range(len(bands["phonon"])):
        for j in range(len(bands["phonon"][i]["band"])):
            if args.units == "meV":
                frequencies.append(bands["phonon"][i]["band"][j]["frequency"]*convert2meV)
            elif args.units == "cm1":
                frequencies.append(bands["phonon"][i]["band"][j]["frequency"]*convert2cm1)
            else:
                frequencies.append(bands["phonon"][i]["band"][j]["frequency"])
    frequencies = np.transpose(np.array(frequencies).reshape(qpoints, len(bands["phonon"][0]["band"])))
    for f in range(len(frequencies)):
        ax[0].plot(distance, frequencies[f], color = myRed, lw = 1)

    for q in range(len(qpoint_ticks)):
        ax[0].plot([qpoint_ticks[q] for _ in range(100)], np.linspace(np.min(frequencies)-100, np.max(frequencies) + 100, 100), 'k-', lw = 0.5)

    ax[0].plot(np.linspace(np.min(distance), np.max(distance), 100), [0 for _ in range(100)], 'k--', lw = 0.5)

    if args.units == "meV":
        ax[0].set_ylabel(r'$\omega$ (meV)', fontsize = 15)
    elif args.units == "cm1":
        ax[0].set_ylabel(r'$\omega$ (cm$^{-1}$)', fontsize = 15)


    if len(args.kpath) == 0:
        ax[0].set_xticks(qpoint_ticks)
        ax[0].set_xlabel(r"${\bf q}$", fontsize = 15)
    else:
        tick_labels = args.kpath
        for t in range(len(tick_labels)):
            if tick_labels[t] == "Gamma":
                tick_labels[t] = "$\Gamma$"
        ax[0].set_xticks(qpoint_ticks)
        ax[0].set_xticklabels(tick_labels)

    if args.xmin == None and args.xmax == None:
        ax[0].set_xlim(np.min(distance), np.max(distance))
    elif ax[0].xmin != None and args.xmax == None:
        ax[0].set_xlim(arg.xmin, np.max(distance))
    elif args.xmin == None and args.xmax != None:
        ax[0].set_ylim(np.min(distance), args.xmax)

    if args.ymin == None and args.ymax == None:
        if args.units == "meV":
            ax[0].set_ylim(np.min(frequencies) - 4, np.max(frequencies) +4)
        elif args.units == "cm1":
            ax[0].set_ylim(np.min(frequencies) -33, np.max(frequencies) + 33)
    else:
        if args.units == "meV":
            if args.ymin != None and args.ymax == None:
                ax[0].set_ylim(args.ymin, np.max(frequencies) + 4)
            elif args.ymin == None and args.ymax != None:
                ax[0].set_ylim(np.min(frequencies) -4 , args.ymax)
            else:
                ax[0].set_ylim(args.ymin, args.ymax)
        else:
            if args.ymin != None and args.ymax == None:
                ax[0].set_ylim(args.ymin, np.max(frequencies) + 33)
            elif args.ymin == None and args.ymax != None:
                ax[0].set_ylim(np.min(frequencies) -33 , args.ymax)
            else:
                ax[0].set_ylim(args.ymin, args.ymax)

    dos = np.loadtxt(args.dos)

    if args.units == "meV":
        omega = [dos[i][0]*convert2meV for i in range(len(dos))]
    elif args.units == "cm1":
        omega = [dos[i][0]*convert2cm1 for i in range(len(dos))]
    else:
        omega = [dos[i][0] for i in range(len(dos))]

    start = 1
    gmax = 0
    for atom in range(len(labels)):
        g = []
        stop = start + count[atom]
        for i in range(len(dos)):
            g.append(sum([dos[i][a] for a in range(start, stop)])/count[atom])
        gmax = np.max(g) if np.max(g) > gmax else gmax
        start = stop
        ax[1].plot(g, omega, color = colors[atom % len(colors)], lw =1, label = labels[atom])

    ax[1].plot(np.linspace(0, np.max(dos), 100), [0 for _ in range(100)], 'k--', lw = 0.5)
    ax[1].set_xlim(0, gmax + 0.025*gmax)
    ax[1].legend()
    if args.units == "meV":
        ax[1].set_xlabel(r"$D(\omega)$ (arb. units)", fontsize = 15)
    elif args.units == "cm1":
        ax[1].set_xlabel(r"$D(\omega)$ (arb. units)", fontsize = 15)
    else:
        ax[1].xlabel(r"$D(\omega)$ (arb. units)", fontsize = 15)
    if args.save != None:
        print("Saving figure to file {}".format(args.save +".pdf"))
        plt.savefig(args.save + ".pdf", format = "pdf")
    else:
        plt.show()

