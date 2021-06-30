#!/usr/bin/python3
################################################################################
#                                   PhononTools                                #
# Program: dos                                                                 #
# Descripton: plot phonon density of states from phonopy                       #
# Author: Harrison LaBollita                                                   #
################################################################################
import matplotlib.pyplot as plt
import sys
import warnings
import numpy as np
import yaml
warnings.filterwarnings("ignore")
try:
    plt.style.use("publish")
except:
    print("The {} mpl-style is not currently available. This is available on Github".format("publish"))


def dos(args):
    convert2cm1 = 33.356
    convert2meV = 4.136

    colors = ['#90AF3D', '#DF9B34', '#5F82AF']

    if args.mesh_file == None:
        print("Need a mesh.yaml file to plot the phonon density of states")
        sys.exit(1)
    plt.figure()

    mesh = yaml.load(open(args.mesh_file), Loader = yaml.FullLoader)
    atom_labels = [mesh["points"][i]["symbol"] for i in range(len(mesh["points"]))]
    labels, count = np.unique(atom_labels, return_counts = True)
    print("atom labels :", labels)
    print("multiplicity :", count)
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
        if args.orientation == "v":
            plt.plot(g, omega, color = colors[atom % len(colors)], lw =1, label = labels[atom])
        else:
            plt.plot(omega, g, color = colors[atom % len(colors)], lw =1, label = labels[atom])

    if args.orientation == "v":
        plt.plot(np.linspace(0, np.max(dos), 100), [0 for _ in range(100)], 'k--', lw = 0.5)
        plt.xlim(0, gmax + 0.025*gmax)
        if args.units == "meV":
            plt.ylabel(r"$\omega$ (meV)", fontsize = 15)
            plt.xlabel(r"$D(\omega)$ (arb. units)", fontsize = 15)
            plt.ylim(np.min(omega) - 4, np.max(omega) + 4)
        elif args.units == "cm1":
            plt.ylabel(r"$\omega$ (cm$^{-1}$)", fontsize = 15)
            plt.xlabel(r"$D(\omega)$ (arb. units)", fontsize = 15)
            plt.ylim(np.min(omega) - 33, np.max(omega) + 33)
        else:
            plt.ylabel(r"$\omega$ (THz)", fontsize = 15)
            plt.xlabel(r"$D(\omega)$ (arb. units)", fontsize = 15)
            plt.ylim(np.min(omega) - 1, np.max(omega) + 1)
    else:
        plt.ylim(0, gmax + 0.025*gmax)

        plt.plot([0 for _ in range(100)], np.linspace(0, np.max(dos), 100), 'k--', lw = 0.5)
        if args.units == "meV":
            plt.xlabel(r"$\omega$ (meV)", fontsize = 15)
            plt.ylabel(r"$D(\omega)$ (arb. units)", fontsize = 15)
            plt.xlim(np.min(omega) - 4, np.max(omega) + 4)
        elif args.units == "cm1":
            plt.xlabel(r"$\omega$ (cm$^{-1}$)", fontsize = 15)
            plt.ylabel(r"$D(\omega)$ (arb. units)", fontsize = 15)
            plt.xlim(np.min(omega) - 33, np.max(omega) + 33)
        else:
            plt.xlabel(r"$\omega$ (THz)", fontsize = 15)
            plt.ylabel(r"$D(\omega)$ (arb. units)", fontsize = 15)
            plt.xlim(np.min(omega) - 1, np.max(omega) + 1)

    plt.legend()
    plt.show()





