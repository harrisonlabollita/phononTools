import numpy as np
import matplotlib.pyplot as plt
import yaml
import sys
import warnings

warnings.filterwarnings("ignore")
try:
    plt.style.use("band_publish")
except:
    print("The {} mpl-style is not currently available. This is available on Github".format("band_publish"))


import argparse

parser = argparse.ArgumentParser()

parser.add_argument("-u", "--units", type=str, choices =["meV", "cm1"], help="units to plot dispersion: meV or cm^-1", default = "meV")
parser.add_argument("-f", "--band_file", type=str, help="a yaml file containing bands from phonopy", default = "band.yaml")
parser.add_argument("-k", "--kpath", nargs="+", default = [])
parser.add_argument("--xmin", type=int, default = None)
parser.add_argument("--ymin", type=int, default = None)
parser.add_argument("--xmax", type=int, default = None)
parser.add_argument("--ymax", type=int, default = None)
parser.add_argument("-s", "--save", default = None)
parser.add_argument("--atoms", nargs = "+", default = [], help = "a list of atoms for which to plot the eigenvectors in the xy plane and z direction")


args = parser.parse_args()

args.atoms = list(map(int, args.atoms))

def main(args):

    plt.figure()

    convert2cm1 = 33.356
    convert2meV = 4.136
    in_colors = ['r', 'b', 'g']
    out_colors = ['orange', 'cyan', 'darkgreen']
    markers = ['o', 's', '^']
    band_file = args.band_file

    print("Processing band file...{0:s}".format(band_file))

    bands = yaml.load(open(band_file), Loader = yaml.FullLoader)

    qpoints = bands["nqpoint"]
    segments = bands["segment_nqpoint"]
    atom_labels = [bands["points"][i]["symbol"] for i in range(len(bands["points"]))]

    labels, count = np.unique(atom_labels, return_counts = True)

    frequencies = []
    distance = []
    normalModes = []

    try:
        test = bands["phonon"][0]["band"][0]["eigenvector"]
    except:
        print("Eigenvectors are not found in file {0:s}".format(band_file))
        sys.exit(1)

    distance = np.array([bands["phonon"][i]["distance"] for i in range(len(bands["phonon"]))])

    qpoint_ticks = [0, ]

    for q in segments:
        qpoint_ticks.append(q + qpoint_ticks[-1])

    qpoint_ticks[-1] -= 1
    qpoint_ticks = [distance[qpoint_ticks[q]] for q in range(len(qpoint_ticks))]

    if len(args.atoms) == 0:
        for i in range(len(bands["phonon"])):
            for j in range(len(bands["phonon"][i]["band"])):
                if args.units == "meV":
                    frequencies.append(bands["phonon"][i]["band"][j]["frequency"]*convert2meV)
                elif args.units == "cm1":
                    frequencies.append(bands["phonon"][i]["band"][j]["frequency"]*convert2cm1)
        frequencies = np.transpose(np.array(frequencies).reshape(qpoints, len(bands["phonon"][0]["band"])))
        for f in range(len(frequencies)):
            plt.plot(distance, frequencies[f], color = 'k', lw =1)
    else:
        for i in range(len(bands["phonon"])):
            for j in range(len(bands["phonon"][i]["band"])):
                if args.units == "meV":
                    frequencies.append(bands["phonon"][i]["band"][j]["frequency"]*convert2meV)
                elif args.units == "cm1":
                    frequencies.append(bands["phonon"][i]["band"][j]["frequency"]*convert2cm1)
        frequencies = np.transpose(np.array(frequencies).reshape(qpoints, len(bands["phonon"][0]["band"])))
        for f in range(len(frequencies)):
            plt.plot(distance, frequencies[f], color = 'k', lw = 1)
        for a in range(len(args.atoms)):
            print("Getting the eigenvector data for atom {0:s} with index {1:d}".format(atom_labels[args.atoms[a]], args.atoms[a]))
            in_plane = []
            out_plane = []
            for i in range(len(bands["phonon"])):
                for j in range(len(bands["phonon"][i]["band"])):
                    eigenVectors = np.array([bands["phonon"][i]["band"][j]["eigenvector"][atom][real][0] for atom in range(bands["natom"]) for real in range(3)]).reshape(bands["natom"], 3)
                    in_plane.append(np.sqrt(eigenVectors[a][0]**2 + eigenVectors[a][1]**2))
                    out_plane.append(eigenVectors[a][2])
            in_plane = np.transpose(np.array(in_plane).reshape(qpoints, len(bands["phonon"][0]["band"])))
            out_plane = np.transpose(np.array(out_plane).reshape(qpoints, len(bands["phonon"][0]["band"])))
            for f in range(len(frequencies)):
                plt.scatter(distance, frequencies[f], 100*in_plane[f], color = in_colors[a % len(in_colors)], marker = markers[a % len(markers)], facecolors = "none", linewidth = 0.5)
                plt.scatter(distance, frequencies[f], 100*out_plane[f], color = out_colors[a % len(out_colors)], marker = markers[a % len(markers)], facecolors = "none", linewidth = 0.5)
            in_plane = []
            out_plane = []


    for q in range(len(qpoint_ticks)):
        plt.plot([qpoint_ticks[q] for _ in range(100)], np.linspace(np.min(frequencies)-100, np.max(frequencies) + 100, 100), 'k-', lw = 0.5)
    plt.plot(np.linspace(np.min(distance), np.max(distance), 100), [0 for _ in range(100)], 'k--', lw = 0.5)

    if args.units == "meV":
        plt.ylabel(r'$\omega$ (meV)', fontsize = 15)
    elif args.units == "cm1":
        plt.ylabel(r'$\omega$ (cm$^{-1}$)', fontsize = 15)


    if len(args.kpath) == 0:
        plt.xlabel(r"${\bf q}$", fontsize = 15)
    else:
        tick_labels = args.kpath
        for t in range(len(tick_labels)):
            if tick_labels[t] == "Gamma":
                tick_labels[t] = "$\Gamma$"
        plt.xticks(qpoint_ticks, tick_labels, fontsize = 15)

    if args.xmin == None and args.xmax == None:
        plt.xlim(np.min(distance), np.max(distance))
    elif args.xmin != None and args.xmax == None:
        plt.xlim(arg.xmin, np.max(distance))
    elif args.xmin == None and args.xmax != None:
        plt.ylim(np.min(distance), args.xmax)

    if args.ymin == None and args.ymax == None:
        if args.units == "meV":
            plt.ylim(np.min(frequencies) - 4, np.max(frequencies) +4)
        elif args.units == "cm1":
            plt.ylim(np.min(frequencies) -33, np.max(frequencies) + 33)
    else:
        if args.units == "meV":
            if args.ymin != None and args.ymax == None:
                plt.ylim(args.ymin, np.max(frequencies) + 4)
            elif args.ymin == None and args.ymax != None:
                plt.ylim(np.min(frequencies) -4 , args.ymax)
            else:
                plt.ylim(args.ymin, args.ymax)
        else:
            if args.ymin != None and args.ymax == None:
                plt.ylim(args.ymin, np.max(frequencies) + 33)
            elif args.ymin == None and args.ymax != None:
                plt.ylim(np.min(frequencies) -33 , args.ymax)
            else:
                plt.ylim(args.ymin, args.ymax)


    if args.save != None:
        plt.savefig(args.save + ".pdf", format = "pdf")
    else:
        plt.show()






if __name__ == "__main__":
    main(args)
