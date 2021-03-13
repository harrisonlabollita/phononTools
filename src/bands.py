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
parser.add_argument("-f", "--band_file", type=str, help="a yaml file containing bands from phonopy")
parser.add_argument("-k", "--kpath", nargs="+", default = [])
parser.add_argument("--xmin", type=int, default = None)
parser.add_argument("--ymin", type=int, default = None)
parser.add_argument("--xmax", type=int, default = None)
parser.add_argument("--ymax", type=int, default = None)
parser.add_argument("-s", "--save", default = None)
args = parser.parse_args()


def main(args):

    convert2cm1 = 33.356
    convert2meV = 4.136

    band_file = args.band_file

    print("Processing band file...{0:s}".format(band_file))

    bands = yaml.load(open(band_file), Loader = yaml.FullLoader)

    qpoints = bands["nqpoint"]
    segments = bands["segment_nqpoint"]

    labels, count = np.unique([bands["points"][i]["symbol"] for i in range(len(bands["points"]))], return_counts = True)

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

    La_ip = []
    La_op = []
    Ni_ip = []
    Ni_op = []
    O_ip = []
    O_op = []
    for i in range(len(bands["phonon"])):
        for j in range(len(bands["phonon"][i]["band"])):
            if args.units == "meV":
                frequencies.append(bands["phonon"][i]["band"][j]["frequency"]*convert2meV)
            elif args.units == "cm1":
                frequencies.append(bands["phonon"][i]["band"][j]["frequency"]*convert2cm1)

            eigenVectors = np.array([bands["phonon"][i]["band"][j]["eigenvector"][atom][real][0] for atom in range(bands["natom"]) for real in range(3)]).reshape(bands["natom"], 3)
            La_ip.append(sum([np.sqrt(eigenVectors[k][0]**2 + eigenVectors[k][1]**2) for k in [0]]))
            La_op.append(sum([eigenVectors[k][2] for k in [0]]))
            Ni_ip.append(sum( [np.sqrt(eigenVectors[k][0]**2 + eigenVectors[k][1]**2) for k in [4]]))
            Ni_op.append(sum([eigenVectors[k][2] for k in [4]]))
            O_ip.append(sum([np.sqrt(eigenVectors[k][0]**2 + eigenVectors[k][1]**2) for k in [9]]))
            O_op.append(sum([ eigenVectors[k][2] for k in [9]]))


    frequencies = np.transpose(np.array(frequencies).reshape(qpoints, len(bands["phonon"][0]["band"])))
    La_ip       = np.transpose(np.array(La_ip).reshape(qpoints, len(bands["phonon"][0]["band"])))
    La_op       = np.transpose(np.array(La_op).reshape(qpoints, len(bands["phonon"][0]["band"])))
    Ni_ip       = np.transpose(np.array(Ni_ip).reshape(qpoints, len(bands["phonon"][0]["band"])))
    Ni_op       = np.transpose(np.array(Ni_op).reshape(qpoints, len(bands["phonon"][0]["band"])))
    O_ip       = np.transpose(np.array(O_ip).reshape(qpoints, len(bands["phonon"][0]["band"])))
    O_op       = np.transpose(np.array(O_op).reshape(qpoints, len(bands["phonon"][0]["band"])))


    plt.figure()

    for f in range(len(frequencies)):
        plt.plot(distance, frequencies[f], color = 'k', lw = 1)
        plt.scatter(np.delete(distance, np.arange(0, distance.size, 2)), np.delete(frequencies[f], np.arange(0, frequencies[f].size, 2)), 100*np.delete(La_ip[f], np.arange(0, La_ip[f].size, 2)), marker = "o", color = 'red', linewidth = 0.5)
        plt.scatter(np.delete(distance, np.arange(0, distance.size, 2)), np.delete(frequencies[f], np.arange(0, frequencies[f].size, 2)), 100*np.delete(La_op[f], np.arange(0, La_op[f].size, 2)), marker = "o", color = 'orange', linewidth = 0.5)
        plt.scatter(np.delete(distance, np.arange(0, distance.size, 2)), np.delete(frequencies[f], np.arange(0, frequencies[f].size, 2)), 100*np.delete(Ni_ip[f], np.arange(0, Ni_ip[f].size, 2)), marker = "^", color = 'blue',  linewidth = 0.5)
        plt.scatter(np.delete(distance, np.arange(0, distance.size, 2)), np.delete(frequencies[f], np.arange(0, frequencies[f].size, 2)), 100*np.delete(Ni_op[f], np.arange(0, Ni_op[f].size, 2)), marker = "^", color = 'cyan',  linewidth = 0.5)
        plt.scatter(np.delete(distance, np.arange(0, distance.size, 2)), np.delete(frequencies[f], np.arange(0, frequencies[f].size, 2)), 100*np.delete(O_ip[f], np.arange(0, O_ip[f].size, 2)), marker = "s", color = 'green',  linewidth = 0.5)
        plt.scatter(np.delete(distance, np.arange(0, distance.size, 2)), np.delete(frequencies[f], np.arange(0, frequencies[f].size, 2)), 100*np.delete(O_op[f], np.arange(0, O_op[f].size, 2)), marker = "s", color = 'limegreen', linewidth = 0.5)

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
