import sys
import os
try:
    import numpy as np
except:
    print("Numpy is one of my dependencies")
try:
    import yaml
except:
    print("Yaml is one of my dependencies")
try:
    import matplotlib.pyplot as plt
except:
    print("Matplotlib is one of my dependencies")
try:
    plt.style.use("publish")
except:
    print("Couldn't find the publish mpl style. This is available in the GitHub repo")

if len(sys.argv) == 1:
    print("Aborting....need at least one input file")
    sys.exit(1)

elif len(sys.argv) == 2:
    if 'mesh' in sys.argv[1]:
        mesh_file = sys.argv[1]
        irreps_file = None
        band_file = None
    elif 'irrep' in sys.argv[1]:
        mesh_file = None
        irreps_file = sys.argv[1]
        band_file = None
    else:
        mesh_file = None
        irreps_file = None
        band_file = sys.argv[1]

elif len(sys.argv) == 3:
    if 'mesh' in sys.argv[1]:
        mesh_file = sys.argv[1]
        if 'irrep' in sys.argv[2]:
            irreps_file = sys.argv[2]
            band_file = None
        else:
            band_file = sys.argv[2]
            irreps_file = None

    if 'irrep' in sys.argv[1]:
        irreps_file = sys.argv[1]
        if 'mesh' in sys.argv[2]:
            mesh_file = sys.argv[2]
            band_file = None
        else:
            band_file = sys.argv[2]
            mesh_file = None

    if 'band' in sys.argv[1]:
        band_file = sys.argv[1]
        if 'mesh' in sys.argv[2]:
            mesh_file = sys.argv[2]
            irreps_file = None
        else:
            irreps_file = sys.argv[2]
            band_file = None

else:
    if 'mesh' in sys.argv[1]:
        mesh_file = sys.argv[1]
        if 'irrep' in sys.argv[2]:
            irreps_file = sys.argv[2]
            band_file = sys.argv[3]
        else:
            band_file = sys.argv[2]
            irreps_file = sys.argv[3]

    if 'irrep' in sys.argv[1]:
        irreps_file = sys.argv[1]
        if 'mesh' in sys.argv[2]:
            mesh_file = sys.argv[2]
            band_file = sys.argv[3]
        else:
            band_file = sys.argv[2]
            mesh_file = sys.argv[3]

    if 'band' in sys.argv[1]:
        band_file = sys.argv[1]
        if 'mesh' in sys.argv[2]:
            mesh_file = sys.argv[2]
            irreps_file = sys.argv[3]
        else:
            irreps_file = sys.argv[2]
            band_file = sys.argv[3]

convert2cm1 = 33.35641

def plotParams(style):
    params = {}
    params["xlimits"]    = None
    params["ylimits"]    = None
    params["ticklabels"] = None
    params["color"]      = None
    params["save"]       = None
    with open(style) as f:
        for _, line in enumerate(f):
            params[line.split()[0][:-1]] = line.split()[1:]
    return params


def processMesh(mesh_file):
    print("Processing mesh file...{0:s}".format(mesh_file))
    mesh = yaml.load(open(mesh_file), Loader=yaml.FullLoader)
    print("Getting eigenvectors at q-point {0} for {1} atoms".format(mesh["phonon"][0]["q-position"], mesh["natom"]))
    phononFrequencies = []
    phononNormalModes = []
    for i in range(len(mesh["phonon"][0]["band"])):
        eigenValue = mesh["phonon"][0]["band"][i]["frequency"]*convert2cm1
        eigenVectors = [mesh["phonon"][0]["band"][i]["eigenvector"][atom][real][0] for atom in range(mesh["natom"]) for real in range(3)]
        phononFrequencies.append(eigenValue)
        phononNormalModes.append(eigenVectors)
    phononNormalModes = np.array(phononNormalModes).reshape(len(mesh["phonon"][0]["band"]), mesh["natom"], 3)
    return phononFrequencies, phononNormalModes

def processIrreps(irreps_file):
    print("Processing irreps file...{0:s}".format(irreps_file))
    irreps = yaml.load(open(irreps_file), Loader=yaml.FullLoader)
    print("Crystal point group: {0:s}".format(irreps["point_group"]))
    print("Found {0:d} normal modes".format(len(irreps["normal_modes"])))
    characters = [irreps["normal_modes"][ch]["ir_label"] for ch in range(len(irreps["normal_modes"]))]
    irr, counts = np.unique(characters, return_counts = True)
    string = "Irreducible representation:"
    for i in range(len(irr)):
        if i == 0:
            string += " {0:d}{1:s}".format(counts[i], irr[i])
        else:
            string += " + {0:d}{1:s}".format(counts[i], irr[i])
    print(string)
    frequency_table = {}
    for char in irr:
        frequencies = []
        for fre in range(len(irreps["normal_modes"])):
            if irreps["normal_modes"][fre]["ir_label"] == char:
                frequencies.append(irreps["normal_modes"][fre]["frequency"]*convert2cm1)
        frequency_table[char] = frequencies
    return frequency_table

def processBands(band_file):
    print("Processing band file...{0:s}".format(band_file))
    bands = yaml.load(open(band_file), Loader=yaml.FullLoader)
    qpoints = bands["nqpoint"]
    segments = bands["segment_nqpoint"]
    frequencies = []
    distance = []
    for i in range(len(bands["phonon"])):
        distance.append(bands["phonon"][i]["distance"])
    qpoint_ticks = [0, ]
    for q in segments:
        qpoint_ticks.append(q + qpoint_ticks[-1])
    qpoint_ticks[-1] -= 1
    qpoint_ticks = [distance[qpoint_ticks[q]] for q in range(len(qpoint_ticks))]
    for i in range(len(bands["phonon"])):
        for j in range(len(bands["phonon"][i]["band"])):
            frequencies.append(bands["phonon"][i]["band"][j]["frequency"]*convert2cm1)
    frequencies = np.transpose(np.array(frequencies).reshape(qpoints, len(bands["phonon"][0]["band"])))
    try:
        style = open("dispersion.style")
    except:
        plt.figure()
        myRed = '#FF503D'
        for f in range(len(frequencies)):
            plt.plot(distance, frequencies[f], color = myRed, lw = 1)
        plt.xlabel(r'q', fontsize = 15)
        plt.ylabel(r'$\omega$ (cm$^{-1}$)', fontsize = 15)
        plt.xlim(np.min(distance), np.max(distance))
        plt.ylim(np.min(frequencies) - 33, np.max(frequencies) + 33)
        plt.plot(np.linspace(np.min(distance), np.max(distance), 100), [0 for _ in range(100)], 'k--', lw = 0.5)
        plt.savefig('phonon_disperion.pdf', format = 'pdf')
    parameters = plotParams("dispersion.style")
    myRed = '#FF503D'
    plt.figure()
    plt.plot(np.linspace(np.min(distance), np.max(distance), 100), [0 for _ in range(100)], 'k--', lw = 0.5)
    if parameters["color"] != None:
        for f in range(len(frequencies)):
            plt.plot(distance, frequencies[f], color = params["color"])
    else:
        for f in range(len(frequencies)):
            plt.plot(distance, frequencies[f], color = myRed)
    if parameters["ticklabels"] != None:
        for lab in range(len(parameters["ticklabels"])):
            if parameters["ticklabels"][lab] == "Gamma":
                parameters["ticklabels"][lab] == r"$\Gamma$"
        plt.xticks(qpoint_ticks, parameters["ticklabels"])
        for q in range(len(qpoint_ticks)):
            plt.plot([qpoint_ticks[q] for _ in range(100)], np.linspace(np.min(frequencies)-100, np.max(frequencies) + 100, 100), 'k-', lw = 0.5)
    else:
        for qt in qpoint_ticks:
            plt.plot([qt for _ in range(100)], np.linspace(np.min(frequencies)-100, np.max(frequencies) + 100, 100), 'k-', lw = 0.5)

    if parameters["xlimits"] != None:
        plt.xlim(parameters["xlimits"])
    else:
        plt.xlim(np.min(distance), np.max(distance))
    if parameters["ylimits"] != None:
        plt.ylim(parameters["ylimits"])
        plt.ylabel(r'$\omega$ (cm$^{-1}$)', fontsize = 15)
    else:
        plt.ylim(np.min(frequencies) - 33, np.max(frequencies) + 33)
        plt.ylabel(r'$\omega$ (cm$^{-1}$)', fontsize = 15)
    if parameters["save"] != None:
        plt.savefig(parameters["save"]+".pdf", format = 'pdf')
    else:
        plt.savefig('phonon_disperion.pdf', format = 'pdf')

def main(mesh_file, irreps_file, band_file):
    if mesh_file != None and irreps_file != None:
        phononFrequencies, phononNormalModes = processMesh(mesh_file)
        frequency_table = processIrreps(irreps_file)
        return phononFrequencies, phononNormalModes, frequency_table
    if mesh_file != None and irreps_file == None:
        phononFrequencies, phononNormalModes = processMesh(mesh_file)
        return phononFrequencies, phononNormalModes
    if mesh_file == None and irreps_file != None:
        frequency_table = processIrreps(band_file)
        return frequency_table
    if mesh_file == None and irreps_file == None and band_file != None:
        processBands(band_file)
    if mesh_file == None and irreps_file == None and band_file == None:
        print("Aborting....")

if __name__ == "__main__":
    main(mesh_file, irreps_file, band_file)
