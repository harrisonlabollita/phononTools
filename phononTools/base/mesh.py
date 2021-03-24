#!/usr/bin/python3
################################################################################
#                                   PhononTools                                #
# Program: mesh                                                                #
# Description: Gather the phonon frequencies and eigenvectors at a given       #
# qpoint                                                                       #
# Author: Harrison LaBollita                                                   #
################################################################################
import yaml
import numpy as np
def mesh(mesh_file, qpoint, unit):
    # Description: Routine to process a mesh.yaml file
    # Inputs: mesh.yaml file, qpoint to get eigenvectors and eigenvalues, unit
    # to convert phonon frequency to.
    # Output: frequencies and eigenvectors at the given qpoint
    convert2cm1 = 33.356
    convert2meV = 4.136

    print("Processing mesh file...{0:s}".format(mesh_file))
    mesh = yaml.load(open(mesh_file), Loader=yaml.FullLoader)
    print("Getting eigenvectors at q-point {0} for {1} atoms".format(mesh["phonon"][qpoint]["q-position"], mesh["natom"]))
    phononFrequencies = []
    phononNormalModes = []
    for i in range(len(mesh["phonon"][qpoint]["band"])):
        if unit == "meV":
            eigenValue = mesh["phonon"][qpoint]["band"][i]["frequency"]*convert2meV
        elif unit == "cm1":
            eigenValue = mesh["phonon"][qpoint]["band"][i]["frequency"]*convert2cm1
        elif unit == "THz":
            eigenValue = mesh["phonon"][qpoint]["band"][i]["frequency"]
        else:
            eigenValue = mesh["phonon"][qpoint]["band"][i]["frequency"]*convert2meV
        eigenVectors = [mesh["phonon"][qpoint]["band"][i]["eigenvector"][atom][real][0] for atom in range(mesh["natom"]) for real in range(3)]
        phononFrequencies.append(eigenValue)
        phononNormalModes.append(eigenVectors)
    phononNormalModes = np.array(phononNormalModes).reshape(len(mesh["phonon"][qpoint]["band"]), mesh["natom"], 3)
    return phononFrequencies, phononNormalModes
