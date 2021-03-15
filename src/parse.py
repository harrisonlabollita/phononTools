#!/usr/bin/python3
################################################################################
#                                   PhononTools                                #
# Program: parse                                                               #
# Descripton: parse the command line to drive the phononTools scripts          #
# Author: Harrison LaBollita                                                   #
################################################################################

def phononTools_parse():

    import argparse
    parser = argparse.ArgumentParser()

    # general
    parser.add_argument("-u", "--units", type=str, choices = ["meV", "cm1", "THz"], help="units to use: meV, cm^-1, or THz", default = "meV")
    parser.add_argument("-bf", "--band_file", default = None,  help="a yaml file from phonopy")
    parser.add_argument("-mf", "--mesh_file", default = None,  help="a yaml file from phonopy")
    parser.add_argument("-if", "--irreps_file", default = None, help="a yaml file from phonopy")

    # mesh
    parser.add_argument("-qp", "--qpoint", type=int, default = 0, help="index of q point")

    # plotting phonon dispersion
    parser.add_argument("-k", "--kpath", nargs="+", default = [], help="a list of kpoints for labelling x axis of phonon dispersion")
    parser.add_argument("--xmin", type=int, default = None, help = "x-axis minimum")
    parser.add_argument("--ymin", type=int, default = None, help = "y-axis minimum")
    parser.add_argument("--xmax", type=int, default = None, help = "x-axis maximum")
    parser.add_argument("--ymax", type=int, default = None, help = "y-axis maximum")
    parser.add_argument("-s", "--save", default = None, help = "file name for saving graph")
    parser.add_argument("--atoms", nargs = "+", default = [], help = "a list of atoms for which to plot the eigenvectors in the xy plane and z direction")
    parser.add_argument("--sum", default = False, help = "gather the sum of the eigenvectors for each atomic species")
    parser.add_argument("--weight", type = int, default = 100, help = "weight to scale the size of the modulus of the eigenvectors")
    parser.add_argument("--colors", nargs = "+", default = [], help = "a list of colors used for eigenvector plotting")
    # plotting dos
    parser.add_argument("--dos", type = str, default = None, help ="file for plotting density of states")
    parser.add_argument("--orientation", default = "h", choices = ["v", "h"], help = "define the axes that you would like to plot the dos")
    args = parser.parse_args()

    args.atoms = list(map(int, args.atoms))

    return args

