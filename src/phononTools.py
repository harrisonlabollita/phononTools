#!/usr/bin/python3
################################################################################
#                                   PhononTools                                #
# Author: Harrison LaBollita                                                   #
################################################################################
import os
import sys
import numpy as np
import yaml
import matplotlib.pyplot as plt

from parse import phononTools_parse
from mesh import phononTools_mesh
from irreps import phononTools_irreps
from plot import phononTools_plot


def main():

    args = phononTools_parse()

    if args.mesh_file != None:
        pFreq, pVecs = phononTools_mesh(args.mesh_file, args.qpoint, args.units)
        print("We suggst running an interactive python session where one types python -i")
        print("You will be able to access the frequencies and vectors by typing pFreq and pVecs")
    else:
        pFreq = None
        pVecs = None

    if args.irreps_file != None:
        pFreqTable = phononTools_irreps(args.irreps_file, args.units)
        print("We suggst running an interactive python session where one types python -i")
        print("You will be able to access the character table by typing pFreqTable")
    else:
        pFreqTable = None

    if args.band_file != None or args.dos != None:
        phononTools_plot(args)

    if args.band_file == None and pFreq == None and pVecs == None and pFreqTable == None:
        print("PhononTools requires at least one yaml file. For more information, try python phononTools.py -h")
    return pFreq, pVecs, pFreqTable

if __name__ == "__main__":
    pFreq, pVecs, pFreqTable = main()
