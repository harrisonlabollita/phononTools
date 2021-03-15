import sys
import os
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
    else:
        pFreq = None
        pVecs = None

    if args.irreps_file != None:
        pFreqTable = phononTools_irreps(args.irreps_file)
    else:
        pFreqTable = None

    if args.band_file != None:
        phononTools_plot(args)

    if args.band_file == None and pFreq == None and pVecs == None and pFreqTable == None:
        print("PhononTools requires at least one yaml file. For more information, try python phononTools.py -h")
    return pFreq, pVecs, pFreqTable


if __name__ == "__main__":
    pFreq, pVecs, pFreqTable = main()
