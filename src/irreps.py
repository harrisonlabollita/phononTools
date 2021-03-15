################################################################################
#                                   PhononTools                                #
# Program: irrep                                                               #
# Descripton: generate a character table with corresponding frequencies        #
# Author: Harrison LaBollita                                                   #
################################################################################
import yaml

def phononTools_irreps(irreps_file):
    # Description: Process a irreps.yaml file
    # Input: yaml file from phonopy
    # Output: A character table with the corresponding frequncies

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
