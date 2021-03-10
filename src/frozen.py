# setup a frozen phonon calculation for a cluster
# after running phonopy -d --dim=" X X X" -c POSCAR
# this script will setup the directories for the cluster
import sys
import os

try:
    maxDisplacement = int(sys.argv[1])
except:
    print("Need the number of displacement files given as an argument")
    sys.exit(1)

try:
    slurm = sys.argv[2]
except:
    print("Need the cluster submission script name")
    sys.exit(1)

for i in range(1, maxDisplacement + 1):
    if i < 10:
        os.mkdir("disp-00{}".format(i))
        os.rename("POSCAR-00{}".format(i), "disp-00{}/POSCAR".format(i))
        os.system("cp INCAR POTCAR KPOINTS {} disp-00{}".format(slurm, i))
    if i >= 10 and i < 100:
        os.mkdir("disp-0{}".format(i))
        os.rename("POSCAR-0{}".format(i), "disp-0{}/POSCAR".format(i))
        os.system("cp INCAR POTCAR KPOINTS {} disp-0{}".format(slurm, i))
    if i >= 100:
        os.mkdir("disp-{}".format(i))
        os.rename("POSCAR-{}".format(i), "disp-{}/POSCAR".format(i))
        os.system("cp INCAR POTCAR KPOINTS {} disp-{}".format(slurm, i))
