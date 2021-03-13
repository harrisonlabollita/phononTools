import numpy as np
import matplotlib.pyplot as plt
import yaml
import sys

try:
    plt.style.use("band_publish")
except:
    print("The {} mpl-style is not currently available. This is available on Github".format("band_publish"))

convert2cm1 = 33.35641
convert2meV = 4.132231405

band_file = sys.argv[1]

print("Processing band file...{0:s}".format(band_file))

bands = yaml.load(open(band_file), Loader = yaml.FullLoader)

qpoints = bands["nqpoint"]
segments = bands["segment_nqpoint"]

labels, count = np.unique([bands["points"][i]["symbol"] for i in range(len(bands["points"]))], return_counts = True)


frequencies = []
distance = []
normalModes = []

try:
    test = bands["phonon"][0]["band"][0]["frequency"]
except:
    print("Eigenvectors are not found in file {0:s}".format(band_file))

for i in range(len(bands["phonon"])):
    distance.append(bands["phonon"][i]["distance"])
distance = np.array(distance)
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
        frequencies.append(bands["phonon"][i]["band"][j]["frequency"]*convert2meV)
        eigenVectors = np.array([bands["phonon"][i]["band"][j]["eigenvector"][atom][real][0] for atom in range(bands["natom"]) for real in range(3)]).reshape(bands["natom"], 3)
        La_ip.append(sum([np.sqrt(eigenVectors[k][0]**2 + eigenVectors[k][1]**2) for k in range(0, 4) ])/count[0])
        La_op.append(sum([eigenVectors[k][2] for k in range(0, 4)])/count[0])
        Ni_ip.append(sum( [np.sqrt(eigenVectors[k][0]**2 + eigenVectors[k][1]**2) for k in range(4, 7)])/count[1])
        Ni_op.append(sum([eigenVectors[k][2] for k in range(4, 7)])/count[1])
        O_ip.append(sum([np.sqrt(eigenVectors[k][0]**2 + eigenVectors[k][1]**2) for k in range(7, 15)])/count[2])
        O_op.append(sum([ eigenVectors[k][2] for k in range(7, 15)])/count[2])


frequencies = np.transpose(np.array(frequencies).reshape(qpoints, len(bands["phonon"][0]["band"])))
La_ip       = np.transpose(np.array(La_ip).reshape(qpoints, len(bands["phonon"][0]["band"])))
La_op       = np.transpose(np.array(La_op).reshape(qpoints, len(bands["phonon"][0]["band"])))
Ni_ip       = np.transpose(np.array(Ni_ip).reshape(qpoints, len(bands["phonon"][0]["band"])))
Ni_op       = np.transpose(np.array(Ni_op).reshape(qpoints, len(bands["phonon"][0]["band"])))
O_ip       = np.transpose(np.array(O_ip).reshape(qpoints, len(bands["phonon"][0]["band"])))
O_op       = np.transpose(np.array(O_op).reshape(qpoints, len(bands["phonon"][0]["band"])))

tick_labels = ["$\Gamma$", "M", "N", "P", "X"]

plt.figure()

for f in range(len(frequencies)):
    plt.plot(distance, frequencies[f], color = 'k', lw = 1)
    plt.scatter(np.delete(distance, np.arange(0, distance.size, 2)), np.delete(frequencies[f], np.arange(0, frequencies[f].size, 2)), 100*np.delete(La_ip[f], np.arange(0, La_ip[f].size, 2)), marker = "o", color = 'red', facecolors = "none", linewidth = 0.5)
    plt.scatter(np.delete(distance, np.arange(0, distance.size, 2)), np.delete(frequencies[f], np.arange(0, frequencies[f].size, 2)), 100*np.delete(La_op[f], np.arange(0, La_op[f].size, 2)), marker = "o", color = 'orange', facecolors = "none", linewidth = 0.5)
    plt.scatter(np.delete(distance, np.arange(0, distance.size, 2)), np.delete(frequencies[f], np.arange(0, frequencies[f].size, 2)), 100*np.delete(Ni_ip[f], np.arange(0, Ni_ip[f].size, 2)), marker = "^", color = 'blue', facecolors = "none", linewidth = 0.5)
    plt.scatter(np.delete(distance, np.arange(0, distance.size, 2)), np.delete(frequencies[f], np.arange(0, frequencies[f].size, 2)), 100*np.delete(Ni_op[f], np.arange(0, Ni_op[f].size, 2)), marker = "^", color = 'cyan', facecolors = "none", linewidth = 0.5)
    plt.scatter(np.delete(distance, np.arange(0, distance.size, 2)), np.delete(frequencies[f], np.arange(0, frequencies[f].size, 2)), 100*np.delete(O_ip[f], np.arange(0, O_ip[f].size, 2)), marker = "s", color = 'green', facecolors = "none", linewidth = 0.5)
    plt.scatter(np.delete(distance, np.arange(0, distance.size, 2)), np.delete(frequencies[f], np.arange(0, frequencies[f].size, 2)), 100*np.delete(O_op[f], np.arange(0, O_op[f].size, 2)), marker = "s", color = 'limegreen', facecolors = "none", linewidth = 0.5)

for q in range(len(qpoint_ticks)):
    plt.plot([qpoint_ticks[q] for _ in range(100)], np.linspace(np.min(frequencies)-100, np.max(frequencies) + 100, 100), 'k-', lw = 0.5)

#plt.ylabel(r'$\omega$ (cm$^{-1}$)', fontsize = 15)
plt.ylabel(r'$\omega$ (meV)', fontsize = 15)
plt.xticks(qpoint_ticks, tick_labels, fontsize = 15)
plt.xlim(np.min(distance), np.max(distance))
plt.ylim(np.min(frequencies) - 4, np.max(frequencies) +4)
plt.show()
