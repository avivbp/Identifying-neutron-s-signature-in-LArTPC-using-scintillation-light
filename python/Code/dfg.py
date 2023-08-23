import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LogNorm

from functions import *

dat = pd.read_csv("../phot_eKin.csv")
energies = dat.electron_eKin
plt.hist(energies)
plt.show()

dat = pd.read_csv("../numPhotSec_3keV.csv")
num = dat.numPhot
print(max(num))
h, bins = np.histogram(num, np.linspace(0, 30, 15))
bins_centers = (bins[1:] + bins[:-1]) / 2
plt.step(bins_centers, h, label='3 keV secondary electron')
plt.xlabel('number of photons')
plt.ylabel('frequency')
plt.title('number of optical photons produced by secondary 3 keV electrons')
plt.show()

dat = pd.read_csv("../Scintillation_Data/scint_yield_checks/scintDist_56keV_electron_0.89.csv")
num = dat.numPhotons
h, bins = np.histogram(num, np.linspace(0, 2000, 50), density=True)
bins_centers = (bins[1:] + bins[:-1]) / 2
plt.step(bins_centers, h, label='56 keV primary electron')

dat = pd.read_csv("../Scintillation_Data/scint_yield_checks/numPhotSec_56keV.csv")
num = dat.numPhot[dat.numPhot > 0]
# print(len(num[num == 0]))
# print(len(num))
h, bins = np.histogram(num, np.linspace(0, 2000, 50), density=True)
bins_centers = (bins[1:] + bins[:-1]) / 2
plt.step(bins_centers, h, label='56 keV secondary electron')

dat = pd.read_csv("../Scintillation_Data/scint_yield_checks/scintDist_56keV_gamma_0.89.csv")
num = dat.numPhotons
h, bins = np.histogram(num, np.linspace(0, 2000, 50), density=True)
bins_centers = (bins[1:] + bins[:-1]) / 2
plt.step(bins_centers, h, label='56 keV $\gamma$')

plt.legend()
plt.xlabel('number of photons')
plt.ylabel('frequency')
plt.title('number of optical photons produced by primary vs secondary 56 keV electrons vs 56 keV $\gamma$')
plt.show()

dat = pd.read_csv("../scintDist.csv")
num_photons = dat.numPhotons
num_photons = num_photons[num_photons > 0]
h1, bins1 = np.histogram(num_photons, 100, density=True)
bins_centers = (bins1[1:] + bins1[:-1]) / 2
plt.step(bins_centers, h1, label='60 keV $\gamma$ , infinite detector')
ymax1 = 1.1 * np.max(h1)
xmax1 = bins1[np.where(h1 == np.max(h1))[0][0]]
plt.plot([xmax1, xmax1], [0, ymax1], '--', color='blue')

dat = pd.read_csv("../scintDist_e.csv")
num_photons = dat.numPhotons
num_photons = num_photons[num_photons > 0]
h, bins = np.histogram(num_photons, 100, density=True)
bins_centers = (bins[1:] + bins[:-1]) / 2
plt.step(bins_centers, h, label='60 keV electron , infinite detector')
ymax = 1.1 * np.max(h)
xmax = bins[np.where(h == np.max(h))[0][0]]
plt.plot([xmax, xmax], [0, ymax], '--', color='orange')
xticks = list(np.linspace(0, 2000, 11))
xticks.extend([xmax, xmax1])
xticks.remove(2000)
xticks.append(2150)
plt.xticks(xticks)
plt.ylim([0, 0.0031])
plt.legend()
plt.show()

# ymax = 1.1 * np.max(h)
# plt.plot([2731, 2731], [0, ymax], '--r')
# plt.xlabel('Number of generated photons')
# plt.ylabel('Frequency [a.u.]')
# plt.title('Number of generated photons from 60 keV $\gamma$ ray in LAr')
# plt.show()
#
# h1, bins1 = np.histogram(num_photons, 100)
# bins1_centers = (bins1[1:] + bins1[:-1]) / 2
#
# print(bins1_centers[65])
# print(bins1_centers[90])
#
# parameters = curve_fit(gauss, bins1_centers[65:90], h1[65:90],
#                        p0=[2300, 400, 2750, 100])
# miu = parameters[0][2]
# d_miu = parameters[1][2][2]
# sigma = parameters[0][3]
# d_sigma = parameters[1][3][3]
#
# start = 70
# end = 100
# bins1_centers = bins1_centers[start:end]
# h1 = h1[start:end]
# plt.step(bins1_centers, h1, where='mid')
# plt.plot(bins1_centers, gauss(bins1_centers, parameters[0][0], parameters[0][1], miu, sigma),
#          '--',
#          label='Fit result', )
# plt.xlabel('num photons created')
# plt.ylabel('count')
# plt.title('num photons created by 59.5keV gamma ray in LAr')
# plt.show()
#
# print("gaussian H + A * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))\n with params H = ", parameters[0][0], ", A = ",
#       parameters[0][1], ", $\mu$ = ", miu, ", $\sigma$ = ", sigma)

dat2 = pd.read_csv("../elasticPhot.csv")
num_photons = dat2.numPhotonsCreated

h1, bins1 = np.histogram(num_photons, np.linspace(200, 600, 60))
mean_N_ph = np.mean(num_photons)
std_N_ph = np.std(num_photons)
plt.plot([mean_N_ph, mean_N_ph], [0, np.max(h1)], '--r')
plt.text(mean_N_ph + 2 * std_N_ph, np.max(h1) * 0.8, "mean =%.1f, std = %.1f" % (mean_N_ph, std_N_ph), fontsize=10,
         color='red')
bins1_centers = (bins1[1:] + bins1[:-1]) / 2

# parameters = curve_fit(gauss, bins1_centers, h1,
#                        p0=[0, 350, 420, 15])
# miu = parameters[0][2]
# d_miu = parameters[1][2][2]
# sigma = parameters[0][3]
# d_sigma = parameters[1][3][3]
#
# bins1_centers = bins1_centers
# h1 = h1
plt.step(bins1_centers, h1, where='mid')
# plt.plot(bins1_centers, gauss(bins1_centers, parameters[0][0], parameters[0][1], miu, sigma),
#          '--',
#          label='Fit result', )
plt.xlabel('num photons created')
plt.ylabel('count')
plt.title('num photons created by 45-55 keV Ar recoil nucleus')
plt.show()
#
# print("gaussian H + A * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))\n with params H = ", parameters[0][0], ", A = ",
#       parameters[0][1], ", $\mu$ = ", miu, ", $\sigma$ = ", sigma)

# dat = pd.read_csv("../dEdx.csv")
# # in keV
# init_energy = 1000 * dat.initialEnergy
# dE_dx = dat.dEdx

# energies = np.linspace(0, 10000, 100001)
# avg_dE_dx = np.zeros(100001)
# num_entries = np.zeros(100001)
#
# for i in range(len(dE_dx)):
#     energy = init_energy[i]
#     idx = int(energy // 100)
#     avg_dE_dx[idx] += dE_dx[i]
#     num_entries[idx] += 1
#
# for j in range(len(avg_dE_dx)):
#     avg_dE_dx[j] = avg_dE_dx[j] / num_entries[j]

# energies = energies[avg_dE_dx < 10]
# avg_dE_dx = avg_dE_dx[avg_dE_dx < 10]
# plt.plot(energies, avg_dE_dx, '.')
# plt.xlabel('initial energy[MeV]')
# plt.ylabel('dE/dx[MeV/cm]')
# plt.yscale('log')
# plt.show()

# two_d_plot(init_energy, dE_dx, "title", "initial energy[MeV]", "dE/dx [MeV/cm]", 0, 0, 10)
# plt.hist2d(init_energy, dE_dx, 1000, density=True)
# plt.plot(init_energy, dE_dx, '.')
# plt.yscale('log')
# plt.show()

# plt.hist2d(init_energy, dE_dx, bins=[np.arange(0, 10000, 100), np.arange(0, 10, 0.1)])
# plt.hist2d(init_energy, dE_dx, bins=(np.linspace(0, 1.e7, 50), np.linspace(0, 5000, 50)), cmap='hot_r', norm=LogNorm())
# plt.hist2d(init_energy, dE_dx, bins=(np.linspace(0, 1.e5, 50), np.linspace(0, 500, 50)), cmap='hot_r', norm=LogNorm())
# # plt.xscale('log')
# # plt.yscale('log')
# plt.xlabel('initial energy [MeV]')
# plt.ylabel('dE/dx [MeV/cm]')
# plt.title('dE/dx as a function of proton energy')
# plt.show()


dat = pd.read_csv("../photLoc.csv")
event_ids = dat.eventID
num_photons = np.zeros(1000)
low_events = []
for i in range(len(event_ids)):
    num_photons[event_ids[i]] += 1

for i in range(1000):
    if num_photons[i] < 1000:
        low_events.append(i)
num_high = 1000 - len(low_events)
#
# plt.hist(num_photons, 50)
# plt.show()

xPos = dat.xPos
yPos = dat.yPos
zPos = dat.zPos
# array of mean number of p.e that reached the specific area in the detector
num_pes = [[0 for _ in range(num_high)] for i in range(10)]

# look only at 1 side for now
x = 50
ys = np.linspace(-50, 50, 10)
zs = np.linspace(-50, 50, 10)

k = 0
curr_id = 0
for i in range(len(event_ids)):
    if event_ids[i] not in low_events:
        if event_ids[i] != curr_id:
            k += 1
            curr_id = event_ids[i]
        if i % 100000 == 0:
            print("event ", i)
        if xPos[i] == 50:
            event_id = event_ids[i]
            # y and z positions relative to -50,-50 corner
            y = yPos[i] + 50
            z = zPos[i] + 50

            for j in range(min(int(y / 10), int(z / 10)), 10):
                num_pes[j][k] += 1

coverages = [i ** 2 / 6 for i in range(1, 11)]
num_photoelectrons = [np.mean(num_pes[i]) for i in range(len(num_pes))]
print(num_photoelectrons)
std_photoelectrons = [np.std(num_pes[i]) for i in range(len(num_pes))]
print(std_photoelectrons)

plt.errorbar(coverages, num_photoelectrons, yerr=std_photoelectrons)
plt.xlabel('coverage[%]', fontsize=12)
plt.ylabel('# p.e', fontsize=12)
plt.title('number of p.e as a function of coverage for 5.3MeV neutrons in 1x1x1m LAr cell', fontsize=12)
plt.show()
