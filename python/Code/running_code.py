import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LogNorm

from functions import *

# dat = pd.read_csv("../numPE.csv")
# numPE = dat.numPE
# plt.hist(numPE, np.arange(0, max(numPE), 1))
# plt.show()
#
# dat = pd.read_csv("../lightYield.csv")
# nup = dat.numUp
# ndown = dat.numDown
# ntot = nup + ndown
# plt.hist(ntot, np.arange(0, max(ntot), 2))
# plt.show()

# dat = pd.read_csv("../dEdx.csv")
# dEdx = dat.dEdx  # MeV / mm
# plt.hist(dEdx, np.linspace(min(dEdx) * 0.8, np.max(dEdx), 100))
# plt.xlabel('$\\frac{dE}{dx} [\\frac{MeV}{mm}]$')
# plt.ylabel('count')
# plt.title('$\\frac{dE}{dx}$ of Ar nuclear recoils $E_i < 232 [keV]$')
# plt.show()

# dat = pd.read_csv("../dEdx_43keV.csv")  # this gives dEdx in keV/cm !!!!!!!!!!!!!!!!!!!!
# dEdx = dat.dEdx * 10 ** (-2)  # now dEdx is in MeV/mm !!!!!!!!!!!!!!!!!!!!!!!!
# plt.hist(dEdx, np.linspace(0.1, 400, 100))
# plt.xlabel('$\\frac{dE}{dx} [\\frac{MeV}{mm}]$', fontsize=15)
# plt.ylabel('count', fontsize=15)
# plt.show()
# dEdx = dEdx[0]

dat = pd.read_csv("../dEdX_e.csv")
pathLength = dat.pathLength
numSteps = dat.numSteps
totEnergy = dat.totEnergy

print("mean path length = " + str(np.mean(pathLength)) + " cm")
print("mean number of steps taken by e- = " + str(np.mean(numSteps)))

dx = np.mean(pathLength) * 10 / np.mean(numSteps)  # in mm
print("mean step length = " + str(dx) + " mm")
# print("total energy = " + str(totEnergy) + " keV")
# dEdx = totEnergy[0] / (np.mean(numSteps) * dx * 10 ** 3)  # MeV / mm
dEdx = np.mean(totEnergy / pathLength) / 10
# dEdx = 35
# dEdx = 11.5 / (10 ** 3 * np.mean(pathLength))
dEdx = 0.6
print("dE/dX = " + str(dEdx) + " MeV / mm")
# print(dEdx)
# print(dx)
S = 5.1 * 10 ** 4  # photons / MeV
l = 50
l_att = 50
cr_E = [11.5, 28.5, 43.4, 60.5, 119.5]
electron_pes = [155, 428, 680, 980, 2028]
cr_pe = [12, 25, 49, 70, 150]
sim_E = [11.5, 28.5, 43.4, 119.5]
sim_pe = [11.5, 24, 29, 50]
QEs = [1]
kBs = np.linspace(0.1, 1, 5)
N_steps = np.linspace(0, 1000, 1000)
for QE in QEs:
    for kB in kBs:
        L = dx * N_steps * S * (dEdx / (1 + kB * dEdx))
        N_pe = L * np.e ** (-l / l_att) * QE
        delta_E = N_steps * dEdx * dx * 10 ** 3  # in keV
        plt.plot(delta_E, N_pe, label='kB = ' + str(kB))
    # plt.plot(cr_E, cr_pe, '.', label='Creus data')
    # plt.plot(sim_E, sim_pe, '.', label='current simulation with kB = 0.1')
    plt.plot(cr_E, electron_pes, '.', color='royalblue', markersize=15, label='current simulation with kB = 0.1')
    plt.xlabel('$\Delta_E$ [keV]')
    plt.ylabel('$N_{PE}$')
    plt.title('Calculation results for an (over?) simplified model, QE = ' + str(QE * 100) + " %")
    plt.legend()
    plt.xlim([0, 200])
    plt.ylim([0, np.max(electron_pes) * 1.5])
    # plt.yscale('log')
    plt.show()

# plt.hist(totEnergy / numSteps, 100)
# plt.xlabel('E [MeV]')
# plt.title('average energy transfer to optical photons per step')
# plt.show()
#
# plt.hist(pathLength / numSteps, 100)
# plt.xlabel('distance [cm]')
# plt.title('average path length of Ar nucleus per step')
# plt.show()
#
# plt.hist(numSteps)
# plt.xlabel('#')
# plt.title('number of steps Ar nucleus takes')
# plt.show()

data = pd.read_csv("../Scintillation_Data/Creus/90deg_data.csv")
pes_data = data.num_photons
counts_data = data.counts
bin_size = pes_data[1] - pes_data[0]
I_tot = np.sum(counts_data * bin_size)
counts_data = counts_data
plt.plot(pes_data, counts_data, '.', color='red', label='creus et al. data')

dat = pd.read_csv("../lightYield_90deg_optical_tof.csv")

num_up = dat.numUp
num_down = dat.numDown
num_pes = num_up + num_down
tof = dat.timeOfFlight

h, bins = np.histogram(num_pes, np.arange(0, max(num_pes), bin_size))
I = sum([h[i] * (bins[1] - bins[0]) for i in range(len(h))])
h = h * I_tot / (bin_size * sum(h[0:]))
plt.step((bins[1:] + bins[:-1]) / 2, h, color='red',
         label='simulation with optical photon tof measurement')

# step_from_dist(dists=t1 - t0, bins=np.arange(0, 1000, 1), xlabel="time-of-flight [ns]", ylabel="count", title='')
# step_from_dist([t0, t1], [np.arange(0, 1000, 1) for _ in range(2)], "time [ns]", "count", "",
#                labels=["argon cell signal timing", "liquid scintillator signal timing"], colors=["red", "blue"],
#                yscale=1)
#
dat = pd.read_csv("../lightYield_0.1_birks_90deg_1_fano.csv")
num_up = dat.numUp
num_down = dat.numDown
num_pes = num_up + num_down

h, bins = np.histogram(num_pes, np.arange(0, max(num_pes), bin_size))
h = h * I_tot / (bin_size * sum(h[0:]))
plt.step((bins[1:] + bins[:-1]) / 2, h, color='orange',
         label='simulation using $0.1\\frac{mm}{MeV}$ Birks and 1 fano')

plt.xlabel('# p.e')
plt.ylabel('count')
# plt.xlim([0, 42])
plt.legend()
plt.title('number of p.e following 90 degree neutron scatter off LAr')
plt.show()
#

# n = 5
# t_A = 40
# t_0 = 100
# tau_2 = 1600  # supposed to be measured decay time of slow component on event by event basis
# # for simplicity first use 1600 ns for every event
# t_i = [t_0 + tau_2 * np.log((n + 1) / (n + 1 - i)) for i in range(0, 6)]
#
# sigma_ped = np.array([i for i in range(1, 6)])
# sigma_i = np.array([5 * sigma_ped[i-1] * t_i[i] - t_i[i - 1] for i in range(1, 6)])  # integration error for each bin
# W_i = 1 / sigma_i ** 2  # inverse square of sigma_i
# L_i = np.array([i for i in range(1, 6)])  # integration over each bin
# L = (n + 1) * sum(L_i * W_i) / sum(W_i)
# print(L)
