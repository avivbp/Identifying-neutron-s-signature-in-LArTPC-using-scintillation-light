from fileinput import filename

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
from matplotlib.colors import LogNorm
import os
from functions import *
import scipy as sc

# os.chdir("../current/ARIS_neutron")
# for file in os.listdir("."):
#     dat = pd.read_csv(file)
#     print(file)
#     N_ph = dat.numPhotons
#     plt.hist(N_ph[N_ph != 0], 100)
#     plt.show()

# dat = pd.read_csv("../current/numPH_25deg.csv")
# N_ph = dat.numPhotons
# N_ph = N_ph[N_ph != 0]
# plt.hist(N_ph, 100)
# plt.show()

# os.chdir("../current/")
# for file in os.listdir("."):
#     dat = pd.read_csv(file)
#     print(file)
#     N_ph = dat.numPhotons
#     plt.hist(N_ph[N_ph != 0], 70)
#     plt.show()

dat = pd.read_csv("../current/lightYield.csv")
num_pe = dat.numDown + dat.numUp
num_elastic_scatter = dat.numElasticSensitive
inelastic_scatter = dat.scatteredInelastically
external_scatter = dat.scatteredNotSensitive
detector = dat.detector
recoil_energy = dat.nucleusRecoilEnergy
#
num_pe = num_pe[num_elastic_scatter == 1]
num_pe = num_pe[inelastic_scatter == 0]
num_pe = num_pe[external_scatter == 0]

N_ph = dat.numPhotonsSensitive
num_pe = num_pe[N_ph != 0]
N_ph = N_ph[N_ph != 0]
N_ph = N_ph[num_elastic_scatter == 1]
N_ph = N_ph[inelastic_scatter == 0]
N_ph = N_ph[external_scatter == 0]

N_ph0 = list(N_ph[detector == "A0"])
N_ph1 = list(N_ph[detector == "A1"])
N_ph2 = list(N_ph[detector == "A2"])
N_ph3 = list(N_ph[detector == "A3"])
N_ph4 = list(N_ph[detector == "A4"])

A0_pe = list(num_pe[detector == "A0"])
A1_pe = list(num_pe[detector == "A1"])
A2_pe = list(num_pe[detector == "A2"])
A3_pe = list(num_pe[detector == "A3"])
A4_pe = list(num_pe[detector == "A4"])

As = [A0_pe, A1_pe, A2_pe, A3_pe, A4_pe]
N_phs = [N_ph0, N_ph1, N_ph2, N_ph3, N_ph4]

energies = [11.5, 28.5, 43.3, 60.5, 119.5]
# path_list = []
# path = "../current/"
# for file in os.listdir(path):
#     if not file.startswith("light"):
#         path_list.append(path + file)
# compare_gaussian(path_list, [[450, 650], [1200, 1600], [1900, 2400], [2700, 3300], [5500, 6500]],
#                  [[0, 500, 567, 20], [0, 500, 1415, 40], [0, 500, 2160, 70], [0, 500, 3022, 90], [0, 500, 6048, 160]], 30,
#                  "numPhotons")
N_ph_star = [567, 1415, 2160, 3022, 6048]
mius = np.load("../current/mius.npy")
sigmas = np.load("../current/sigmas.npy")

for i in range(len(mius)):
    plt.errorbar(energies[i], mius[i], yerr=sigmas[i], fmt='.',capsize=5)
    plt.ylabel('$N_{PH}$',fontsize=16)
    plt.xlabel('energy[$keV_{ee}$]',fontsize=16)
    plt.title('$\\mu_{PH}^*$ = $E [keV_{ee}]\\times{Scint yield [\\frac{phot}{keV}]\\times{f(kB)}}$ and $\\sigma_{PH}^*$ as a function of photon energy',fontsize=16)
    plt.legend()
a, b = np.polyfit(energies, mius, 1)
xs = np.linspace(0, 120, 100)
ys = a * xs + b
plt.plot(xs, ys, '--', linewidth=1)
plt.show()

for i in range(len(As)):
    for j in range(len(As[i])):
        N_ph_star = np.random.normal(loc=mius[i], scale=sigmas[i])
        As[i][j] = As[i][j] * N_ph_star * 0.28 / N_phs[i][j]

filename0 = "../Scintillation_Data/Creus/25deg_data.csv"
filename1 = "../Scintillation_Data/Creus/40deg_data.csv"
filename2 = "../Scintillation_Data/Creus/50deg_data.csv"
filename3 = "../Scintillation_Data/Creus/60deg_data.csv"
filename4 = "../Scintillation_Data/Creus/90deg_data.csv"

data_filenames = [filename0, filename1, filename2, filename3, filename4]
errors_filenames = []

for i in range(len(data_filenames)):
    err_down_file = data_filenames[i].replace("data", "err_down")
    err_up_file = data_filenames[i].replace("data", "err_up")
    errors_filenames.append((err_down_file, err_up_file))

# for i in range(len(As)):
#     # print(As[i])
#     plot_numPE(As[i], data_filenames[i], filenames_errors=errors_filenames[i])

# filename = "../ARIS_fit/Na/numPE_gamma_0.511_E_0.03_birks_150_absLen_0.495_TPBEfic.csv"
# dat = pd.read_csv(filename)
# numPE = dat.numPE
# numPE = numPE[numPE != 0]
# if filename.split("_")[2] == "0.511":
#     h, bins = np.histogram(numPE[numPE > 3100], 25)
#     bins_center = (bins[1:] + bins[:-1]) / 2
#     parameters = curve_fit(gauss, bins_center, h,
#                            p0=[0, 45, 3300, 60])
#     miu = parameters[0][2]
#     H = parameters[0][0]
#     A = parameters[0][1]
#     d_miu = parameters[1][2][2]
#     sigma = parameters[0][3]
#     d_sigma = parameters[1][3][3]
#     # print("miu = " + str(miu)+" +- " + str(d_miu) +", sigma = " + str(sigma))
#     plt.plot(bins_center, gauss(bins_center, H, A, miu, sigma),
#              '--',
#              label='Fit result : $\\mu = $ ' + "{0:.2f}".format(miu) + ", $\\sigma = $ " + "{0:.2f}".format(sigma))
#     plt.step((bins[1:] + bins[:-1]) / 2, h, label='Geant4 simulation')
#     xs = [3300 for _ in range(np.max(h))]
#     ys = np.linspace(0, np.max(h), np.max(h))
#     plt.plot(xs, ys, '--', label="$^{22} Na$ p.e. spectrum")
#     plt.xlabel('S1 [p.e]')
#     plt.ylabel('counts')
#     plt.title('$^{22}Na $ p.e. spectrum')
#     plt.legend()
#     plt.show()
# plot_numPE(numPE, "../ARIS_fit/Am/Am_data.csv", label='$^{241} Am$ p.e. spectrum',
#            gaussian_init_guess=[0, 3500, 365, 30])
#
# my_pes = [374, 3322]
# my_sigmas = [20, 47]
# ARIS_pes = [365, 3296]
# ARIS_sigmas = [34, 30]
# a, b = np.polyfit(ARIS_pes, my_pes, 1)
# xs = np.linspace(0, 3500, 3500)
# plt.plot(xs, a * xs + b, '--', linewidth=3, label='linear fit')
# plt.errorbar(ARIS_pes, my_pes, yerr=my_sigmas, xerr=ARIS_sigmas, fmt='o')
# plt.xlabel('ARIS data mean p.e.')
# plt.ylabel('Geant4 simulation mean p.e.')
# plt.title('comparison of photopeaks between Geant4 simulation and ARIS data')
# plt.xticks(np.linspace(0, 3500, 10))
# plt.legend()
# plt.show()

path_list = []
path = "../ARIS_fit/neutrons/sim/"
for file in os.listdir(path):
    if not file.startswith("numPE_neutron") and not file.startswith("miu") and not file.startswith("sigmas"):
        path_list.append(path + file)
compare_gaussian(path_list, [[250, 350], [500, 700], [750, 950], [900, 1140], [1600, 2300], [2700, 3500], [4500, 5200],
                             [5400, 6200]],
                 [[0, 500, 296, 20], [0, 500, 620, 20], [0, 500, 828, 30], [0, 500, 1020, 40], [0, 500, 1973, 60],
                  [0, 500, 3226, 70], [0, 500, 4864, 100], [0, 500, 5832, 150]], 30, "numPhotons")
# compare_gaussian(["../Darkside_fit/Na22/0.44tpb/numPE_Na22_0.007_3000_100k.csv",
#                   "../Darkside_fit/Na22/0.44tpb/numPE_Na22_0.007_3000_100k.csv"], [[3800, 5000], [10000, 12000]],
#                  [[0.3, 0.25, 4486, 19], [0.1, 0.3, 11000, 300]], 200)
# compare_gaussian(["../Darkside_fit/Co57/0.44tpb/numPE_Co57_0.007_3000.csv"], [[950, 1200]], [[0.06, 0.06, 1100, 60]],
#                  40)
# compare_gaussian(["../Darkside_fit/Cs137/0.44tpb/numPE_Cs137_0.007_3000.csv"], [[5000, 6000]], [[0.05, 0.1, 5650, 100]],
#                  200)


# sim_mius = [4379, 10978, 1033, 5667]
# sim_sigmas = [72, 120, 34, 81]
# data_mius = [4486, 10096, 1082, 6010]
# data_sigmas = [153, 318, 57, 186]
# a, b = np.polyfit(data_mius, sim_mius, 1)
#
# plt.errorbar(data_mius, sim_mius, xerr=data_sigmas, yerr=sim_sigmas, fmt='o', )
# xs = np.linspace(0, 11000, 1000)
# plt.plot(xs, a * xs + b, '--', linewidth=3, label='linear fit')
# plt.xlabel('Darkside data mean p.e.')
# plt.ylabel('Geant4 simulation mean p.e.')
# plt.title('comparison of photopeaks between Geant4 simulation and Darkside data')
# plt.xticks(np.linspace(0,11000,12))
# plt.legend()
# plt.show()

#################################################################################################################
sim_file = "../ARIS_fit/neutrons/sim/numPE_neutron_0.03_birks_150_absLen_0.495_TPBEfic.csv"
dat = pd.read_csv(sim_file)
num_pe = dat.numDown + dat.numUp
num_elastic_scatter = dat.numElasticSensitive
inelastic_scatter = dat.scatteredInelastically
detector = dat.detector
N_ph = dat.numPhotons

# num_pe = num_pe[num_elastic_scatter == 1]
# num_pe = num_pe[inelastic_scatter == 0]
#
# N_ph = N_ph[num_elastic_scatter == 1]
# N_ph = N_ph[inelastic_scatter == 0]

N_ph0 = list(N_ph[detector == "A0"])
N_ph1 = list(N_ph[detector == "A1"])
N_ph2 = list(N_ph[detector == "A2"])
N_ph3 = list(N_ph[detector == "A3"])
N_ph4 = list(N_ph[detector == "A4"])
N_ph5 = list(N_ph[detector == "A5"])
N_ph6 = list(N_ph[detector == "A6"])
N_ph7 = list(N_ph[detector == "A7"])

N_phs = [N_ph0, N_ph1, N_ph2, N_ph3, N_ph4, N_ph5, N_ph6, N_ph7]

A0_pe = list(num_pe[detector == "A0"])
A1_pe = list(num_pe[detector == "A1"])
A2_pe = list(num_pe[detector == "A2"])
A3_pe = list(num_pe[detector == "A3"])
A4_pe = list(num_pe[detector == "A4"])
A5_pe = list(num_pe[detector == "A5"])
A5_mult_elastic = (num_pe[detector == "A5"])[num_elastic_scatter > 1]
A6_pe = list(num_pe[detector == "A6"])
A6_mult_elastic = (num_pe[detector == "A6"])[num_elastic_scatter > 1]
A7_pe = list(num_pe[detector == "A7"])
A7_mult_elastic = (num_pe[detector == "A7"])[num_elastic_scatter > 1]

# print("percentage of multiple elastic scatters for detector A5 = ", 100 * len(A5_mult_elastic) / len(A5_pe))
# print("percentage of multiple elastic scatters for detector A6 = ", 100 * len(A6_mult_elastic) / len(A6_pe))
# print("percentage of multiple elastic scatters for detector A7 = ", 100 * len(A7_mult_elastic) / len(A7_pe))
As = [A0_pe, A1_pe, A2_pe, A3_pe, A4_pe, A5_pe, A6_pe, A7_pe]

L_effs = [0.243, 0.258, 0.252, 0.268, 0.286, 0.304, 0.331, 0.349]
energies = [7.1, 13.7, 17.8, 21.7, 40.5, 65.4, 98.1, 117.8]
mius = np.load("../ARIS_fit/neutrons/sim/mius.npy")
sigmas = np.load("../ARIS_fit/neutrons/sim/sigmas.npy")

for i in range(len(mius)):
    plt.errorbar(energies[i], mius[i], yerr=sigmas[i], fmt='.', capsize=5)
    plt.ylabel('$N_{PH}$', fontsize=16)
    plt.xlabel('energy[$keV_{ee}$]', fontsize=16)
    plt.title(
        '$\\mu_{PH}^*$ = $E [keV_{ee}]\\times{Scint yield [\\frac{phot}{keV}]\\times{f(kB)}}$ and $\\sigma_{PH}^*$ as a function of photon energy',
        fontsize=16)
    plt.legend()
a, b = np.polyfit(energies, mius, 1)
xs = np.linspace(0, 120, 100)
ys = a * xs + b
plt.plot(xs, ys, '--')
plt.show()

for i in range(len(As)):
    for j in range(len(As[i])):
        N_ph_star = np.random.normal(loc=mius[i], scale=sigmas[i])
        As[i][j] = As[i][j] * N_ph_star * L_effs[i] / N_phs[i][j]

filename0 = "../ARIS_fit/neutrons/data/A0_data.csv"
filename1 = "../ARIS_fit/neutrons/data/A1_data.csv"
filename2 = "../ARIS_fit/neutrons/data/A2_data.csv"
filename3 = "../ARIS_fit/neutrons/data/A3_data.csv"
filename4 = "../ARIS_fit/neutrons/data/A4_data.csv"
filename5 = "../ARIS_fit/neutrons/data/A5_data.csv"
filename6 = "../ARIS_fit/neutrons/data/A6_data.csv"
filename7 = "../ARIS_fit/neutrons/data/A7_data.csv"

init_guesses = [[40, 70, 10, 2], [15, 40, 20, 3], [30, 140, 25, 3], [60, 200, 35, 5], [60, 130, 65, 6],
                [30, 60, 120, 15], [40, 110, 200, 20], [20, 200, 185, 30]]
data_filenames = [filename0, filename1, filename2, filename3, filename4, filename5, filename6, filename7]

for i in range(len(As)):
    plot_numPE(As[i], data_filenames[i])
