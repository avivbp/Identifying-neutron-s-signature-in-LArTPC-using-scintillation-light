from cProfile import label
from collections import defaultdict
from fileinput import filename
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
from matplotlib.colors import LogNorm
import os

from Code.functions import L_eff
from functions import *
import scipy as sc
from collections import Counter

dat = pd.read_csv("../current/scatterStats.csv")
ExtScatter = dat.ExtScatter
CryoScatter = dat.CryoScatter
innerLayerScatter = dat.innerLayerScatter
numElasticSensitive, numInelasticSensitive = dat.numElasticSensitive, dat.numInelasticSensitive
outerCellScatterAngle, CryoScatterAngle, outerCellEDep, CryoEDep = (dat.outerCellScatterAngle,
                                                                    dat.CryoScatterAngle, dat.outerCellEDep,
                                                                    dat.CryoEDep)
SensitiveInteraction = [0 for i in range(len(numInelasticSensitive))]
for i in range(len(numInelasticSensitive)):
    if numInelasticSensitive[i] > 0 or numElasticSensitive[i] > 0:
        SensitiveInteraction[i] = 1

x,y = discrete_count(SensitiveInteraction)

plt.bar(x, y, width=0.2, align='center', edgecolor='black')
plt.xlabel('scattered')
plt.ylabel('count')
plt.title('number of events with scatter in fiducial volume')
plt.show()

plt.hist(numElasticSensitive, color='blue',label='number of elastic scatters')
plt.hist(numInelasticSensitive, color='red',label='number of inelastic scatters')
plt.xlabel('number of interactions')
plt.ylabel('count')
plt.title('number of elastic/inelastic interactions in fiducial volume')
plt.legend()
plt.show()

x1,y1 = discrete_count(ExtScatter)
x2,y2 = discrete_count(CryoScatter)
x3,y3 = discrete_count(innerLayerScatter)

plt.bar(x1, y1, width=0.5, align='center', edgecolor='black', label='outer vessel scatter')
plt.bar(x2, y2, width=0.3, align='center', edgecolor='black', label='Cryo scatter')
plt.bar(x3, y3, width=0.15, align='center', edgecolor='black', label='aluminium scatter')
plt.xticks(np.arange(0, np.max(ExtScatter)))
plt.xlabel('scattered')
plt.ylabel('count')
plt.yscale('log')
plt.title('number of events with external scatter')
plt.legend()
plt.show()

plt.hist(outerCellEDep, label='outer cell energy deposit')
plt.hist(CryoEDep, label='Cryo energy deposit')
plt.xlabel('energy deposit [keV]')
plt.ylabel('count')
plt.title("neutron energy deposit not in fiducial volume")
plt.show()

energies = [2.5, 14, 30, 35, 40]
lengths = [150, 250, 350, 450, 550, 650, 1000, 3000]

with open("../current/batch_runs.mac", "w") as f:
    f.write("/control/verbose 0\n/run/verbose 0\n/event/verbose 0\n/tracking/verbose 0\n/run/initialize\n\n")
    for E in energies:
        f.write(f"/gun/energy {E} MeV\n")
        for L in lengths:
            f.write(f"/mysim/setAbsLength {L} cm\n")
            f.write("/run/beamOn 10000\n")

# dat = pd.read_csv("../current/electron_pe_timing_small_det.csv")
# timings = dat.tAbsorbed
# events = dat.eventID
# for event in set(events):
#     print(max(timings[events==event]))
pfs, tot_int = plot_prompt_fraction_graph("../current/photLoc.csv", "muon")


#
# filename1 = "../current/prompt_study/50_50/7_1500_ns/photTimes_mu_100x200micro_1side.csv"
# filename2 = "../current/prompt_study/50_50/7_1500_ns/photTimes_e_100x200micro_1side.csv"
# filename3 = "../current/prompt_study/50_50/7_1500_ns/photTimes_p_100x200micro_1side.csv"
# filename4 = "../current/prompt_study/50_50/7_1500_ns/photTimes_n_100x200micro_1side.csv"
#
# file_list = [filename1, filename2, filename3, filename4]
# particles = [file_list[i].split("/")[5].split("_")[1] for i in range(len(file_list))]
# plot_prompt_fraction_grid(file_list, particles)


#
# pfs, tot_int = plot_prompt_fraction_graph("../current/prompt_study/electron_pe_timing_small_det.csv", "electron")
# pfs2, tot_int2 = plot_prompt_fraction_graph("../current/prompt_study/proton_pe_timing_small_det.csv", "proton")
# pfs3, tot_int3 = plot_prompt_fraction_graph("../current/prompt_study/electron_pe_timing_protoDUNE.csv", "electron")
# pfs4, tot_int4 = plot_prompt_fraction_graph("../current/prompt_study/proton_pe_timing_protoDUNE.csv", "proton")

# avg_pf_electron = np.average(pfs)
# sigma_pf_electron = np.std(pfs)
# avg_pf_proton = np.average(pfs2)
# sigma_pf_proton = np.std(pfs2)


def time_distribution(t, As, At):
    # Example: Bi-exponential (fast and slow components)
    fast_time = 7
    slow_time = 1500
    fast = As / fast_time * np.exp(-t / fast_time)
    slow = At / slow_time * np.exp(-t / slow_time)
    return fast + slow


# Sampling parameters
# t_min = 0
# t_max = 10000  # ns
# fig = plt.figure(figsize=(8, 6))
# bin_width = 10
#
# t_values = np.arange(t_min, t_max, 10)
# y_values = time_distribution(t_values, 0.8, 1 - 0.8)
# plt.plot(t_values, y_values, '.', label='singlet amplitude = 0.8')
# y_values = time_distribution(t_values, 0.1, 1 - 0.1)
# plt.plot(t_values, y_values, '.', color='black', label='singlet amplitude = 0.1')
# equation = r"$\ell(t) = \frac{A_S}{\tau_S} \exp\left(-\frac{t}{\tau_S}\right) + \frac{A_T}{\tau_T} \exp\left(-\frac{t}{\tau_T}\right)$"
# plt.text(0.5, 0.65, equation, fontsize=22, ha='center', va='bottom', transform=plt.gca().transAxes)
# plt.yscale('log')
# plt.xlabel("time [ns]",fontsize=20)
# plt.ylabel("amplitude [a.u.]",fontsize=15)
# plt.yscale('log')
# plt.title("")
# plt.legend()
# fig.savefig("example_dist.png")
# plt.close(fig)
# #
# prompt_window = (-300, 500)  # relative to the peak
#
# amplitudes = np.linspace(0, 1, 100)
# prompt_fractions = []
# for amp in amplitudes:
#     # 1. Sample the function
#     t_values = np.arange(t_min, t_max, bin_width)
#     y_values = time_distribution(t_values, amp, 1 - amp)
#
#     # 2. Normalize the distribution
#     area = simpson(y_values, x=t_values)
#     y_values_normalized = y_values / area
#
#     # 3. Find peak time
#     peak_index = np.argmax(y_values_normalized)
#     peak_time = t_values[peak_index]
#
#     # 4. Define integration windows
#     prompt_start = peak_time + prompt_window[0]
#     prompt_end = peak_time + prompt_window[1]
#
#     # 5. Define masks for integration
#     in_prompt = (t_values >= prompt_start) & (t_values <= prompt_end)
#     in_total = (t_values >= prompt_start)
#
#     # 6. Compute prompt fraction
#     numerator = simpson(y_values_normalized[in_prompt], x=t_values[in_prompt])
#     denominator = simpson(y_values_normalized[in_total], x=t_values[in_total])
#     prompt_fraction = numerator / denominator if denominator > 0 else np.nan
#     prompt_fractions.append(prompt_fraction)
#
# fig = plt.figure(figsize=(8, 6))
# plt.plot(amplitudes, prompt_fractions, '.', label='analytic calculation of singlet amplitude vs prompt fraction')
# a, b = np.polyfit(amplitudes, prompt_fractions, 1)
# print(f'line equation = {a}x + {b}')
# plt.plot(amplitudes, a * amplitudes + b, '-',
#          label='linear fit')
# plt.errorbar(0.8, avg_pf_electron, yerr=sigma_pf_electron, fmt='.', capthick=2, capsize=2, markeredgecolor='k',
#              label="Simulated prompt fraction for electron with singlet amplitude = 0.8")
# plt.errorbar(0.1, avg_pf_proton, yerr=sigma_pf_proton, fmt='.', capthick=2, capsize=2, markeredgecolor='k',
#              label="Simulated prompt fraction for proton with singlet amplitude = 0.1")
# plt.xlabel('singlet relative amplitude',fontsize=20)
# plt.ylabel('prompt fraction',fontsize=20)
# plt.title("")
# plt.legend()
# fig.savefig("analytic_pf.png", dpi=300, bbox_inches='tight')
# plt.close(fig)

# plot_prompt_fraction_graph(filename,"muon")
# file_list = []
# directory = "../current/prompt_study"
# for entry in os.scandir(directory):
#     if entry.is_file():  # check if it's a file
#         print(entry.path)
#     elif entry.is_dir():
#         for dire in os.scandir(entry.path):
#             for entry2 in os.scandir(dire.path):
#                 if entry2.is_file():
#                     # print(entry2.path)
#                     if len(entry2.path.split("/")[2].split("\\")[3].split("_")) == 3:
#                         print(entry2.path)
#                         file_list.append(entry2.path)
#                     particle = entry2.path.split("/")[2].split('\\')[3].split("_")[1]
# if particle =="mu":
#     particle_name = "muon"
# else:
#     particle_name = "neutron"
# plot_prompt_fraction_graph(entry2.path,particle_name)

# plot_prompt_fraction_grid(file_list,"idk",use_name=True)
#
# dat = pd.read_csv("../current/Bi-207_spect.csv")
# energy = dat.energy
# counts = dat.counts
# plt.plot(energy,counts)
# plt.xlabel("channel")
# plt.ylabel('counts')
# plt.show()

# sim_file = "../ARIS_fit/neutrons/sim/numPE_neutron_0.03_birks_150_absLen_0.495_TPBEfic.csv"
# dat = pd.read_csv(sim_file)
# num_pe = dat.numDown + dat.numUp
# num_elastic_scatter = dat.numElasticSensitive
# inelastic_scatter = dat.scatteredInelastically
# detector = dat.detector
# N_ph = dat.numPhotons

# num_pe = num_pe[num_elastic_scatter == 1]
# num_pe = num_pe[inelastic_scatter == 0]
# N_ph = N_ph[num_elastic_scatter == 1]
# N_ph = N_ph[inelastic_scatter == 0]
#
# N_ph0 = list(N_ph[detector == "A0"])
# N_ph1 = list(N_ph[detector == "A1"])
# N_ph2 = list(N_ph[detector == "A2"])
# N_ph3 = list(N_ph[detector == "A3"])
# N_ph4 = list(N_ph[detector == "A4"])
# N_ph5 = list(N_ph[detector == "A5"])
# N_ph6 = list(N_ph[detector == "A6"])
# N_ph7 = list(N_ph[detector == "A7"])

#
# N_phs = [N_ph0, N_ph1, N_ph2, N_ph3, N_ph4, N_ph5, N_ph6, N_ph7]

#
# A0_pe = list(num_pe[detector == "A0"])
# A1_pe = list(num_pe[detector == "A1"])
# A2_pe = list(num_pe[detector == "A2"])
# A3_pe = list(num_pe[detector == "A3"])
# A4_pe = list(num_pe[detector == "A4"])
# A5_pe = list(num_pe[detector == "A5"])
# A5_mult_elastic = (num_pe[detector == "A5"])[num_elastic_scatter > 1]
# A6_pe = list(num_pe[detector == "A6"])
# A6_mult_elastic = (num_pe[detector == "A6"])[num_elastic_scatter > 1]
# A7_pe = list(num_pe[detector == "A7"])
# A7_mult_elastic = (num_pe[detector == "A7"])[num_elastic_scatter > 1]
#
# # print("percentage of multiple elastic scatters for detector A5 = ", 100 * len(A5_mult_elastic) / len(A5_pe))
# # print("percentage of multiple elastic scatters for detector A6 = ", 100 * len(A6_mult_elastic) / len(A6_pe))
# # print("percentage of multiple elastic scatters for detector A7 = ", 100 * len(A7_mult_elastic) / len(A7_pe))
# As = [A0_pe, A1_pe, A2_pe, A3_pe, A4_pe, A5_pe, A6_pe, A7_pe]
#
L_effs = [0.243, 0.258, 0.252, 0.268, 0.286, 0.304, 0.331, 0.349]
energies = [7.1, 13.7, 17.8, 21.7, 40.5, 65.4, 98.1, 117.8]
mius = np.load("../ARIS_fit/neutrons/sim/mius.npy")
sigmas = np.load("../ARIS_fit/neutrons/sim/sigmas.npy")


#
# for i in range(len(mius)):
#     plt.errorbar(energies[i], mius[i], yerr=sigmas[i], fmt='.', capsize=5)
#     plt.ylabel('$N_{PH}$', fontsize=16)
#     plt.xlabel('energy[$keV_{ee}$]', fontsize=16)
#     plt.title(
#         '$\\mu_{PH}^*$ = $E [keV_{ee}]\\times{Scint yield [\\frac{phot}{keV}]\\times{f(kB)}}$ and $\\sigma_{PH}^*$ as a function of photon energy',
#         fontsize=16)


def linear(x, m, n):
    return m * x + n


a, b = np.polyfit(energies, mius, 1)
popt = curve_fit(linear, energies, mius, p0=[a, b], sigma=sigmas)
a = popt[0][0]
b = popt[0][1]

data_dir = "../current/hdProtoDUNE/pes/everything"
pattern = re.compile(r"numPE_([0-9.]+)_cm_absLen_([0-9.]+)_MeV")

# energy -> list of (absLen, filepath)
energy_groups = defaultdict(list)

os.chdir(data_dir)
for filename in os.listdir("."):
    match = pattern.match(filename)
    if match:
        abs_len = float(match.group(1))  # 3000
        energy = float(match.group(2))  # 40
        # full_path = os.path.join(data_dir, filename)
        energy_groups[energy].append((abs_len, filename))

# Sort energies for row order
sorted_energies = sorted(energy_groups.keys())

# Determine column count (max absLen entries per energy)
n_rows = len(sorted_energies)
n_cols = max(len(files) for files in energy_groups.values())

# Create subplot grid
fig, axes = plt.subplots(n_rows, n_cols, figsize=(4 * n_cols, 3 * n_rows), squeeze=False)

for row, energy in enumerate(sorted_energies):
    # Sort files by absLen
    sorted_group = sorted(energy_groups[energy], key=lambda x: x[0])
    for col, (abs_len, filepath) in enumerate(sorted_group):
        # print(filepath)
        ax = axes[row][col]
        dat = pd.read_csv(filepath)
        eventID = dat.eventID
        pes = dat.numPE
        eDep = dat.eDep
        photons = dat.numPhotons
        inelastic = dat.numInelastic
        num_elastic = dat.numElasticSensitive

        # eventID = eventID[inelastic == 0]
        # pes = np.array(pes[inelastic == 0])
        # photons = np.array(photons[inelastic == 0])
        # eDep = np.array(eDep[inelastic == 0])
        # num_elastic = np.array(num_elastic[inelastic == 0])

        eventID = eventID[photons != 100000]
        pes = np.array(pes[photons != 100000])
        eDep = np.array(eDep[photons != 100000])
        num_elastic = np.array(num_elastic[photons != 100000])
        photons = np.array(photons[photons != 100000])

        # print(sorted(list(photons))[-3])

        for j in range(len(pes)):
            if num_elastic[j]:
                avg_energy_deposit = eDep[j] / num_elastic[j]
                if avg_energy_deposit > energies[-1]:
                    index = -1
                else:
                    energies.append(avg_energy_deposit)
                    index = sorted(energies).index(avg_energy_deposit)
                    energies.pop(index)
                N_ph_star = a * avg_energy_deposit + b
                if index == -1:
                    pes[j] = pes[j] * N_ph_star * L_effs[index] / photons[j]
                else:
                    pes[j] = pes[j] * N_ph_star * (L_effs[index - 1] + L_effs[index]) / (2 * photons[j])

        h, bins = np.histogram(pes, bins=np.max(pes) // 2, density=True)
        bins_centers = (bins[1:] + bins[:-1]) / 2
        bins_centers = np.insert(bins_centers, 0, bins[0])  # Left edge of first bin
        h = np.insert(h, 0, 0)  # Zero count before the first bin
        min_idx = np.argwhere(h[1:] <= 0.000001)[0][0] + 1
        min_val = bins[min_idx]
        ax.step(bins_centers, h)
        ax.set_xlim([-1, min_val])
        # ax.set_yscale('log')
        ax.set_xlabel("number of p.e.")
        ax.set_ylabel("counts (normalized)")
        ax.set_title(filepath.split("_")[1] + " cm absLength, " + filepath.split("_")[4] + " MeV neutron")

plt.suptitle("number of p.e. expected in hdProtoDUNE from neutron elastic scatters", fontsize=14)
plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig("numPE_hdProtoDUNE.png", dpi=300)
plt.close(fig)

# xs = np.linspace(0, 120, 100)
# ys = a * xs + b
# plt.plot(xs, ys, '--',label=f"line = {a}x + {b}")
# plt.legend()
# plt.show()
#
# for i in range(len(As)):
#     for j in range(len(As[i])):
#         N_ph_star = np.random.normal(loc=mius[i], scale=sigmas[i])
#         As[i][j] = As[i][j] * N_ph_star * L_effs[i] / N_phs[i][j]
#
# filename0 = "../ARIS_fit/neutrons/data/A0_data.csv"
# filename1 = "../ARIS_fit/neutrons/data/A1_data.csv"
# filename2 = "../ARIS_fit/neutrons/data/A2_data.csv"
# filename3 = "../ARIS_fit/neutrons/data/A3_data.csv"
# filename4 = "../ARIS_fit/neutrons/data/A4_data.csv"
# filename5 = "../ARIS_fit/neutrons/data/A5_data.csv"
# filename6 = "../ARIS_fit/neutrons/data/A6_data.csv"
# filename7 = "../ARIS_fit/neutrons/data/A7_data.csv"
#
# init_guesses = [[40, 70, 10, 2], [15, 40, 20, 3], [30, 140, 25, 3], [60, 200, 35, 5], [60, 130, 65, 6],
#                 [30, 60, 120, 15], [40, 110, 200, 20], [20, 200, 185, 30]]
# data_filenames = [filename0, filename1, filename2, filename3, filename4, filename5, filename6, filename7]
#
# for i in range(len(As)):
#     plot_numPE(As[i], data_filenames[i])
