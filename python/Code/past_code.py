from functions import *

# files = ["free_path_dist_1MeV.csv", "free_path_dist_5MeV.csv", "free_path_dist_10MeV.csv", "free_path_dist_20MeV.csv",
#          "free_path_dist_50MeV.csv", "free_path_dist_100MeV.csv"]
# mean_free_paths = []
# stds = []
# for file in files:
#     distances = pd.read_csv(file).distance
#     mfp = np.mean(distances)
#     mean_free_paths.append(mfp)
#     std = np.std(distances)
#     stds.append(std)
# energies = [1, 5, 10, 20, 50, 100] # still need to upload runs of 10000 events with 20,50,100 MeV
# print(mean_free_paths)
# print(stds)
# plt.errorbar(energies, mean_free_paths, yerr=stds)
# plt.xlabel("photon energy [MeV]")
# plt.ylabel("mean free path [cm]")
# plt.title("mean free path vs. energy plot")
# plt.show()


# example of extraction of attenuation coefficient

# xs = np.array([0.2, 0.3, 0.25, 0.275, 0.225, 0.175, 0.15])
# ys = np.array([0.359, 0.2162, 0.2801, 0.2453, 0.3186, 0.4104, 0.4671])
# x = np.linspace(0, 1, 100)
# y = model_func(x, 5.098)
# plt.xlabel("Absorber thickness [cm]")
# plt.ylabel("N(transmitted)/N0")
# plt.title("percentage of transmitted particles vs absorber thickness - 10keV")
# plt.plot(xs, ys, '.', label='sampled thicknesses')
# plt.plot(x, y, '-', label='fitted function')
# straight_line = np.linspace(0, 1, 100)
# straight_line_y = [1 / np.e for i in range(100)]
# plt.plot(straight_line, straight_line_y, label='$1/\mu$ fit results')
# plt.legend()
# plt.show()

# filename = "references/Default_Dataset_diff.csv"
# points = pd.read_csv(filename)
# xs = points.x
# ys = points.y
# #
# filename2 = "Default_Dataset2.csv"
# points2 = pd.read_csv(filename2)
# xs2 = points2.x
# ys2 = points2.y
#
# for i in range(len(xs2)):
#     xs2[i] *= 1000

# in keV
# energies = [10, 20, 50, 75, 100, 250, 500, 1000]
# mius = [5.098, 0.7053, 0.186, 0.156, 0.143, 0.1098, 0.0830, 0.0627]
# stds = np.array([0.00583, 0.00194, 0.00197, 0.00160, 0.00299, 0.000598, 0.00148, 0.00040])

# mius = [5.104, 0.733, 0.202, 0.171, 0.157, 0.117, 0.0916, 0.0672]
# stds = np.array([0.00643, 0.00195, 0.00113, 0.00134, 0.000823, 0.00077, 0.000409, 0.000123])
#
# plt.plot(xs, ys, '-',
#          label='$\mu = \mu_\\tau + \mu_\sigma + \mu_{\sigma R} + \mu_\kappa $ [Shehata et al. ,ME Conf.Proc.4(2007)]')
# plt.errorbar(energies, mius, yerr=stds, fmt='.', capthick=2, capsize=2, markeredgecolor='k',
#              label='VALIDATION - Geant4 simulated $\gamma$ attenuation coefficient in water')
# plt.plot(xs2, ys2, '-', label='Grupen, Tum. Ther. Part. Beams. 538.(2000)')
# # plt.legend(loc='upper right', fontsize=9)
# plt.grid()
# plt.xscale("log")
# plt.yscale("log")
# plt.xlabel("Photon energy [keV]")
# plt.ylabel("Linear attenuation coefficient [1/cm]")
# plt.title("Linear attenuation coefficient vs photon energy")
# # plt.show()
#
# arr = [1.00000E-03, 4.078E+03
#     , 1.50000E-03, 1.376E+03
#     , 2.00000E-03, 6.173E+02
#     , 3.00000E-03, 1.929E+02
#     , 4.00000E-03, 8.278E+01
#     , 5.00000E-03, 4.258E+01
#     , 6.00000E-03, 2.464E+01
#     , 8.00000E-03, 1.037E+01
#     , 1.00000E-02, 5.329E+00
#     , 1.50000E-02, 1.673E+00
#     , 2.00000E-02, 8.096E-01
#     , 3.00000E-02, 3.756E-01
#     , 4.00000E-02, 2.683E-01
#     , 5.00000E-02, 2.269E-01
#     , 6.00000E-02, 2.059E-01
#     , 8.00000E-02, 1.837E-01
#     , 1.00000E-01, 1.707E-01
#     , 1.50000E-01, 1.505E-01
#     , 2.00000E-01, 1.370E-01
#     , 3.00000E-01, 1.186E-01
#     , 4.00000E-01, 1.061E-01
#     , 5.00000E-01, 9.687E-02
#     , 6.00000E-01, 8.956E-02
#     , 8.00000E-01, 7.865E-02
#     , 1.00000E+00, 7.072E-02
#     , 1.25000E+00, 6.323E-02
#     , 1.50000E+00, 5.754E-02
#     , 2.00000E+00, 4.942E-02
#     , 3.00000E+00, 3.969E-02
#     , 4.00000E+00, 3.403E-02
#     , 5.00000E+00, 3.031E-02
#     , 6.00000E+00, 2.770E-02
#     , 8.00000E+00, 2.429E-02
#     , 1.00000E+01, 2.219E-02
#     , 1.50000E+01, 1.941E-02
#     , 2.00000E+01, 1.813E-02]
#
# energies = []
# atten_coeffs = []
# for i in range(len(arr)):
#     if i % 2 == 0:
#         energies.append(arr[i] * 1000)
#     else:
#         atten_coeffs.append(arr[i])
#
# plt.plot(energies, atten_coeffs, '-', label='NIST Standard Reference Database 126')
# plt.errorbar(energies2, mius, yerr=3 * stds, fmt='.', capthick=2, capsize=2, markeredgecolor='k',
#              label='This work - simulated $\gamma$ attenuation')
# plt.legend(loc='upper right', fontsize=8)
# plt.grid()
# plt.xscale("log")
# plt.yscale("log")
# plt.xlabel("Photon energy [keV]")
# plt.ylabel("Linear attenuation coefficient [1/cm]")
# plt.title("Linear attenuation coefficient vs photon energy")
# plt.xlim(5,1200)
# plt.show()

# in keV
# energies_lAr = np.array([10, 20, 50, 75, 100, 250, 500, 1000, 10000, 50000, 100000, 1000000])
# energies_lAr = energies_lAr/1000
# mius_lAr = [87.154, 11.765, 0.862, 0.372, 0.247, 0.136, 0.103, 0.0761, 0.0333, 0.0375, 0.0423, 0.0533]
# stds_lAr = [0.0735, 0.0172, 0.00239, 0.00370, 0.00132,
#             0.000706, 0.000417, 0.000220, 8.25e-05, 0.000171, 0.000183, 0.000638]
# plt.errorbar(energies_lAr, mius_lAr, yerr=stds_lAr, fmt='.', capthick=2, capsize=2, markeredgecolor='k',
#              color='blue', label='Geant4 simulated $\gamma$ attenuation coefficient in lAr')
# plt.legend(loc='upper right', fontsize=7)
# plt.xlabel("Photon energy [MeV]")
# plt.ylabel("Linear attenuation coefficient [1/cm]")
# plt.grid()
# plt.xscale("log")
# plt.yscale("log")
# plt.title("Linear attenuation coefficient vs photon energy in lAr")
# plt.show()

# file_name = "table.csv"
# data = pd.read_csv(file_name, sep=';')
# energies = data.energy
# cross_sections = data.cx
# cross_sections_err = data.cx_err
#
# for i in range(len(energies)):
#     energies[i] = float(energies[i][:-3]) / 10 ** 6
#     cross_sections[i] = float(cross_sections[i][:-2])
#     cross_sections_err[i] = float(cross_sections_err[i][:-2])
#
# plt.errorbar(energies[1500:1700], cross_sections[1500:1700], fmt='-', yerr=cross_sections_err[1500:1700],
#              label='Winters et al., Phys. Rev. C 43(1991)')

# in barns and MeV
# sampled_energies = [0.855, 0.86, 0.865, 0.87, 0.875, 0.88, 0.885, 0.89, 0.895, 0.9]
# sampled_cross_sections = [0.944, 3.810, 1.727, 1.030, 4.389, 2.192, 1.642, 1.784, 1.039, 3.995]
# sampled_stds = [0.0670, 0.135, 0.0907, 0.070, 0.145, 0.102, 0.0884, 0.0922, 0.0703, 0.138]

# plt.errorbar(sampled_energies, sampled_cross_sections, fmt='.', yerr=sampled_stds,
#              label='this work - total cross section', capthick=2,
#              capsize=2, markeredgecolor='k')
# plt.xlabel("Neutron Energy [MeV]")
# plt.ylabel("cross section [barn]")
# plt.title("Neutron energy vs cross section")
# plt.legend()
# plt.xlim([0.8, 1])
# plt.ylim([0, 5])
# plt.show()
# elastic_energies = [1, 9, 10, 11, 20, 30, 40, 50]
# elastic_cross_section = [3.156, 1.062, 0.916, 0.890, 0.755, 1.341, 1.415, 1.385]
# elastic_std = [0.0388, 0.0225, 0.0209, 0.0206, 0.0189, 0.0253, 0.0259, 0.0257]
#
# n_p_energies = [9, 10, 11, 20]
# n_p_cross_section = [0.001898, 0.00807, 0.00949, 0.0285]
# n_p_std = [0.000949, 0.00196, 0.00212, 0.00368]
#
# n_n_prime_energies = [9, 10, 11, 20]
# n_n_prime_cross_section = [1.597, 1.574, 1.442, 0.236]
# n_n_prime_std = [0.0276, 0.0274, 0.0262, 0.0106]
#
# n_2n_energies = [9, 11, 20, 30, 40, 50]
# n_2n_cross_section = [0, 0.137, 0.969, 0.235, 0.183, 0.143]
# n_2n_std = [0, 0.00807, 0.0215, 0.0106, 0.00931, 0.00823]
#
# n_3n_energies = [30, 40, 50]
# n_3n_cross_section = [0.531, 0.370, 0.262]
# n_3n_std = [0.0159, 0.0133, 0.0112]
#
# more_energies = [1, 2, 9, 10, 11, 20, 30, 40, 50]
# more_cross = [2.899, 4.873, 2.674, 2.522, 2.494, 2.002, 2.358, 2.296, 2.163]
# more_std = [0.118, 0.153, 0.0357, 0.0347, 0.0345, 0.0977, 0.0335, 0.0331, 0.0321]
# plt.errorbar(sampled_energies + more_energies, sampled_cross_sections + more_cross, fmt='--',
#              yerr=sampled_stds + more_std,
#              label='Geant4 $\sigma_{total}$', capthick=2,
#              capsize=2,
#              markeredgecolor='k')
# plt.errorbar(elastic_energies, elastic_cross_section, fmt='--',
#              yerr=elastic_std,
#              label='Geant4 $\sigma_{elastic}$', capthick=2,
#              capsize=2,
#              markeredgecolor='k')
# plt.errorbar(n_n_prime_energies, n_n_prime_cross_section, fmt='--',
#              yerr=n_n_prime_std,
#              label="Geant4 $\sigma_{^{40}Ar(n,n')^{40}*Ar}$ ", capthick=2,
#              capsize=2,
#              markeredgecolor='k')
# plt.errorbar(n_2n_energies, n_2n_cross_section, fmt='--',
#              yerr=n_2n_std,
#              label="Geant4 $\sigma_{^{40}Ar(n,2n)^{39}*Ar}$ ", capthick=2,
#              capsize=2,
#              markeredgecolor='k')
# plt.errorbar(n_3n_energies, n_3n_cross_section, fmt='--',
#              yerr=n_3n_std,
#              label="Geant4 $\sigma_{^{40}Ar(n,3n)^{38}*Ar}$", capthick=2,
#              capsize=2,
#              markeredgecolor='k')
# plt.errorbar(n_p_energies, n_p_cross_section, fmt='--',
#              yerr=n_p_std,
#              label="Geant4 $\sigma_{^{40}Ar(n,p)^{40}*Cl}$", capthick=2,
#              capsize=2,
#              markeredgecolor='k')
# xs = [9.895 for i in range(100)]
# ys = np.linspace(0, 8, 100)
# plt.plot(xs, ys, '--',
#          label='$S_n$, $S_p$ of $^{40}Ar$ - JUN CHEN Publication cut-off:\n 30-Sep-2015 ENSDF insertion: 2017-02 Publication:\n Nuclear Data Sheets 140, 1 (2017)',
#          color="black")
#
# xs = [12.529 for i in range(100)]
# ys = np.linspace(0, 8, 100)
# plt.plot(xs, ys, '--', color="black")
#
# plt.xlabel("Neutron Energy [MeV]")
# plt.ylabel("cross section [barn]")
# plt.title("n - $^{40}$Ar Cross Section Decomposition")
# plt.xscale('log')
# plt.yscale('log')
# plt.legend(fontsize=6)
# plt.show()

# in MeV
# energies_lAr = np.array([0.1, 0.15, 0.25, 0.255, 0.26, 0.275, 0.3, 0.4, 0.5, 0.75, 0.85, 0.858, 0.875, 0.9, 1, 10, 50])
# mius_lAr = [0.045, 0.0103, 0.00689, 0.156, 0.101, 0.0577, 0.0410, 0.0331, 0.0264,
#             0.0251, 0.0218, 0.0552, 0.0892, 0.0822, 0.0642, 0.0512, 0.0429]
# stds_lAr = [0.000145, 0.0000492, 0.0000184, 0.00112, 0.000923, 0.000244,
#             0.000182, 0.000153, 0.000159, 0.0000757, 0.0000713,
#             0.000257, 0.000437, 0.000483,
#             0.000198, 0.000322, 0.000376]
# plt.errorbar(energies_lAr, mius_lAr, yerr=stds_lAr, fmt='.', capthick=2, capsize=2, markeredgecolor='k')
# plt.xlabel("Neutron energy [MeV]")
# plt.ylabel("Linear attenuation coefficient [1/cm]")
# plt.grid()
# extra_ticks = [0.05]
# plt.yticks(list(plt.yticks()[0]) + extra_ticks)
# plt.xscale("log")
# plt.yscale("log")
# plt.title("Linear attenuation coefficient vs neutron energy in lAr")
# plt.show()

# reactions = pd.read_csv("reactions_20MeV.csv")
# neutron_num = reactions.neutronNum
# proton_num = reactions.protonNum
# gamma_num = reactions.gammaNum
# Z = reactions.Z
# A = reactions.A
#
# min_A = min(A)
# max_A = max(A)
# bins_y = np.arange(min_A, max_A + 2) - 0.5
# min_Z = min(Z)
# max_Z = max(Z)
# bins_x = np.arange(min_Z, max_Z + 2) - 0.5
# fig, ax = plt.subplots()
# hist, xbins, ybins, im = ax.hist2d(Z, A, bins=[bins_x, bins_y], cmap=plt.cm.jet)
# plt.xticks( np.arange(min_Z, max_Z + 1, 1))
# for i in range(3):
#     for j in range(6):
#         ax.text(xbins[i] + 0.5, ybins[j] + 0.5, s="{:.1f}".format(hist[i, j] * 100 / len(Z)) + "%"
#                 , color="w", ha="center",
#                 va="center", fontweight="bold")
# plt.title("2D histogram of Z vs A")
# plt.xlabel("Z of heavy nucleus")
# plt.ylabel("A of heavy nucleus")
# plt.show()

# protons_when_zero_neutrons = []
# protons_when_one_neutron = []
# protons_when_two_neutrons = []
# protons_when_three_neutrons = []
# gammas_when_zero_neutrons = []
# gammas_when_one_neutron = []
# gammas_when_two_neutrons = []
# gammas_when_three_neutrons = []
#
# for i in range(len(neutron_num)):
#     if neutron_num[i] == 0:
#         protons_when_zero_neutrons.append(proton_num[i])
#         gammas_when_zero_neutrons.append(gamma_num[i])
#
#     elif neutron_num[i] == 1:
#         protons_when_one_neutron.append(proton_num[i])
#         gammas_when_one_neutron.append(gamma_num[i])
#
#     elif neutron_num[i] == 2:
#         protons_when_two_neutrons.append(proton_num[i])
#         gammas_when_two_neutrons.append(gamma_num[i])
#
#     elif neutron_num[i] == 3:
#         protons_when_three_neutrons.append(proton_num[i])
#         gammas_when_three_neutrons.append(gamma_num[i])

# print(list(neutron_num).count(3))
# plt.hist(protons_when_zero_neutrons)
# plt.xlabel("number of ejected protons")
# plt.ylabel("number of events")
# plt.title("number of ejected protons when zero neutrons ejected")
# plt.show()
# plt.hist(protons_when_one_neutron)
# plt.xlabel("number of ejected protons")
# plt.ylabel("number of events")
# plt.title("number of ejected protons when one neutron ejected")
# plt.show()
# plt.hist(protons_when_two_neutrons)
# plt.xlabel("number of ejected protons")
# plt.ylabel("number of events")
# plt.title("number of ejected protons when two neutrons ejected")
# plt.show()
# plt.hist(protons_when_three_neutrons)
# plt.xlabel("number of ejected protons")
# plt.ylabel("number of events")
# plt.title("number of ejected protons when three neutrons ejected")
# plt.show()
#
# plt.hist(gammas_when_zero_neutrons)
# plt.xlabel("number of ejected photons")
# plt.ylabel("number of events")
# plt.title("number of ejected photons when zero neutrons ejected")
# plt.show()
# plt.hist(gammas_when_one_neutron)
# plt.xlabel("number of ejected photons")
# plt.ylabel("number of events")
# plt.title("number of ejected photons when one neutron ejected")
# plt.show()
# plt.hist(gammas_when_two_neutrons)
# plt.xlabel("number of ejected photons")
# plt.ylabel("number of events")
# plt.title("number of ejected photons when two neutrons ejected")
# plt.show()
# plt.hist(gammas_when_three_neutrons)
# plt.xlabel("number of ejected photons")
# plt.ylabel("number of events")
# plt.title("number of ejected photons when three neutrons ejected")
# plt.show()

# process_list = list(pd.read_csv("reactions_9MeV.csv").reaction)
#
# processes = ["40Ar(n,n')40*Ar", "40Ar(n,2n)39*Ar", "40Ar(n,3n)38*Ar", "elastic", "40Ar(n,p)40*Cl",
#              "40Ar(n,Ngamma)40*Ar", "other", "40Ar(n,n' p)39*Cl", "40Ar(n,$\\alpha$)37*S", "40Ar(n,2n p)38*Cl",
#              "40Ar(n,3n p)37*Cl", "40Ar(n,4n)37*Ar"]
#
# frequencies = np.array([process_list.count("nnprime"), process_list.count("nTwoN"), process_list.count("nThreeN"),
#                         process_list.count("elastic"), process_list.count("np"), process_list.count("nNGamma"),
#                         process_list.count("other"), process_list.count("nnprimeP"),
#                         process_list.count("alphaKO"), process_list.count("nTwoNP"), process_list.count("nThreeNP"),
#                         process_list.count("nFourN")])
#
# for i in range(len(processes)):
#     print(processes[i], " count :", frequencies[i])
#
# cross_section(sum(frequencies), 10 ** 6, 0.1, string="total ")
# for i in range(len(frequencies)):
#     if frequencies[i] > 200:
#         cross_section(frequencies[i], 10 ** 6, 0.1, string=processes[i])
#
# n_p_idx = processes.index("40Ar(n,p)40*Cl")
# cross_section(frequencies[n_p_idx], 10 ** 6, 0.1, string="(n,p)")

# plt.hist(steps_length_dist, 40)
# plt.xlim([0, 400])
# plt.title("2.5MeV neutrons total path length in lAr")
# plt.ylabel("number of events")
# plt.xlabel("total path length [cm]")
# plt.show()
#
# xbins = np.arange(0, 3) - 0.5
# plt.xticks(np.arange(0, 2, 1))
# plt.hist(proton_num_dist, bins=xbins, rwidth=0.2)
# plt.title("number of protons ejected with $E_k$ > 0.5 MeV")
# plt.ylabel("number of events")
# plt.xlabel("number of protons")
# plt.show()
#
# xbins = np.arange(0, 3) - 0.5
# plt.xticks(np.arange(0, 2, 1))
# plt.hist(neutron_num_dist, bins=xbins, rwidth=0.2)
# plt.title("number of new neutrons ejected in an event")
# plt.ylabel("number of events")
# plt.xlabel("number of neutrons")
# plt.show()

# num_interactions_dist = dat6.numInteractions
# tot_path_length = dat6.sumStepLengths
#
# hist, xbins, ybins, im = two_d_plot(num_interactions_dist, tot_path_length,
#                                     "2.5 MeV neutrons total path length vs number of interactions",
#                                     "number of interactions", "total path length[cm]", 0.5, 7, 10, xlim=[-0.5, 9.5],
#                                     ylim=[-0.5, 204],
#                                     integer_xDist=True)

# print(photon_energy_dist[photon_energy_dist > 2.5])
# for energy in set(photon_energy_dist):
#     count = list(photon_energy_dist).count(energy)
#     if count > 100:
#         print("energy = ", energy, " MeV")
#         print("number of $\gamma$ = ", count, "\n")
# plt.xticks(np.arange(0, 6.2, 0.2))
# plt.hist(photon_energy_dist, bins=np.arange(0, 6.2, 0.05), rwidth=0.4)
# plt.title("$\gamma$ energy distribution across all events")
# plt.ylabel("number of $\gamma$")
# plt.xlabel("$E_{\gamma} [MeV]$ ")
# plt.yscale("log")
# plt.show()

# two_d_plot(distance_between_interactions_dist, time_between_interactions_dist,
#            "distance between interactions vs time between interactions 2D hist", "distance[cm]", "time[ns]")

# plt.hist(distance_between_interactions_dist, 50)
# plt.xlabel("distance[cm]")
# plt.ylabel("number of interactions")
# plt.title("distance between interactions histogram")
# plt.show()

# counts, bins, blabla = plt.hist(distance_between_interactions_dist, 100)
# plt.xlabel("distance between consequent interactions[cm]")
# plt.ylabel("count")
# plt.title("distance between interactions histogram")
#
# plt.plot(bins, max(counts) * np.e ** (-0.080228 * bins))
# plt.show()
#
# for i in range(len(bins) - 1):
#     bins[i] = (bins[i] + bins[i + 1]) / 2
#
# sig = np.sqrt(counts[:53]) / sum(counts[:53])
# counts = counts / sum(counts[:53])

# print(sig)
# print(counts)
# print(bins)

# fit = curve_fit(model_func, bins[1:54], counts[:53], 0.06, sig)
#
# print("\nestimated miu = ", fit[0][0])
# print("estimated std of miu = ", np.sqrt(fit[1][0][0]))
# print("chi square and p = ", chisquare(counts, max(counts) * np.e ** (-0.080228 * bins[1:])))

# plt.hist(tot_path_length,25)
# plt.xticks(np.arange(0, max(tot_path_length), 20))
# plt.xlabel('total path length[cm]')
# plt.ylabel('number of events')
# plt.title('2.5 MeV neutrons total path length in LAr')
# plt.show()
#
# plt.hist(neutronEnergy,50,rwidth=0.7)
# plt.xticks(np.arange(0, 2.6, 0.2))
# plt.xlabel('neutron escape energy[MeV]')
# plt.ylabel('count')
# plt.title('neutron escape energy distribution')
# plt.show()
#
# print(len(neutronEnergy))
# neutronEnergy = list(neutronEnergy)
# ids = list(ids)
# for i in range(len(ids)):
#     if ids[i] != i:
#         ids.insert(i, i)
#         neutronEnergy.insert(i, 0)

# plt.plot(tot_path_length, neutronEnergy, '.')
# plt.xticks(np.arange(0, max(tot_path_length) + 1, 15), fontsize=8)
# plt.xlabel('total path length[cm]')
# plt.ylabel('energy[MeV]')
# plt.title('neutron escape energy vs total path length')
# plt.show()

# two_d_plot(tot_path_length, neutronEnergy, 'neutron escape energy vs total path length', 'total path length[cm]',
#            'energy[MeV]', 6, 0.04, 8, xlim=[-0.5, 160])

# hist, xbins, ybins, im = plt.hist2d(distance_between_interactions_dist, time_between_interactions_dist, 50)
# plt.show()
#
# plt.plot(xbins, ybins, '-', label='average bin value slope')
# plt.plot(distance_between_interactions_dist, time_between_interactions_dist, '.', label='distance vs time distribution')
# plt.legend()
# plt.xlabel("distance[cm]")
# plt.ylabel("time[ns]")
# plt.title("distance vs time between interactions")
# plt.show()
#
# print(np.polyfit(xbins, ybins, 1))

# all_gammas_deposit = dat4.totEnergyDeposit
# secondary_gammas_deposit = dat4.secondaryGammasDeposit
#
# plt.hist(secondary_gammas_deposit, 50,align="left",rwidth=0.5)
# plt.xticks(np.arange(0, 8500, 250), fontsize=6)
# plt.xlabel('Energy Deposit[keV]')
# plt.ylabel('count')
# plt.yscale('log')
# plt.title('secondary $\gamma$ energy deposit per event')
# plt.show()

# neutronEnergy = np.array(neutronEnergy)
# neutronEnergytmp = neutronEnergy[num_interactions_dist > 1]
# print(len(neutronEnergy))
#
#
# h, bins = np.histogram(neutronEnergytmp, np.linspace(0,2.5,100))
# plt.step((bins[1:]+bins[:-1])/2, h,label='num interactions > 1')
#
# plt.xlabel('energy[MeV]')
# plt.ylabel('count')
# plt.title('n escape energy vs amount of n-nucleus interactions')
# neutronEnergy = neutronEnergy[num_interactions_dist == 1]
# h,bins = np.histogram(neutronEnergy, np.linspace(0,2.5,100))
# plt.step((bins[1:]+bins[:-1])/2, h, label='num interactions = 1')
# plt.plot([2.375 for i in range(100)],np.linspace(0,500,100),label='expected average neutron energy after one elastic scatter')
# plt.legend()
# plt.show()

# plt.hist(n_prime_energies, 50, rwidth=0.7)
# plt.xlabel("n' energy[MeV]")
# plt.ylabel("count")
# plt.title("n' energies after creation")
# plt.show()

# plt.hist(tot_path_length, 25)
# plt.xticks(np.arange(0, max(tot_path_length), 35))
# plt.xlabel('total path length[cm]')
# plt.ylabel('number of events')
# plt.title('2.5 MeV neutrons total path length in LAr')
# plt.show()
#
# plt.hist(neutronEnergy, 50, rwidth=0.7)
# plt.xticks(np.arange(0, 2.6, 0.2))
# plt.xlabel('neutron escape energy[MeV]')
# plt.ylabel('count')
# plt.title('neutron escape energy distribution')
# plt.show()

# plt.plot(tot_path_length, neutronEnergy, '.')
# plt.xticks(np.arange(0, max(tot_path_length) + 1, 25), fontsize=8)
# plt.xlabel('total path length[cm]')
# plt.ylabel('energy[MeV]')
# plt.title('neutron escape energy vs total path length')
# plt.show()
#
# two_d_plot(tot_path_length, neutronEnergy, 'neutron escape energy vs total path length', 'total path length[cm]',
#            'energy[MeV]', 10, 0.04, 10, xlim=[-0.5, 160])

# photon_energy_dist = photon_energy_dist[photon_energy_dist < 1.45]
# photon_energy_dist = photon_energy_dist[photon_energy_dist > 0.1]
# for energy in set(photon_energy_dist):
#     count = list(photon_energy_dist).count(energy)
#     if count > 20:
#         print("energy = ", energy, " MeV")
#         print("number of $\gamma$ = ", count, "\n")

# plt.xticks(np.arange(0, 6.2, 0.2))
# y,binEdges = np.histogram(photon_energy_dist,bins=np.arange(0, 6.2, 0.02))
# bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
# menStd     = np.sqrt(y)
# width      = 0.05
# plt.bar(bincenters, y, width=width, color='b', yerr=menStd)
# plt.title("$\gamma$ energy distribution across all events")
# plt.ylabel("number of $\gamma$")
# plt.xlabel("$E_{\gamma} [MeV]$ ")
# plt.yscale("log")
# plt.show()

# distance_between_interactions_dist = dat3.distanceDiff
# time_between_interactions_dist = dat3.timeDiff

# hist, xbins, ybins, im = plt.hist2d(distance_between_interactions_dist, time_between_interactions_dist, 50)
# plt.show()
#
# xs = np.linspace(0, 120, 100)
# ys1 = xs / 30
# ys2 = xs / (30 * 0.07281169519)
# plt.plot(xs, ys1, '-',color='yellow', label='speed of light')
# plt.plot(xs, ys2, '-',color='purple', label='speed of 2.5MeV neutrons')
# plt.plot(xbins, ybins, '-', label='average bin value slope')
# plt.plot(distance_between_interactions_dist, time_between_interactions_dist, '.', label='distance vs time distribution')
# plt.legend()
# plt.xlabel("distance[cm]")
# plt.ylabel("time[ns]")
# plt.title("distance vs time between interactions")
# plt.show()

# neutronEnergy12 = list(neutronEnergy12)
# ids12 = list(ids12)
# for i in range(len(ids12)):
#     if ids12[i] != i:
#         ids12.insert(i, i)
#         neutronEnergy12.insert(i, 0)
#
# neutronEnergyHalf = list(neutronEnergyHalf)
# idsHalf = list(idsHalf)
# for i in range(len(idsHalf)):
#     if idsHalf[i] != i:
#         idsHalf.insert(i, i)
#         neutronEnergyHalf.insert(i, 0)
#
# neutronEnergy1 = list(neutronEnergy1)
# ids1 = list(ids1)
# for i in range(len(ids1)):
#     if ids1[i] != i:
#         ids1.insert(i, i)
#         neutronEnergy1.insert(i, 0)
#
# h, bins = np.histogram(neutronEnergy1, np.linspace(0, 2.5, 100))
# plt.step((bins[1:] + bins[:-1]) / 2, h, label='detector size 1x1x1 meter')
#
# h, bins = np.histogram(neutronEnergy12, np.linspace(0, 2.5, 100))
# plt.step((bins[1:] + bins[:-1]) / 2, h, label='detector size 1x2x2 meter')
#
# h, bins = np.histogram(neutronEnergyHalf, np.linspace(0, 2.5, 100))
# plt.step((bins[1:] + bins[:-1]) / 2, h, label='detector size 0.5x0.5x0.5 meter')
#
# plt.legend()
# plt.xlabel('energy[MeV]')
# plt.ylabel('count')
# plt.title('n escape energy of different detector sizes')
# plt.show()

# numElasticInteractions = dat6.numElastic
#
# nPrimes = dat9.nPrime
# neutronEnergies = dat9.neutronEscapeEnergy
# ids1 = dat9.eventID
#
# nPrimeEnergies = neutronEnergies[nPrimes == 1]
# nPrimaryEnergies = neutronEnergies[nPrimes == 0]
#
# elasticEnergies = dat10.energy
#
# h, bins = np.histogram(neutronEnergies[numElasticInteractions == 1], np.linspace(0, 2.5, 400))
# plt.step((bins[1:] + bins[:-1]) / 2, h, label='neutrons escape energies which performed 1 elastic scatter')
#
# h, bins = np.histogram(neutronEnergies[numElasticInteractions == 2], np.linspace(0, 2.5, 400))
# plt.step((bins[1:] + bins[:-1]) / 2, h, label='neutrons escape energies which performed 2 elastic scatter')
#
# h, bins = np.histogram(neutronEnergies[numElasticInteractions == 3], np.linspace(0, 2.5, 400))
# plt.step((bins[1:] + bins[:-1]) / 2, h, label='neutrons escape energies which performed 3 elastic scatter')
#
# h, bins = np.histogram(elasticEnergies, np.linspace(0, 2.5, 800))
# plt.step((bins[1:] + bins[:-1]) / 2, h, label='n energies after performing one elastic scatter')
#
# h, bins = np.histogram(neutronEnergies, np.linspace(0, 2.5, 100))
# plt.step((bins[1:] + bins[:-1]) / 2, h, label='all neutrons escape energies')
#
# h, bins = np.histogram(nPrimeEnergies, np.linspace(0, 2.5, 100))
# plt.step((bins[1:] + bins[:-1]) / 2, h, label='neutrons which performed n nPrime escape energies')
#
# h, bins = np.histogram(nPrimaryEnergies, np.linspace(0, 2.5, 100))
# plt.step((bins[1:] + bins[:-1]) / 2, h, label='neutrons which didn't perform n nPrime escape energies')

# plt.legend()
# plt.xlabel('energy[MeV]')
# plt.ylabel('count')
# plt.title('n escape energies')
# plt.show()

# distance_between_interactions_dist = dat3.distanceDiff
#
# counts, bins, whatever = plt.hist(distance_between_interactions_dist, 50)
# xs = np.linspace(0, 120, 120)
# ys = max(counts) * np.e ** (-0.06772304510116327 * xs)
# plt.plot(xs, ys)
# plt.xlabel("distance[cm]")
# plt.ylabel("number of interactions")
# plt.title("distance between interactions histogram")
# plt.show()
#
# print(1 / 0.06772304510116327)
# print(0.0007965963367187171 / 0.06772304510116327 ** 2)

# def model(x, lamda, maxValue):
#     return maxValue * np.e ** (-x * lamda)
#
#
# xs = (bins[1:] + bins[:-1]) / 2
# ys = counts
# miu0 = 0.12
#
# for i in range(len(counts)):
#     counts[i] = float(counts[i])
#
# fit = curve_fit(model, xs, ys, [miu0, max(ys)], sigma=np.sqrt(ys))
# print("\nestimated miu = ", fit[0][0])
# print("estimated std of miu = ", np.sqrt(fit[1][0][0]))
# print(stats.ttest_ind(counts, max(counts) * np.e ** (-fit[0][0] * bins[1:])))

# recoil_angles = dat13.recoilAngle
# plt.hist(recoil_angles, 50)
# plt.xlabel("recoil angle [rad]")
# plt.ylabel("count")
# plt.title("n-Ar elastic scatter nucleus recoil angle distribution")
# plt.show()
#
# plt.hist(np.cos(recoil_angles), 50)
# plt.xlabel("cos (recoil angle)")
# plt.ylabel("count")
# plt.title("n-Ar elastic scatter nucleus recoil angle distribution")
# plt.show()

# formula_arr = [0.24012853645 * (np.cos(recoil_angles[i]) ** 2) for i in range(len(recoil_angles))]
# computed = np.array([2.5 - elasticEnergies[i] for i in
#                      range(len(elasticEnergies))])
# diff = abs(formula_arr - computed)
#
# plt.hist(diff,500)
# plt.xlabel("difference [MeV]")
# plt.ylabel("count")
# plt.title("difference between Geant4 calculated energy loss vs $\\frac{4Mm_n}{(M+m_n)^2}\\bullet cos^2(\\beta)$")
# plt.show()

# plt.hist(deposited_energy, 50, rwidth=0.5)
# plt.xticks(np.arange(0, max(deposited_energy), 200), fontsize=6)
# plt.xlabel('Energy Deposit[keV]')
# plt.ylabel('count')
# plt.yscale('log')
# plt.title('secondary $\gamma$ energy deposit per $\gamma$')
# plt.show()

# num_interactions_dist = dat6.numInteractions
# tot_path_length = dat6.sumStepLengths
#
# num_interactions_dist2 = dat.numInteractions
# tot_path_length2 = dat.sumStepLengths

# plt.hist(num_interactions_dist)
# plt.show()
# hist, xbins, ybins, im = two_d_plot(num_interactions_dist, tot_path_length,
#                                     "2.5 MeV neutrons total path length vs number of interactions",
#                                     "number of interactions", "total path length[cm]", 0.5, 7, 10, xlim=[-0.5, 9.5],
#                                     ylim=[-0.5, 204],
#                                     integer_xDist=True)

# CORNER_KWARGS = dict(
#     label_kwargs=dict(fontsize=15),
#     bins=[len(set(num_interactions_dist)), 25],
#     title_kwargs=dict(fontsize=16),
#     truth_color="tab:orange",
#     quantiles=[0.16, 0.84],
#     levels=(1 - np.exp(-0.5), 1 - np.exp(-2), 1 - np.exp(-9 / 2.0)),
#     plot_density=False,
#     plot_datapoints=False,
#     fill_contours=True,
#     max_n_ticks=20,
#     top_ticks=True,
#     verbose=False,
#     use_math_text=True,
# )
#
# data = pd.DataFrame(dict(numberOfInteractions=num_interactions_dist2, TotalPathLength=tot_path_length2))
# fig = corner.corner(data, **CORNER_KWARGS)
# axes = fig.get_axes()
# for i in range(len(axes)):
#     if i == 0:
#         axes[i].set_ylabel("number of events")
#         axes[i].set_yticks(np.arange(0, 3000, 400))
#         axes[i].set_yticklabels(np.arange(0, 3000, 400))
#     if i == 3:
#         axes[i].set_xlabel("total path length[cm]")
#         axes[i].set_yticks(np.arange(0, 4000, 400))
#         axes[i].set_yticklabels(np.arange(0, 4000, 400))
#         axes[i].set_ylabel("number of events")
#     axes[i].tick_params(axis='both', labelsize=6)
# fig.show()
# plt.plot(1, 2)
# plt.show()

# neutron_angle = dat17.neutronAngle
# h, bins = np.histogram(neutron_angle, np.linspace(0, np.pi, 50), density=True)
# # plt.step((bins[1:] + bins[:-1]) / 2, h, label="40Ar", color='purple')
# # plt.hist(neutron_angle, 50, label="40Ar")
# neutron_angle_H = dat18.neutronAngle
# h, bins = np.histogram(neutron_angle_H, np.linspace(0, np.pi, 50), density=True)
# plt.step((bins[1:] + bins[:-1]) / 2, h, label="${^{1}_{1}H}$", color='blue')
# # plt.hist(neutron_angle_H, 50, label="2H")
# neutron_angle_O = dat19.neutronAngle
# h, bins = np.histogram(neutron_angle_O, np.linspace(0, np.pi, 50), density=True)
# # plt.step((bins[1:] + bins[:-1]) / 2, h, label="16O", color='red')
# # plt.hist(neutron_angle_O, 50, label="16O")
# neutron_angle_Si = dat20.neutronAngle
# h, bins = np.histogram(neutron_angle_Si, np.linspace(0, np.pi, 50), density=True)
# # plt.step((bins[1:] + bins[:-1]) / 2, h, label="28Si", color="orange")
# # plt.hist(neutron_angle, 50, label="28Si")
# neutron_angle_Na = dat21.neutronAngle
# h, bins = np.histogram(neutron_angle_Na, np.linspace(0, np.pi, 50), density=True)
# # plt.step((bins[1:] + bins[:-1]) / 2, h, label="23Na", color='yellow')
# # plt.hist(neutron_angle, 50, label="23Na")
# neutron_angle_deut = dat23.neutronAngle
# h, bins = np.histogram(neutron_angle_deut, np.linspace(0, np.pi, 50), density=True)
# # plt.step((bins[1:] + bins[:-1]) / 2, h, label="2Deut", color='green')
# plt.legend()
# neutron_angle_Xe = dat24.neutronAngle
# h, bins = np.histogram(neutron_angle_Xe, np.linspace(0, np.pi, 50), density=True)
# # plt.step((bins[1:] + bins[:-1]) / 2, h, label="132Xe", color='black',alpha=0.2)
# plt.title(
#     "$ \Delta E = E_i \\bullet \\frac{4Mm_n}{(M+m_n)^2}\\bullet cos^2(\\beta)$ - neutron elastic scattering angle distribution for different materials ")
#
# plt.xlabel("neutron scattering angle[rad]")
# plt.ylabel("count (normalized)")
# # plt.yscale("log")
# plt.show()

# neutron_angle_H = dat18.neutronAngle
# cos_neutron = np.cos(neutron_angle_H)
# plt.hist(neutron_angle_H, 50)
# plt.show()
# plt.hist(cos_neutron)
# plt.show()
#
# neutron_angle_O = dat19.neutronAngle
# plt.hist(neutron_angle_O, 50)
# plt.title("oxygen")
# plt.show()
#
#
# plt.hist(np.cos(neutron_angle_O))
# plt.title("oxygen")
# plt.show()
#
# neutron_angle = dat20.neutronAngle
# plt.hist(neutron_angle, 50)
# plt.show()
#
# neutron_angle = np.cos(neutron_angle)
# plt.hist(neutron_angle)
# plt.show()
#
# neutron_angle = dat21.neutronAngle
# plt.hist(neutron_angle, 50)
# plt.title("Na")
# plt.show()
#
# neutron_angle = np.cos(neutron_angle)
# plt.hist(neutron_angle)
# plt.title("Na")
# plt.show()

# h, bins = np.histogram(neutronEnergies_1m, np.linspace(0, 2.5, 500))
# # plt.step((bins[1:] + bins[:-1]) / 2, h, label='detector size 1x1x1 meter')
# plt.errorbar((bins[1:] + bins[:-1]) / 2, h, fmt='.', yerr=np.sqrt(h), label='detector size 1x1x1 meter')
#
# h, bins = np.histogram(neutronEnergies_halfm, np.linspace(0, 2.5, 500))
# # plt.step((bins[1:] + bins[:-1]) / 2, h, label='detector size 1x2x2 meter')
# plt.errorbar((bins[1:] + bins[:-1]) / 2, h, fmt='.', yerr=np.sqrt(h), label='detector size 1x2x2 meter')
#
# h, bins = np.histogram(neutronEnergies_1m_2m, np.linspace(0, 2.5, 500))
# # plt.step((bins[1:] + bins[:-1]) / 2, h, label='detector size 0.5x0.5x0.5 meter')
# plt.errorbar((bins[1:] + bins[:-1]) / 2, h, fmt='.', yerr=np.sqrt(h), label='detector size 0.5x0.5x0.5 meter')

# plt.legend()
# plt.xlabel('energy[MeV]')
# plt.ylabel('count')
# plt.title('n escape energy of different detector sizes')
# plt.xticks(np.arange(0, 2.5, 0.1))
# plt.show()

# photon_energy = dat2.Energy
# emission_time = dat2.time
#
# #
# # print(photon_energy[photon_energy>6])
# print(emission_time)
# # h, bins = np.histogram(emission_time, 50)
# # plt.step((bins[1:] + bins[:-1]) / 2, h)
# plt.hist(emission_time, np.arange(0, 50, 0.2))
# plt.title("photon emission time")
# plt.xlabel("time[ns]")
# plt.ylabel("count")
# plt.yscale("log")
# plt.show()
# # plt.hist(photon_energy, 50)
# # plt.yscale("log")
# # plt.show()
# two_d_plot(photon_energy, emission_time, "photon energy vs emission time", "photon energy [MeV]", "emission time[ns]",
#            0.035, 0, 8, xlim=[0, 1.5])

# primary_flags = dat22.primary
# total_path_lengths = dat22.sumStepLengths
# secondary_n_tot_path = total_path_lengths[primary_flags == 0]
# primary_n_tot_path = total_path_lengths[primary_flags == 1]
#
# h, bins = np.histogram(primary_n_tot_path, 50)
# plt.step((bins[1:] + bins[:-1]) / 2, h, label="primary neutrons total path length")
#
# h, bins = np.histogram(secondary_n_tot_path, 50)
# plt.step((bins[1:] + bins[:-1]) / 2, h, label="secondary neutrons total path length")
#
# plt.xlabel("total path length[cm]")
# plt.ylabel("number of events")
# plt.title("neutrons total path length")
#
# # distance_until_inelastic = []
# # for i in range(len(total_path_lengths)-1):
# #     if primary_flags[i] == 1 and primary_flags[i + 1] == 0:
# #         distance_until_inelastic.append(total_path_lengths[i])
# #
# # print(min(distance_until_inelastic))
# # plt.hist(distance_until_inelastic, 50, label='distance until inelastic interaction')
# plt.legend()
# plt.show()
#
# primary_flags = dat5.nPrime
# primaries_escape_energies = dat5.neutronEscapeEnergy[primary_flags == 0]
# primaries_escape_energies = primaries_escape_energies[primaries_escape_energies != 0]
# secondary_escape_energies = dat5.neutronEscapeEnergy[primary_flags == 1]
# print(len(secondary_escape_energies))
# secondary_escape_energies = secondary_escape_energies[secondary_escape_energies != 0]
#
# h, bins = np.histogram(primaries_escape_energies, 100)
#
# # plt.bar((bins[1:] + bins[:-1]) / 2, h, width=0.02, label="primary neutrons escape energy",
# #         yerr=np.sqrt(h))
# plt.step((bins[1:] + bins[:-1]) / 2, h, label='primary neutrons escape energy')
#
# h, bins = np.histogram(secondary_escape_energies, 100)
# plt.step((bins[1:] + bins[:-1]) / 2, h, label='secondary neutrons escape energy')
# # plt.bar((bins[1:] + bins[:-1]) / 2, h, width=0.02,
# #         label="secondary neutrons escape energy", yerr=np.sqrt(h))
#
# plt.xlabel('escape energy[MeV]')
# plt.xticks(np.linspace(0, 2.5, 15))
# plt.ylabel('count')
# plt.title('neutron escape energy distribution')
# # plt.yscale('log')
# plt.legend()
# plt.show()
#
# plt.plot(primary_n_tot_path, primaries_escape_energies, '.', label='primary neutrons')
# plt.plot(secondary_n_tot_path, secondary_escape_energies, '.', label='secondary neutrons')
# plt.legend()
# plt.xlabel('total path length[cm]')
# plt.ylabel('neutron escape energy[MeV]')
# plt.title('neutron escape energy vs total path length')
# plt.show()

# primary_flags_1 = dat7.nPrime
# primary_neutronEnergy1 = dat7.neutronEscapeEnergy[primary_flags_1 == 0]
# primary_neutronEnergy1 = primary_neutronEnergy1[primary_neutronEnergy1 != 0]
# secondary_neutronEnergy1 = dat7.neutronEscapeEnergy[primary_flags_1 == 1]
# secondary_neutronEnergy1 = secondary_neutronEnergy1[secondary_neutronEnergy1 != 0]
#
# primary_flags_half = dat8.nPrime
# primary_neutronEnergy_half = dat8.neutronEscapeEnergy[primary_flags_half == 0]
# primary_neutronEnergy_half = primary_neutronEnergy_half[primary_neutronEnergy_half != 0]
# secondary_neutronEnergy_half = dat8.neutronEscapeEnergy[primary_flags_half == 1]
# secondary_neutronEnergy_half = secondary_neutronEnergy_half[secondary_neutronEnergy_half != 0]
#
# primary_flags_1_2 = dat5.nPrime
# primary_neutronEnergy1_2 = dat5.neutronEscapeEnergy[primary_flags_1_2 == 0]
# primary_neutronEnergy1_2 = primary_neutronEnergy1_2[primary_neutronEnergy1_2 != 0]
# secondary_neutronEnergy1_2 = dat5.neutronEscapeEnergy[primary_flags_1_2 == 1]
# secondary_neutronEnergy1_2 = secondary_neutronEnergy1_2[secondary_neutronEnergy1_2 != 0]
#
# h, bins = np.histogram(primary_neutronEnergy1, np.linspace(0, 2.5, 50))
# plt.step((bins[1:] + bins[:-1]) / 2, h, label='primary neutrons,detector size 1x1x1 meter', color='red')
# # plt.bar(bins[:-1], h, width=0.02, yerr=np.sqrt(h), label='primary neutrons,detector size 1x1x1 meter')
#
# h, bins = np.histogram(secondary_neutronEnergy1, np.linspace(0, 2.5, 50))
# plt.step((bins[1:] + bins[:-1]) / 2, h, '--', label='secondary neutrons,detector size 1x1x1 meter', color='red')
# # plt.bar(bins[:-1], h, width=0.02, yerr=np.sqrt(h), label='secondary neutrons,detector size 1x1x1 meter')
#
# h, bins = np.histogram(primary_neutronEnergy1_2, np.linspace(0, 2.5, 50))
# plt.step((bins[1:] + bins[:-1]) / 2, h, label='primary neutrons,detector size 1x2x2 meter', color='blue')
#
# # plt.bar(bins[:-1], h, width=0.02, yerr=np.sqrt(h), label='primary neutrons,detector size 1x2x2 meter')
#
# h, bins = np.histogram(secondary_neutronEnergy1_2, np.linspace(0, 2.5, 50))
# plt.step((bins[1:] + bins[:-1]) / 2, h, '--', label='secondary neutrons,detector size 1x2x2 meter', color='blue')
# # plt.bar(bins[:-1], h, width=0.02, yerr=np.sqrt(h), label='secondary neutrons,detector size 1x2x2 meter')
#
# h, bins = np.histogram(primary_neutronEnergy_half, np.linspace(0, 2.5, 50))
# plt.step((bins[1:] + bins[:-1]) / 2, h, label='primary neutrons,detector size 0.5x0.5x0.5 meter', color='black')
# # plt.bar(bins[:-1], h, width=0.02, yerr=np.sqrt(h), label='primary neutrons,detector size 0.5x0.5x0.5 meter')
#
# h, bins = np.histogram(secondary_neutronEnergy_half, np.linspace(0, 2.5, 50))
# plt.step((bins[1:] + bins[:-1]) / 2, h, '--', label='secondary neutrons,detector size 0.5x0.5x0.5 meter', color='black')
# # plt.bar(bins[:-1], h, width=0.02, yerr=np.sqrt(h), label='secondary neutrons,detector size 0.5x0.5x0.5 meter')
#
# plt.legend()
# plt.xlabel('energy[MeV]')
# plt.ylabel('count')
# plt.title('n escape energy of different detector sizes')
# plt.yscale('log')
# plt.show()

# filename = "neutronData/elasticScatters/elasticEnergies_neutron_angle.csv"
# dat = pd.read_csv(filename)
#
# elasticEnergies = dat.energy
# recoil_angles = dat.recoilAngle
# nucleus_energies = 2.5 - elasticEnergies
#
# filename2 = "inelasticEnergies_1m.csv"
# dat2 = pd.read_csv(filename2)
#
# inelastic_recoil_angles = dat2.recoilAngle
# inelastic_nucleus_energies = dat2.energy
#
# filename3 = "nNPrime_scat_ang_dist.csv"
# dat3 = pd.read_csv(filename3)
# emitted_n_angle = dat3.neutronAngle
#
# plt.hist(np.cos(emitted_n_angle), 50)
# plt.ylabel('count')
# plt.title("cosine of emitted n' angle following (n,n') process")
# plt.show()
# # h, bins = np.histogram(nucleus_energies * 1000, 50)
# # plt.step((bins[1:] + bins[:-1]) / 2, h, label='following 2.5MeV neutron elastic scatter')
# #
# # h, bins = np.histogram(inelastic_nucleus_energies, 100)
# # plt.step((bins[1:] + bins[:-1]) / 2, h, label='following (n,nPrime) process')
# #
# # plt.legend()
# # plt.xlabel("energy[keV]")
# # plt.ylabel("count")
# # plt.title("${^{40}}Ar$ nucleus recoil energies")
# # plt.show()
# #
# plt.hist(inelastic_recoil_angles, 50)
# plt.xlabel("recoil angle[rad]")
# plt.ylabel("count")
# plt.title("${^{40}}Ar$ nucleus recoil angles following (n,n') process")
# plt.show()
#
# plt.plot(inelastic_nucleus_energies, inelastic_recoil_angles, '.')
# plt.xlabel("energy[keV]")
# plt.ylabel("recoil angle[rad]")
# plt.title("${^{40}}Ar$ nucleus recoil energies vs recoil angle following (n,n') process")
# plt.show()

# filename2 = "inelasticEnergies_1m.csv"
# dat2 = pd.read_csv(filename2)
# inelastic_recoil_angles = dat2.recoilAngle
#
# plt.hist(inelastic_recoil_angles, 50)
# plt.xlabel("recoil angle[rad]")
# plt.ylabel("count")
# plt.title("${^{40}}Ar$ nucleus recoil angles following (n,n') process")
# plt.show()

#####################################################################################
# Creus et al. paper simulation code

# dat = pd.read_csv("C:\\Users\\Aviv\\PycharmProjects\\pythonProject3\\NuclearRecoils\\25\elasticEnergies.csv")
# dat2 = pd.read_csv("C:\\Users\\Aviv\\PycharmProjects\\pythonProject3\\NuclearRecoils\\25\inelasticEnergies.csv")
#
# # nucleus_energy_experimental = dat3.experimentalRecoilEnergy
# # neutron_time_of_flight_experimental = dat3.timeOfFlight
# # nucleus_energy_experimental = nucleus_energy_experimental[neutron_time_of_flight_experimental > 42]
# # nucleus_energy_experimental = nucleus_energy_experimental[neutron_time_of_flight_experimental < 49]
#
# # print("total number of detected events = ", len(nucleus_energy_experimental))
#
# # plt.hist(nucleus_energy_experimental, 50)
# # plt.xlabel('recoil energy[keV]')
# # plt.ylabel('count')
# # plt.title('experimental nucleus recoil energy reconstruction')
# # plt.show()
# #
# # plt.hist(neutron_time_of_flight_experimental, 50)
# # plt.xlabel('time[ns]')
# # plt.ylabel('count')
# # plt.title('experimental neutron time of flight')
# # plt.show()
#
# nucleus_energy = dat.energy
# neutron_time_of_flight_elastic = dat.timeOfFlight
# nucleus_energy = nucleus_energy[neutron_time_of_flight_elastic > 42]
# nucleus_energy = nucleus_energy[neutron_time_of_flight_elastic < 49]
# nucleus_energy_experimental = [11.5 for _ in range(len(nucleus_energy))]
# numElastic = dat.numElastic
# nucleus_energy_elastic = nucleus_energy[numElastic == 1]
# nucleus_energy_multiple = nucleus_energy[numElastic > 1]
# neutron_angle_elastic = dat.neutronAngle
# external = dat.external
# externals = nucleus_energy[external == 1]
# singles = nucleus_energy_elastic[external == 0]
# print("num of purely external events = ", len(externals[numElastic == 0]))
# print("num of single elastic events = ", len(singles))
# print("num of single elastic+external = ", len(nucleus_energy_elastic))
# print("num of multiple elastics = ", len(nucleus_energy_multiple[external == 0]))
# print("total number of detected events = ",
#       len(externals[numElastic == 0]) + len(nucleus_energy_elastic) + len(nucleus_energy_multiple))
# event_ids_elastic = dat.eventID[neutron_time_of_flight_elastic > 42]
# event_ids_elastic = event_ids_elastic[neutron_time_of_flight_elastic < 49]
#
# nucleus_energy_inelastic = dat2.energy
# neutron_time_of_flight_inelastic = dat2.timeOfFlight
# event_ids_inelastic = dat2.eventID
#
# # time of flight cutoff
# nucleus_energy_inelastic = nucleus_energy_inelastic[neutron_time_of_flight_inelastic > 42]
# nucleus_energy_inelastic = nucleus_energy_inelastic[neutron_time_of_flight_inelastic < 49]
# print("num of inelastic events = ", len(nucleus_energy_inelastic))
# # print(len(nucleus_energy_experimental) - len(nucleus_energy_elastic) - len(nucleus_energy_inelastic))
#
# plt.hist(neutron_angle_elastic[external == 0] * 180 / np.pi, 200)
# plt.xlabel('neutron scattering angle[degrees]')
# plt.ylabel('count')
# plt.title('neutron scattering angle distribution')
# plt.show()
#
# h, bins = np.histogram(nucleus_energy_elastic[external == 0], np.linspace(0, 100, 200))
# plt.step((bins[1:] + bins[:-1]) / 2, h, label='single elastic', color='blue')
#
# h1, bins = np.histogram(nucleus_energy_elastic, np.linspace(0, 100, 200))
# plt.errorbar((bins[1:] + bins[:-1]) / 2, h1, yerr=np.sqrt(h1), fmt='.', label='single elastic + external', color='red')
#
# h, bins = np.histogram((nucleus_energy_multiple, np.linspace(0, 100, len(nucleus_energy_multiple))))
# plt.errorbar((bins[1:] + bins[:-1]) / 2, h, yerr=np.sqrt(h), fmt='.', label='multiple elastic', color='green')
#
# xs = [11.5 for i in range(500)]
# ys = np.linspace(0, max(h1), 500)
#
# h2, bins = np.histogram(nucleus_energy_inelastic, np.linspace(0, 100, 200))
# plt.errorbar((bins[1:] + bins[:-1]) / 2, h2, yerr=np.sqrt(h2), fmt='.', label='inelastic', color='orange')
#
# plt.plot(xs, ys, label='theoretical recoil energy = 11.5 keV', color='black')
# # plt.plot((bins[1:] + bins[:-1]) / 2, h1 + h2, '.', label='total', color='black')
#
# plt.xlabel('energy[keV]')
# plt.ylabel('count')
# plt.yscale('log')
# plt.title('recoil nucleus energy distribution following elastic scatter with $\\theta = 25^{\circ}$')
# plt.legend()
# plt.ylim([1, 4000])
# plt.show()
#
# # print(neutron_time_of_flight_elastic[neutron_time_of_flight_elastic < 0])
# neutron_time_of_flight_elastic = neutron_time_of_flight_elastic[neutron_time_of_flight_elastic > 0]
# neutron_time_of_flight_single = neutron_time_of_flight_elastic[numElastic == 1]
#
# h, bins = np.histogram(neutron_time_of_flight_single[external == 0], np.linspace(0, 150, 300))
# plt.step((bins[1:] + bins[:-1]) / 2, h, label='single elastic', color='blue')
#
# h1, bins = np.histogram(neutron_time_of_flight_single, np.linspace(0, 150, 300))
# plt.errorbar((bins[1:] + bins[:-1]) / 2, h1, yerr=np.sqrt(h1), fmt='.', label='single elastic + external', color='red')
#
# h, bins = np.histogram(
#     neutron_time_of_flight_elastic[numElastic > 1], 100)
# plt.errorbar((bins[1:] + bins[:-1]) / 2, h, yerr=np.sqrt(h), fmt='.', label='multiple elastic', color='green')
#
# xs = [43 for _ in range(500)]
# ys = np.linspace(0, max(h), 500)
# plt.plot(xs, ys, label='paper mean time of flight for single elastic', color='black')
#
# h2, bins = np.histogram(neutron_time_of_flight_inelastic, np.linspace(0, 150, 300))
# plt.errorbar((bins[1:] + bins[:-1]) / 2, h2, yerr=np.sqrt(h2), fmt='.', label='inelastic', color='orange')
#
# # plt.plot((bins[1:] + bins[:-1]) / 2, h1 + h2, '.', label='total', color='black')
#
# plt.xlabel('time of flight[ns]')
# plt.ylabel('count')
# plt.title('neutron time of flight histogram')
# plt.yscale('log')
# plt.legend()
# xs = np.linspace(42, 49, 100)
# ys1 = [0 for _ in range(100)]
# ys2 = [max(h1) for _ in range(100)]
#
# plt.fill_between(xs, ys1, ys2, color='green')
# # plt.xlim([41, 48])
# plt.show()
#
# nucleus_energy_experimental = list(nucleus_energy_experimental)
# nucleus_energy = list(nucleus_energy)
# diff = []
# for i in range(len(nucleus_energy)):
#     diff.append((nucleus_energy_experimental[i] - nucleus_energy[i]) / nucleus_energy[i])
#
# ratios = [nucleus_energy_experimental[i] / nucleus_energy[i] for i in range(len(nucleus_energy_experimental))]
#
# plt.plot(singles, 11.5 / singles, '.', label='single elastic', color='blue')
# # plt.plot(nucleus_energy_elastic[external == 1], 11.5 / nucleus_energy_elastic[external == 1],
# #          '.', label='single elastic+external', color='red')
# # plt.plot(nucleus_energy_multiple, 11.5 / nucleus_energy_multiple, '.', label='multiple elastic', color='green')
# # plt.plot(nucleus_energy, ratios, '.', label='values')
# # plt.plot(nucleus_energy, 11.5 / np.array(nucleus_energy), '.', label='fit',color='magenta')
# plt.xlabel('$E_{True}[keV]$')
# plt.ylabel('$\\frac{E_{experimental}}{E_{True}}$')
# plt.title('$\\frac{E_{experimental}}{E_{True}}$ vs $E_{True}$')
# plt.yscale('log')
# plt.xticks([0, 11.5, 50, 100, 200, 400])
# plt.legend()
# plt.xlim([6, 17])
# # plt.ylim([0, 5])
# plt.show()
#
#
# def model(x, a):
#     return a / x
#
#
# fit = curve_fit(model, nucleus_energy, ratios, 10, sigma=np.sqrt(ratios))
# print("\nestimated a = ", fit[0][0])
# print("estimated std of a = ", np.sqrt(fit[1][0][0]))
# nucleus_energy = np.array(nucleus_energy)
#
#
# def chisquared(f_obs, f_exp):
#     return sum((f_obs[i] - f_exp[i]) ** 2 / f_obs) / len(f_obs)
#
#
# print("chi square = ", chisquared(ratios, 11.5 / np.array(nucleus_energy)))
# print(ttest_ind(ratios, fit[0][0] * np.sqrt(1 / nucleus_energy)))
#
# plt.plot(singles, (11.5 - singles) / singles, '.', label='single elastic', color='blue')
# # plt.plot(nucleus_energy_elastic[external == 1],
# #          (11.5 - nucleus_energy_elastic[external == 1]) / nucleus_energy_elastic[external == 1],
# #          '.', label='single elastic+external', color='red')
# # plt.plot(nucleus_energy_multiple, (11.5 - nucleus_energy_multiple) / nucleus_energy_multiple, '.',
# #          label='multiple elastic', color='green')
# plt.xlabel('$E_{True}[keV]$')
# plt.ylabel('$\\frac{E_{experimental} - E_{True}}{E_{True}}$')
# plt.title('$\\frac{E_{experimental} - E_{True}}{E_{True}}$ vs $E_{True}$')
# plt.legend()
# # plt.yscale('log')
# plt.xlim([6, 17])
# # plt.ylim([-1, 5])
# plt.show()
#
# plt.hist(ratios, np.linspace(0, 20, 200))
# plt.xticks([i for i in range(20)])
# plt.yscale('log')
# plt.xlabel('$\\frac{E_{experimental}}{E_{True}}$')
# plt.ylabel('count')
# plt.title('ratio of reconstructed recoil energy and true recoil energy per event')
# plt.show()
#
# # plt.plot(event_ids_elastic, ratios, '.')
# # plt.yscale('log')
# # plt.xlabel('eventID')
# # plt.ylabel('$\\frac{E_{experimental}}{E_{True}}$')
# # plt.title('ratio of reconstructed recoil energy and true recoil energy per event')
# # plt.show()
#
# my_singles = [0.79, 0.78, 0.66, 0.602, 0.503, 0.497]
# my_single_ext = [0.12, 0.11, 0.16, 0.177, 0.239, 0.223]
# my_mult = [0.09, 0.11, 0.18, 0.22, 0.258, 0.28]
#
# their_singles = [0.677, 0.653, 0.548, 0.497, 0.409, 0.291]
# their_single_ext = [0.235, 0.251, 0.314, 0.343, 0.318, 0.451]
# their_mult = [0.088, 0.096, 0.137, 0.160, 0.205, 0.256]
#
# angles = [25, 30, 40, 50, 60, 90]
#
# plt.plot(angles, my_singles, '-', color='blue')
# plt.plot(angles, my_single_ext, '-', color='red')
# plt.plot(angles, my_mult, '-', color='green')
#
# plt.plot(angles, their_singles, '--', color='blue')
# plt.plot(angles, their_single_ext, '--', color='red')
# plt.plot(angles, their_mult, '--', color='green')
#
# plt.xlabel('angle[degrees]')
# plt.ylabel('percentage of contribution')
# plt.title('Fractional contributions to the recoil spectrum after time-of-flight cut')
# plt.ylim([0,1])
# plt.show()

#########################################################################################

# plt.hist(tot_path_length, 100)
# plt.title("primary e- total path length in LAr")
# plt.xlabel('path length[cm]')
# plt.ylabel('count')
# plt.show()

# dat2 = pd.read_csv("C:\\Users\\Aviv\\PycharmProjects\\pythonProject3\\scintCheck_e-.csv")
# energies = dat2.energy
# times = dat2.time
# thetas = dat2.theta
# phis = dat2.phi
# xpos = dat2.xpos
# ypos = dat2.ypos
# zpos = dat2.zpos
#
# dat3 = pd.read_csv("C:\\Users\\Aviv\\PycharmProjects\\pythonProject3\\primaryPos.csv")
# ppx = dat3.xpos
# ppy = dat3.ypos
# ppz = dat3.zpos
#
# plt.hist(energies, 100)
# plt.title("scintillation photons energy distribution")
# plt.xlabel("energy[eV]")
# plt.ylabel("count")
# plt.show()
#
# plt.hist(times, 5000)
# # plt.xlim([0, 2000])
# plt.title("scintillation photons emission time distribution")
# plt.xlabel("time[ns]")
# plt.ylabel("count")
# plt.yscale('log')
# plt.show()
#
# x1 = len(times[times < 40]) / len(times)
# x2 = len(times[times > 40]) / len(times)
# txt1 = "percentage of scintillation photons emitted before 40ns = {percent:.2f}"
# txt2 = "percentage of scintillation photons emitted after 40ns = {percent:.2f}"
# print(txt1.format(percent=x1))
# print(txt2.format(percent=x2))
# #
# plt.hist(thetas, 100)
# plt.title("scintillation photons theta distribution")
# plt.xlabel("angle[rad]")
# plt.ylabel("count")
# plt.show()
#
# plt.hist(phis, 100)
# plt.title("scintillation photons phi distribution")
# plt.xlabel("angle[rad]")
# plt.ylabel("count")
# plt.show()
#
# plt.plot(ypos, zpos, '.', label='scint photons vertex pos')
# plt.plot(ppy, ppz, '.', label='primary electron pos')
# plt.legend()
# plt.xlabel('ypos[cm]')
# plt.ylabel('zpos[cm]')
# plt.title('scintillation photons ypos vs zpos')
# plt.show()
#
# plt.plot(xpos, ypos, '.', label='scint photons vertex pos')
# plt.plot(ppx, ppy, '.', label='primary electron pos')
# plt.legend()
# plt.xlabel('xpos[cm]')
# plt.ylabel('ypos[cm]')
# plt.title('scintillation photons xpos vs ypos')
# plt.show()
#
# plt.plot(xpos, zpos, '.')
# plt.xlabel('xpos[cm]')
# plt.ylabel('zpos[cm]')
# plt.title('scintillation photons xpos vs zpos')
# plt.show()

# fig = plt.figure()
# ax = fig.add_subplot(projection='3d')
# ax.scatter(ppx, ppy, ppz, label='primary electron pos')
# ax.scatter(xpos, ypos, zpos, label='scint photons vertex pos')
# ax.legend()
# ax.set_xlabel('xpos[cm]')
# ax.set_ylabel('ypos[cm]')
# ax.set_zlabel('zpos[cm]')
# ax.set_title('primary electron position over the course of the event')
# plt.show()

# dat = pd.read_csv("C:\\Users\\Aviv\\PycharmProjects\\pythonProject3\\scintDist_smallBirks.csv")
# n_scint = dat.nScint
# n_photons = dat.numPhotons
# # tot_path_length = dat.sumStepLength
#
# plt.hist(n_photons[n_photons > 50000], 50, align='right')
# plt.title("number of scintillation photons created per event with 1MeV electron primary and small birks constant")
# plt.xlabel("number of photons")
# plt.ylabel("count")
# plt.show()

# datt = pd.read_csv("C:\\Users\\Aviv\\PycharmProjects\\pythonProject3\\PMTS_p.csv")
# nTot = datt.numPhotons
# nUp = datt.numUp
# nDown = datt.numDown
#
# plt.hist(nTot, 200)
# plt.xlabel('number of photons')
# plt.ylabel('count')
# plt.title('number of scintillation photons created per event with 1MeV proton primary')
# plt.show()
#
# plt.hist(100 * nUp / nTot, 100, label='nUp')
# plt.hist(100 * nDown / nTot, 100, label='nDown')
# plt.legend()
# plt.xlabel('percentage of photons that reached the PMTS')
# plt.ylabel('count')
# plt.show()
#
# dat2 = pd.read_csv("C:\\Users\\Aviv\\PycharmProjects\\pythonProject3\\scintCheck_p.csv")
# energies = dat2.energy
# times = dat2.time
# thetas = dat2.theta
# phis = dat2.phi
# xpos = dat2.xpos
# ypos = dat2.ypos
# zpos = dat2.zpos

# dat = pd.read_csv("C:\\Users\\Aviv\\PycharmProjects\\pythonProject3\\NuclearRecoils\\90\elasticEnergies.csv")
#
# nucleus_energy = dat.energy
# neutron_time_of_flight_experimental = dat.timeOfFlight
# nucleus_energy = nucleus_energy[neutron_time_of_flight_experimental > 42]
# nucleus_energy = nucleus_energy[neutron_time_of_flight_experimental < 49]
#
# print("total number of detected events = ", len(nucleus_energy))
#
# # plt.hist(nucleus_energy, 500)
# h1, bins1 = np.histogram(nucleus_energy, np.linspace(0, 200, 100))
# bins1_centers = (bins1[1:] + bins1[:-1]) / 2
#
# parameters = curve_fit(gauss, bins1_centers, h1,
#                        p0=[5, 31, 120, 7])
# miu = parameters[0][2]
# d_miu = parameters[1][2][2]
# sigma = parameters[0][3]
# d_sigma = parameters[1][3][3]
#
# bools = bins1_centers > 110
# bools = bools < 130
# bins1_centers = bins1_centers[bools]
# h1 = h1[bools]
# plt.step(bins1_centers, h1, where='mid')
# plt.plot(bins1_centers, gauss(bins1_centers, parameters[0][0], parameters[0][1], miu, sigma),
#          '--',
#          label='Fit result', )
# plt.xlabel('recoil energy[keV]')
# plt.ylabel('count')
# plt.title('nucleus recoil energy reconstruction')
# plt.show()
#
# print("gaussian H + A * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))\n with params H = ", parameters[0][0], ", A = ",
#       parameters[0][1], ", $\mu$ = ", miu, ", $\sigma$ = ", sigma)
#
# print("miu = ", miu, " +- ", np.sqrt(sigma ** 2 + d_sigma ** 2 + d_miu ** 2))
#
# print("chisquare = ",
#       chisquare(h1, gauss(bins1_centers, parameters[0][0], parameters[0][1], parameters[0][2], parameters[0][3]))[
#           0] / len(h1))
# print("pvalue = ",
#       chisquare(h1, gauss(bins1_centers, parameters[0][0], parameters[0][1], parameters[0][2], parameters[0][3]))[1])

# xs = [11.5, 16.4, 28.5, 43.4, 60.5, 119.5]
# dx1 = [2.1, 2.4, 3.2, 3.8, 4.3, 4.9]
# ys1 = [0.321, 0.276, 0.209, 0.174, 0.151, 0.1112]
# dy1 = [0.093, 0.075, 0.055, 0.043, 0.035, 0.022]
# dx2 = [2.8, 3.9, 5.2, 6.0, 7.1, 7.7]
# ys2 = [0.386, 0.305, 0.285, 0.294, 0.283, 0.301]
# dy2 = [0.032, 0.020, 0.014, 0.023, 0.018, 0.019]
# ys3 = [11.3235 / 159.1653, 13.8651 / 159.1653, 18.309 / 159.1653, 23.108 / 159.1653, 28.097 / 159.1653,
#        41.344 / 159.1653]
# dy3 = [0.021, 0.025, 0.032, 0.038, 0.044, 0.055]
# plt.errorbar(xs, ys1, xerr=dx1, yerr=dy1, color='b', ecolor='b', fmt='-',
#              label='my data,normalized vs equivalent er energy')
# plt.errorbar(xs, ys2, xerr=dx2, yerr=dy2, color='grey', ecolor='grey', fmt='-', label='Creus et al.')
# plt.errorbar(xs, ys3, xerr=dx1, yerr=dy1, color='b', ecolor='b', fmt='--',
#              label='my data, all normalized vs 60keV pp')
# plt.xlabel('nuclear recoil energy[keV]')
# plt.ylabel('L_eff')
# plt.title('$L_{eff}$ as a function of nuclear recoil energy')
# plt.legend()
# plt.show()

# filename = "../references/Default_Dataset_diff.csv"
# points = pd.read_csv(filename)
# xs = points.x
# ys = points.y
# #
# filename2 = "../references/Default_Dataset2.csv"
# points2 = pd.read_csv(filename2)
# xs2 = points2.x
# ys2 = points2.y
#
# for i in range(len(xs2)):
#     xs2[i] *= 1000
#
# # in keV
# energies = [10, 20, 50, 75, 100, 250, 500, 1000]
# mius = [5.098, 0.7053, 0.186, 0.156, 0.143, 0.1098, 0.0830, 0.0627]
# stds = np.array([0.00583, 0.00194, 0.00197, 0.00160, 0.00299, 0.000598, 0.00148, 0.00040])
#
# mius = [5.104, 0.733, 0.202, 0.171, 0.157, 0.117, 0.0916, 0.0672]
# stds = np.array([0.00643, 0.00195, 0.00113, 0.00134, 0.000823, 0.00077, 0.000409, 0.000123])
#
# # plt.plot(xs, ys, '--', color='orange',
# #          label='$\mu = \mu_\\tau + \mu_\sigma + \mu_{\sigma R} + \mu_\kappa $ [Shehata et al. ,ME Conf.Proc.4(2007)]')
# plt.errorbar(energies, mius, color='orange', yerr=stds, fmt='.', capthick=2, capsize=2, markeredgecolor='k',
#              label='VALIDATION - Geant4 simulated $\gamma$ attenuation coefficient in water')
# plt.plot(xs2, ys2, '--', color='orange', label='Grupen, Tum. Ther. Part. Beams. 538.(2000)')
# energies_lAr = np.array([10, 20, 50, 75, 100, 250, 500, 1000, 10000, 50000, 100000, 1000000])
# energies_lAr = energies_lAr
# mius_lAr = [87.154, 11.765, 0.862, 0.372, 0.247, 0.136, 0.103, 0.0761, 0.0333, 0.0375, 0.0423, 0.0533]
# stds_lAr = [0.0735, 0.0172, 0.00239, 0.00370, 0.00132,
#             0.000706, 0.000417, 0.000220, 8.25e-05, 0.000171, 0.000183, 0.000638]
# plt.errorbar(energies_lAr, mius_lAr, yerr=stds_lAr, fmt='.', capthick=2, capsize=2, markersize=10, markeredgecolor='k',
#              color='blue', label='Geant4 simulated $\gamma$ attenuation coefficient in lAr')
# cs = CubicSpline(energies_lAr, mius_lAr)
# plt.plot(energies_lAr, cs(energies_lAr), '--')
# # plt.legend(loc='upper right', fontsize=9)
# plt.grid()
# plt.xscale("log")
# plt.yscale("log")
# plt.xticks(fontsize=20)
# plt.yticks(fontsize=20)
# plt.xlabel("Photon energy [keV]",fontsize=20)
# plt.ylabel("Linear attenuation coefficient [1/cm]",fontsize=20)
# # plt.title("Linear attenuation coefficient vs photon energy")
# plt.show()

# dat = pd.read_csv('../edep_3500.csv')
#
# edep = dat.eDep
# h, bins = np.histogram(edep, 100)
# plt.step((bins[1:] + bins[:-1]) / 2, h)
# plt.xlabel('energy deposit [keV]')
# plt.ylabel('counts')
# plt.title('energy deposit in liquid_scintillator histogram, 3.5 MeV primary neutron')
# plt.show()

# dat = pd.read_csv("C:\\Users\\Aviv\\PycharmProjects\\pythonProject3\\NuclearRecoils\\25\\truth_energies_25_wide.csv")
#
# nucleus_energy_elastic = dat.experimentalRecoilEnergy
# h1, bins = np.histogram(nucleus_energy_elastic, np.arange(0, max(nucleus_energy_elastic), 1))
# plt.errorbar((bins[1:] + bins[:-1]) / 2, h1, yerr=np.sqrt(h1), fmt='.', label='single elastic + external', color='red')
#
# plt.xlabel('energy[keV]')
# plt.ylabel('count')
# plt.yscale('log')
# plt.title('recoil nucleus energy distribution following elastic scatter with $\\theta = 25^{\circ}$')
# plt.legend()
# plt.ylim([1, 4000])
# plt.show()

# plt.hist(coincidence_time, 100)
# plt.xlabel('event coincidence time between PMTs[ns]')
# plt.ylabel('number of events')
# plt.title('coincidence time between PMTs per event')
# plt.show()

# plot_pmts("../Scintillation_Data/electron_recoil_scint\\PMTS_60.csv", 60, 1, 120, [0, 100, 158, 10])
#
# plot_pmts("../Scintillation_Data/electron_recoil_scint\\PMTS_122.csv", 122, 1, 200, [80, 240, 330, 15])
#
# plot_pmts("../Scintillation_Data/electron_recoil_scint\\PMTS_511.csv", 511, 0, 200, [])
#
# plot_pmts("../Scintillation_Data/electron_recoil_scint\\PMTS_662.csv", 662, 0, 200, [])
#
# plot_pmts("../Scintillation_Data/electron_recoil_scint\\PMTS_1274.csv", 1274, 0, 200, [])

# filenames = ["edep_500_new.csv", "edep_1000_new.csv", "edep_1500_new.csv", "edep_2000_new.csv", "edep_2500_new.csv",
#              "edep_3000_new.csv", "edep_3500_new.csv"]
# energies = [0.5, 1, 1.5, 2, 2.5, 3, 3.5]
# percentages = []
#
# for file in filenames:
#     percentages.append(percentage_detected(file))
#
# plt.plot(energies, percentages, '-')
# plt.xlabel('neutron energy')
# plt.ylabel('% detected')
# plt.title('percentage of detected events as a function of initial neutron energy')
# plt.show()

# for R in Rs:
#     for L_eff in L_effs:
#         # creating new arrays which are smeared according to gaussian dist with sigma around original array
#         num_up_new = []
#         num_down_new = []
#         for k in range(len(num_up)):
#             num_up_new.append(np.random.normal(num_up[k], R * np.sqrt(num_up[k])))
#             num_down_new.append(np.random.normal(num_down[k], R * np.sqrt(num_down[k])))
#         num_up_new = np.array(num_up_new)
#         num_down_new = np.array(num_down_new)
#
#         # convolution with a gaussian to account for the gain fluctuation of the PMTs
#         # gaussian = [np.exp(-(num_up_new[m] / (0.4 * num_up_new[m])) ** 2 / 2) for m in range(len(num_up_new))]
#         # num_up_new = np.convolve(num_up_new, gaussian, mode="full")
#         # num_down_new = np.convolve(num_down_new, gaussian, mode="full")
#
#         h1, bins1 = np.histogram(num_up_new * L_eff, np.arange(0, max(num_up), 2))
#         # plt.step((bins1[1:] + bins1[:-1]) / 2, h1, where='mid', label='number of photons absorbed in top PMT')
#
#         h2, bins2 = np.histogram(num_down_new * L_eff, np.arange(0, max(num_up), 2))
#         # plt.step((bins2[1:] + bins2[:-1]) / 2, h2, where='mid', label='number of photons absorbed in bottom PMT')
#         #
#         # plt.xlabel('number of photons')
#         # plt.ylabel('number of events')
#         # plt.title('number of photons absorbed in bottom and top PMT per event')
#         # plt.legend(loc='upper right', fontsize=8)
#         # plt.show()
#
#         h = (h1 + h2) / 2
#         # bin size should be roughly 2
#         I_katom = sum(counts[i] * (num_photons[i + 1] - num_photons[i]) for i in range(len(counts) - 1)) + 2 * counts[
#             -1]
#         # I_katom = sum(counts) * 2
#         h = h / (2 * sum(h)) * I_katom
#         # print(sum(h))
#         # print(sum(counts))
#         chi_square = sum(((counts[i] - h[i]) ** 2) / (errs[i] ** 2 + h[i]) for i in range(len(counts)))
#         chi_squares[i][j] = chi_square
#         if chi_square < min_chi_square:
#             min_chi_square = chi_square
#             min_h = h
#             min_bins = bins2
#             best_scale = L_eff
#             best_sigma = R
#         # counts = counts / sum(counts)
#         i += 1
#     j += 1
#     i = 0

# fit with curve_fit
# bins1_centers = (bins1[1:] + bins1[:-1]) / 2
#
# start = 1
# end = 21
# parameters = curve_fit(gauss, bins1_centers[start:end], h[start:end],
#                        p0=[4, 25, 22, 6])
# light_yield_n = parameters[0][2]
# err_n = parameters[0][3]
# print("gaussian H + A * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))\n with params H = ", parameters[0][0], ", A = ",
#       parameters[0][1], ", $\mu$ = ", parameters[0][2], ", $\sigma$ = ", parameters[0][3])
#
# xs = bins1_centers[start:end]
# plt.plot(xs, gauss(xs, parameters[0][0], parameters[0][1], parameters[0][2], parameters[0][3]), '--',
#          label='Fit result', )
# plt.step(bins1_centers, h, where='mid', label='number of photons absorbed in PMT')
# plt.legend()
# plt.show()
#
# print("chisquare = ",
#       chisquare(h[start:end], gauss(xs, parameters[0][0], parameters[0][1], parameters[0][2], parameters[0][3]))[
#           0] / len(h))
# print("pvalue = ",
#       chisquare(h[start:end], gauss(xs, parameters[0][0], parameters[0][1], parameters[0][2], parameters[0][3]))[1])

# electron_energies = np.array([60, 122, 341, 478, 1061])
# light_yields = np.array([159.17, 341.56, 1072.94, 1479.63, 3328.27])
# light_yields_err = np.array([13.60, 22.75, 47.08, 53.347, 91.525])
#
# # popt = curve_fit(line, electron_energies, light_yields, sigma=light_yields_err, p0=[3, 0])
# popt = curve_fit(lin, electron_energies, light_yields, sigma=light_yields_err, p0=3)
# print(popt)
#
# a = popt[0][0]
# da = np.sqrt(popt[1][0][0])
# # b = popt[0][1]
# # db = np.sqrt(popt[1][1][1])
#
# # print("line equation = ", a, " (+- ", da, ") x", " + ", b, " (+- ", db, ")")
# print("line equation = ", a, " (+- ", da, ") x")
#
# print("chisquare = ",
#       sum(((light_yields[i] - lin(electron_energies[i], a)) ** 2) / lin(electron_energies[i], a) for i in
#           range(len(electron_energies))) / len(electron_energies))
#
# plt.errorbar(electron_energies, light_yields, yerr=light_yields_err, fmt='.', capthick=3)
#
# plt.plot(electron_energies, lin(electron_energies, a), '--')
# plt.grid()
# plt.xlabel("Electron Energies [keV]")
# plt.ylabel("Light Yield [a.u.]")
# plt.title("Light yield as a function of electron recoil energy")
# plt.show()
#
# popt = np.array(curve_fit(lin, electron_energies, light_yields, sigma=light_yields_err))
# pcov = popt[1]
# a = popt[0][0]
# # da = np.sqrt(popt[1][0][0])
# perr = np.sqrt(pcov[0])
#
# # prepare confidence level curves
# nstd = 3.  # to draw 3-sigma intervals
# popt_up = popt + nstd * perr
# a_up = popt_up[0][0]
# popt_dw = popt - nstd * perr
# a_down = popt_dw[0][0]
#
# fit = lin(electron_energies, a)
# fit_up = lin(electron_energies, a_up)
# fit_dw = lin(electron_energies, a_down)
#
# # plot
# fig, ax = plt.subplots(1)
# plt.rcParams['xtick.labelsize'] = 18
# plt.rcParams['ytick.labelsize'] = 18
# plt.rcParams['font.size'] = 20
# plt.errorbar(electron_energies, light_yields, yerr=light_yields_err, fmt='.', capthick=3)
#
# plt.xlabel("Electron Energies [keV]")
# plt.ylabel("Light Yield [a.u.]")
# plt.title("Light yield as a function of electron recoil energy")
# plt.plot(electron_energies, fit, '--', lw=2, label='best fit curve')
# # plt.plot(electron_energies, lin(electron_energies, a), 'k', lw=2, label='True curve')
# ax.fill_between(electron_energies, fit_up, fit_dw, alpha=.25, label='3 - sigma interval')
# plt.legend(loc='lower right', fontsize=18)
# plt.xticks(np.arange(0, 1200, 100))
# plt.grid()
# plt.show()
#
#
# # L_eff(light_yield_n, err_n, line(16.4, a, b), da * 16.4 + db, 16.4)
# # L_eff(light_yield_n, err_n, lin(16.4, a), da * 16.4, 16.4)
#
# # calibrated vs 60keV Photo-peak
# # L_eff(11.3235, 3.2618, 159.1653, np.sqrt(13.59777 ** 2 + (da * 60) ** 2), 11.5, 2.1)
# # L_eff(13.8651, 3.77856, 159.1653, np.sqrt(13.59777 ** 2 + (da * 60) ** 2), 16.4, 2.3)
# # L_eff(18.309, 4.769, 159.1653, np.sqrt(13.59777 ** 2 + (da * 60) ** 2), 28.5, 3.1)
# # L_eff(23.108, 5.649, 159.1653, np.sqrt(13.59777 ** 2 + (da * 60) ** 2), 43.4, 3.7)
# # L_eff(28.097, 6.5, 159.1653, np.sqrt(13.59777 ** 2 + (da * 60) ** 2), 60.5, 4.0)
# # L_eff(41.344, 7.868, 159.1653, np.sqrt(13.59777 ** 2 + (da * 60) ** 2), 119.5, 4.9)
#
# def leff(Nph_n, E_nr):
#     print("L_eff at E = ", E_nr, " keV is equal to ", (Nph_n / 159.16) * 60 / E_nr)
#
#
# # using pencil beam
# # leff(11.3, 11.5)
# # leff(13.9, 16.4)
# # leff(18.3, 28.5)
# # leff(23.1, 43.3)
# # leff(28.1, 60.5)
# # leff(41.3, 119.5)
#
#
# # calibrated vs equivalent er energy, using pencil beam
# # L_eff(11.3235, 3.2618, 3.067 * 11.5, da * 11.5, 11.5, 2.1)
# # L_eff(13.8651, 3.77856, 3.067 * 16.4, da * 16.4, 16.4, 2.3)
# # L_eff(18.309, 4.769, 3.067 * 28.5, da * 28.5, 28.5, 3.1)
# # L_eff(23.108, 5.649, 3.067 * 43.4, da * 43.3, 43.4, 3.7)
# # L_eff(28.097, 6.5, 3.067 * 60.5, da * 60.5, 60.5, 4.0)
# # L_eff(41.344, 7.868, 3.067 * 119.5, da * 119.5, 119.5, 4.9)
#
# # light yield using beam with angular separation, calibrated vs 60keV Photo-peak
# print("using beam with angular separation:\n")
# L_eff(10.18, 4.1, 159.1653, np.sqrt(13.59777 ** 2 + (da * 60) ** 2), 11.5, 2.1)
# L_eff(13.14, 6.18, 159.1653, np.sqrt(13.59777 ** 2 + (da * 60) ** 2), 16.4, 2.3)
# L_eff(17.67, 10, 159.1653, np.sqrt(13.59777 ** 2 + (da * 60) ** 2), 28.5, 3.1)
# L_eff(22.51, 16.13, 159.1653, np.sqrt(13.59777 ** 2 + (da * 60) ** 2), 43.4, 3.7)
# L_eff(26.22, 12.24, 159.1653, np.sqrt(13.59777 ** 2 + (da * 60) ** 2), 60.5, 4.0)
# L_eff(39.34, 26.53, 159.1653, np.sqrt(13.59777 ** 2 + (da * 60) ** 2), 119.5, 4.9)

# L_effs = np.linspace(0.2, 0.4, 100)
# Rs = np.linspace(1, 7, 100)
# L_eff_chi_square("../Scintillation_Data/nucleus_recoil_scint/lightYield_60_wide.csv",
#                  "../Scintillation_Data/Creus/60deg_data.csv", "../Scintillation_Data/Creus/60deg_err_down.csv",
#                  "../Scintillation_Data/Creus/60deg_err_up.csv", 25, 135, 3.75, L_effs, Rs, 60, 5, vmax=300)

# dat = pd.read_csv("../Scintillation_Data/nucleus_recoil_scint/lightYield_90_wide.csv")
# recoil_energies = dat.nucleusRecoilEnergy
# h, bins = np.histogram(recoil_energies, 100)
# plt.step((bins[1:] + bins[:-1]) / 2, h)
# plt.xlabel('recoil energy [keV]')
# plt.ylabel('counts')
# plt.show()

# bins_centers = (bins[1:] + bins[:-1]) / 2
#
# start = 27
# end = 35
# parameters = curve_fit(gauss, bins_centers[start:end], h[start:end],
#                        p0=[5, 35, 120, 6])
# light_yield_n = parameters[0][2]
# err_n = parameters[0][3]
# print("gaussian H + A * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))\n with params H = ", parameters[0][0], ", A = ",
#       parameters[0][1], ", $\mu$ = ", parameters[0][2], ", $\sigma$ = ", parameters[0][3])
#
# xs = bins_centers[start:end]
# plt.plot(xs, gauss(xs, parameters[0][0], parameters[0][1], parameters[0][2], parameters[0][3]), '--',
#          label='Fit result', )
# plt.step(bins_centers, h, where='mid', label='number of photons absorbed in PMT')
# plt.legend()
# plt.show()
#
# print("chisquare = ",
#       chisquare(h[start:end], gauss(xs, parameters[0][0], parameters[0][1], parameters[0][2], parameters[0][3]))[
#           0] / len(h))
# print("pvalue = ",
#       chisquare(h[start:end], gauss(xs, parameters[0][0], parameters[0][1], parameters[0][2], parameters[0][3]))[1])

# 1 is Aviv
# xs1 = [27.65, 42.78, 60.5, 120]
# dx1 = [4.04, 4.09, 5.68, 5.28]
# ys1 = [0.303, 0.333, 0.341, 0.323]
# dy1_up = [0.0197, 0.0524, 0.052, 0.0685]
# dy1_down = [0.081, 0.031, 0.055, 0.0372]
# dy1 = np.array([dy1_down, dy1_up])
#
# # 2 is Creus
# xs2 = [28.5, 43.4, 60.5, 119.5]
# dx2 = [5.2, 6.0, 7.1, 7.7]
# ys2 = [0.285, 0.294, 0.283, 0.301]
# dy2 = [0.014, 0.023, 0.018, 0.019]
#
# plt.errorbar(xs1, ys1, xerr=dx1, yerr=dy1, color='b', ecolor='b', fmt='-', label='Aviv')
# plt.errorbar(xs2, ys2, xerr=dx2, yerr=dy2, color='grey', ecolor='grey', fmt='-', label='Creus et al.')
# plt.xlabel('nuclear recoil energy[keV]', fontsize=15)
# plt.ylabel('L_eff', fontsize=15)
# plt.title('$L_{eff}$ as a function of nuclear recoil energy', fontsize=15)
# plt.legend()
# plt.show()

# dat = pd.read_csv('../Dune_size_5.3MeV_n_escape_E.csv')
# escape_energy_keV = dat.neutronEscapeEnergy
# nPrime = dat.nPrime
#
# tot_events = len(escape_energy_keV)
# print("number of events = ", tot_events)
# escape_energy_keV = escape_energy_keV[escape_energy_keV > 0]
# print("percentage of events with neutrons escaping = ", 100 * len(escape_energy_keV) / tot_events)
# h, bins = np.histogram(escape_energy_keV, np.arange(-5, 550, 5))
# plt.step((bins[1:] + bins[:-1]) / 2, h, where='mid')
# # plt.hist(escape_energy_keV, 100)
# plt.xlabel('energy[keV]', fontsize=15)
# plt.ylabel('count', fontsize=15)
# plt.yscale('log')
# plt.xticks(np.arange(0, 500, 20))
# plt.title('5.3MeV neutrons escape energy starting from center of DUNE', fontsize=15)
# plt.show()
#
# escape_energy_keV = escape_energy_keV[escape_energy_keV < 70]
# escape_energy_keV = escape_energy_keV[escape_energy_keV > 40]
# print("number of escaped neutrons around 60keV = ", len(escape_energy_keV))

# dat = pd.read_csv('../neutronData/5.3MeV_n_E_transfer.csv')
# eventIDs = dat.eventID
# Ar_recoil_transfer = dat.ArRecoil
# gammas_transfer = dat.gammas
# inelasticRecoilTransfer = dat.inelasticRecoilTransfer
#
# theSum = Ar_recoil_transfer + gammas_transfer + inelasticRecoilTransfer
# plt.plot(eventIDs, theSum, '.')
# plt.show()
#
# h, bins = np.histogram(Ar_recoil_transfer, np.arange(-0.05, 5.3, 0.05))
# plt.step((bins[1:] + bins[:-1]) / 2, h, where='mid', label='elastic interaction - energy transfer to nuclear recoils')
#
# h, bins = np.histogram(gammas_transfer, np.arange(-0.05, 5.3, 0.05))
# plt.step((bins[1:] + bins[:-1]) / 2, h, where='mid', label='sum of energy of de-excitation $\gamma$')
#
# h, bins = np.histogram(inelasticRecoilTransfer, np.arange(-0.05, 5.3, 0.05))
# plt.step((bins[1:] + bins[:-1]) / 2, h, where='mid', label='energy transfer to nuclear recoils in inelastic process')
# plt.xlabel('energy transfer[MeV]',fontsize=15)
# plt.ylabel('count',fontsize=15)
# # plt.yscale('log')
# plt.title('histogram of energy transfer contributions',fontsize=15)
# plt.legend()
# plt.show()

# dat = pd.read_csv("../Scintillation_Data/scint_yield_checks/scintDist_976keV_electron_0.414.csv")
# num = dat.numPhotons
# h, bins = np.histogram(num, np.linspace(0, 51300, 500), density=True)
# bins_centers = (bins[1:] + bins[:-1]) / 2
# plt.step(bins_centers, h, label='976 keV electron, 0.414 $[\\frac{mm}{MeV}]$ Birks constant')
# ymax = 1.1 * np.max(h)
# plt.plot([4.41 * 10 ** 4, 4.41 * 10 ** 4], [0, ymax], '--b')
#
# dat = pd.read_csv("../Scintillation_Data/scint_yield_checks/scintDist_976keV_electron_0.696.csv")
# num = dat.numPhotons
# h, bins = np.histogram(num, np.linspace(0, 51300, 500), density=True)
# bins_centers = (bins[1:] + bins[:-1]) / 2
# plt.step(bins_centers, h, label='976 keV electron, 0.696 $[\\frac{mm}{MeV}]$ Birks constant')
# ymax = 1.1 * np.max(h)
# plt.plot([4.175 * 10 ** 4, 4.175 * 10 ** 4], [0, ymax], '--', color='orange')
#
# dat = pd.read_csv("../Scintillation_Data/scint_yield_checks/scintDist_976keV_electron_0.158.csv")
# num = dat.numPhotons
# h, bins = np.histogram(num, np.linspace(0, 51300, 500), density=True)
# bins_centers = (bins[1:] + bins[:-1]) / 2
# plt.step(bins_centers, h, label='976 keV electron, 0.158 $[\\frac{mm}{MeV}]$ Birks constant')
# ymax = 1.1 * np.max(h)
# plt.plot([4.729 * 10 ** 4, 4.729 * 10 ** 4], [0, ymax], '--g')
#
# dat = pd.read_csv("../Scintillation_Data/scint_yield_checks/scintDist_976keV_electron_0.89.csv")
# num = dat.numPhotons
# h, bins = np.histogram(num, np.linspace(0, 51300, 500), density=True)
# bins_centers = (bins[1:] + bins[:-1]) / 2
# plt.step(bins_centers, h, label='976 keV electron, 0.89 $[\\frac{mm}{MeV}]$ Birks constant')
# ymax = 1.1 * np.max(h)
# xmax = bins[np.where(h == np.max(h))[0][0]]
# plt.plot([xmax,xmax], [0, ymax], '--r')
#
# plt.legend()
# plt.xlabel('number of optical photons')
# plt.ylabel('frequency')
# plt.title('number of optical photons produced by 976 keV e-')
# plt.show()
# #
# dat = pd.read_csv("../Scintillation_Data/Calibration_for_coverage/energy_diff.csv")
# energy = dat.energy * 1000
# plt.hist(energy, 100)
# plt.xlabel("energy [eV]")
# plt.ylabel('count')
# plt.title("remaining energy after photoelectric effect caused by 59.5 keV primary $\gamma$")
# plt.show()
# #
# dat = pd.read_csv("../phot_eKin.csv")
# energies = dat.electron_eKin
# plt.hist(energies)
# plt.show()
#
# dat = pd.read_csv("../Scintillation_Data/scint_yield_checks/scintDist_1MeV_gamma_infinite_detector.csv")
# num_photons = dat.numPhotons
# num_photons = num_photons[num_photons > 0]
# h, bins = np.histogram(num_photons, np.linspace(0, 51300, 500), density=True)
# bins_centers = (bins[1:] + bins[:-1]) / 2
# plt.step(bins_centers, h, label='1 MeV gamma, infinite detector')
#
# dat = pd.read_csv("../Scintillation_Data/scint_yield_checks/scintDist_1MeV_electron.csv")
# num_photons = dat.numPhotons
# num_photons = num_photons[num_photons > 0]
# h, bins = np.histogram(num_photons, np.linspace(0, 51300, 500), density=True)
# bins_centers = (bins[1:] + bins[:-1]) / 2
# plt.step(bins_centers, h, label='1 MeV electron')
#
# plt.legend()
# plt.show()

# dat = pd.read_csv("../Scintillation_Data/scint_yield_checks/scintDist_60keV_gamma.csv")
# num_photons = dat.numPhotons
# num_photons = num_photons[num_photons > 0]
# h1, bins1 = np.histogram(num_photons, 100, density=True)
# bins_centers = (bins1[1:] + bins1[:-1]) / 2
# plt.step(bins_centers, h1, label='60 keV $\gamma$ , infinite detector')
# ymax1 = 1.1 * np.max(h1)
# xmax1 = bins1[np.where(h1 == np.max(h1))[0][0]]
# plt.plot([xmax1, xmax1], [0, ymax1], '--', color='blue')
#
# dat = pd.read_csv("../Scintillation_Data/scint_yield_checks/scintDist_60keV_electron.csv")
# num_photons = dat.numPhotons
# num_photons = num_photons[num_photons > 0]
# h, bins = np.histogram(num_photons, 100, density=True)
# bins_centers = (bins[1:] + bins[:-1]) / 2
# plt.step(bins_centers, h, label='60 keV electron , infinite detector')
# ymax = 1.1 * np.max(h)
# xmax = bins[np.where(h == np.max(h))[0][0]]
# plt.plot([xmax, xmax], [0, ymax], '--', color='orange')
# xticks = list(np.linspace(0, 2000, 11))
# xticks.extend([xmax, xmax1])
# xticks.remove(2000)
# xticks.append(2150)
# plt.xticks(xticks)
# plt.ylim([0, 0.0031])
# plt.legend()
# plt.show()

# dat = pd.read_csv("../Scintillation_Data/Calibration_for_coverage/experimentalEnergies_close_detector.csv")
# time_of_flight = dat.timeOfFlight
# success = dat.success
# goal = time_of_flight[success == 1]
# background = time_of_flight[success == 0]
# h, bins = np.histogram(goal, 20)
# bins_center = (bins[1:] + bins[:-1]) / 2
# plt.step(bins_center, h, label='goal')
#
# h, bins = np.histogram(background, 20)
# bins_center = (bins[1:] + bins[:-1]) / 2
# plt.step(bins_center, h, label='background')
#
# plt.xlabel('time [ns]')
# plt.ylabel('frequency')
# plt.legend()
# plt.show()

# dat = pd.read_csv("../neutronData/VitaliyStats.csv")
# event_ids = dat.eventID
# recoil_energies = dat.recoilEnergy
#
# print(min(recoil_energies))
# plt.hist(recoil_energies, 1000)
# plt.xlabel("recoil energy[keV]")
# plt.ylabel("count")
# plt.yscale('log')
# plt.xlim([0, 250])
# plt.title("Ar nucleus recoil energy distribution across 100k events")
# plt.show()
#
# count = 0
# last_id = 0
# num_elastic_scatters = []
# for i in range(len(event_ids)):
#     id = event_ids[i]
#     if id == last_id:
#         count += 1
#     else:
#         num_elastic_scatters.append(count)
#         count = 1
#         last_id = id
# plt.hist(num_elastic_scatters, np.linspace(0, 300, 100))
# plt.xlim([0, 100])
# plt.xlabel("# elastic scatters")
# plt.ylabel("frequency")
# plt.title("number of elastic scatters per events")
# plt.show()
# print(max(num_elastic_scatters))
#
# dat = pd.read_csv("../Scintillation_Data/Calibration_for_coverage/numPhot_option1_gamma_59keV.csv")
# photons = dat.numPhot
# h, bins = np.histogram(photons, 100, density=True)
# bins_centers = (bins[1:] + bins[:-1]) / 2
# plt.step(bins_centers, h, label='$\gamma$')
#
# dat = pd.read_csv("../Scintillation_Data/Calibration_for_coverage/numPhot_59keVgamma_LAr_option4.csv")
# photons = dat.numPhotons
# # photons = photons[photons > 1000]
# h, bins = np.histogram(photons, 100, density=True)
# bins_centers = (bins[1:] + bins[:-1]) / 2
# plt.step(bins_centers, h, label='option4 $\gamma$')
#
# dat = pd.read_csv("../Scintillation_Data/Calibration_for_coverage/numPhot_option1_electron_59keV.csv")
# photons = dat.numPhotons
# # photons = photons[photons > 1000]
# h, bins = np.histogram(photons, 100, density=True)
# bins_centers = (bins[1:] + bins[:-1]) / 2
# # plt.step(bins_centers, h, label='option1 e-')
#
# # dat = pd.read_csv("../Hmmge.csv")
# # photons = dat.numPhotons
# # photons = photons[photons > 1600]
# # h, bins = np.histogram(photons, 100, density=True)
# # bins_centers = (bins[1:] + bins[:-1]) / 2
# # plt.step(bins_centers, h, label='Hmmge')
#
# dat = pd.read_csv("../Scintillation_Data/Calibration_for_coverage/numPhot_option1_electron_59keV.csv")
# photons2 = dat.numPhotons
# photons2 = photons2[photons2 > 1600]
# h2, bins2 = np.histogram(photons2, 100, density=True)
# bins_centers2 = (bins2[1:] + bins2[:-1]) / 2
# # plt.step(bins_centers2, h2, label='e-')
#
# plt.xlabel('number of optical photons produced')
# plt.ylabel('frequency')
# plt.title('number of optical photons produced by 59.5 keV $\gamma$ vs 59.5 keV electron with em_option1 package')
# plt.legend()
# plt.show()
#
# dat = pd.read_csv("../energy_diff.csv")
# energy = dat.energy * 1000
# plt.hist(energy, 100)
# plt.xlabel("energy [eV]")
# plt.ylabel('count')
# plt.title("remaining energy after photoelectric effect caused by 59.5 keV primary $\gamma$")
# plt.show()
#
# dat = pd.read_csv("../Scintillation_Data/Calibration_for_coverage/missing_e_single_pe_29.5_init.csv")
# energies = dat.energy * 1000
# plt.hist(energies, 100, label='29.5 keV initial $\gamma$')
# dat = pd.read_csv("../Scintillation_Data/Calibration_for_coverage/missing_e_single_pe_59.5_init.csv")
# energies = dat.energy * 1000
# plt.hist(energies, 100, label='59.5 keV initial $\gamma$')
# plt.xlabel('energy [eV]')
# plt.ylabel('count')
# plt.title('missing energy following p.e effect with 29.5,59.5 keV initial $\gamma$ , single photoelectron')
# plt.legend()
# plt.show()

# dat = pd.read_csv("../phot_eKin.csv")
# energies = dat.electron_eKin
# high_energy_num = list(energies).count(56.3224)
# lower_energy_num = list(energies).count(59.1866)
# even_lower_energy_num = list(energies).count(59.2529)
# lowerer_energy_num = list(energies).count(59.2506)
# print("number of events with 3.1776 keV difference = ", high_energy_num)
# print("number of events with 313.4 eV difference = ", lower_energy_num)
# print("number of events with 249.4 eV difference =  ", lowerer_energy_num)
# print("number of events with 247.1 eV difference =  ", even_lower_energy_num)
#
# print("ratio = ", lower_energy_num / high_energy_num)
#
# plt.hist(energies, 200)
# plt.xticks(np.linspace(56.2, 59.5, 12))
# plt.xlabel('k.e [keV]')
# plt.ylabel('count')
# plt.yscale('log')
# plt.title('energetic electron kinetic energy following p.e effect with 59.5 keV initial $\gamma$')
# plt.show()
#
# dat = pd.read_csv("../numPhotSec_3keV.csv")
# num = dat.numPhot
# print(max(num))
# h, bins = np.histogram(num, np.linspace(0, 30, 15))
# bins_centers = (bins[1:] + bins[:-1]) / 2
# plt.step(bins_centers, h, label='3 keV secondary electron')
# plt.xlabel('number of photons')
# plt.ylabel('frequency')
# plt.title('number of optical photons produced by secondary 3 keV electrons')
# plt.show()
#
# dat = pd.read_csv("../Scintillation_Data/scint_yield_checks/scintDist_56keV_electron_0.89.csv")
# num = dat.numPhotons
# h, bins = np.histogram(num, np.linspace(0, 2000, 50), density=True)
# bins_centers = (bins[1:] + bins[:-1]) / 2
# plt.step(bins_centers, h, label='56 keV primary electron')
#
# dat = pd.read_csv("../Scintillation_Data/scint_yield_checks/numPhotSec_56keV.csv")
# num = dat.numPhot[dat.numPhot > 0]
# # print(len(num[num == 0]))
# # print(len(num))
# h, bins = np.histogram(num, np.linspace(0, 2000, 50), density=True)
# bins_centers = (bins[1:] + bins[:-1]) / 2
# plt.step(bins_centers, h, label='56 keV secondary electron')
#
# dat = pd.read_csv("../Scintillation_Data/scint_yield_checks/scintDist_56keV_gamma_0.89.csv")
# num = dat.numPhotons
# h, bins = np.histogram(num, np.linspace(0, 2000, 50), density=True)
# bins_centers = (bins[1:] + bins[:-1]) / 2
# plt.step(bins_centers, h, label='56 keV $\gamma$')
#
# plt.legend()
# plt.xlabel('number of photons')
# plt.ylabel('frequency')
# plt.title('number of optical photons produced by primary vs secondary 56 keV electrons vs 56 keV $\gamma$')
# plt.show()
#
# dat = pd.read_csv("../Scintillation_Data/Calibration_for_coverage/scintDist_gamma.csv")
# num_photons = dat.numPhotons
# num_photons = num_photons[num_photons > 0]
# h1, bins1 = np.histogram(num_photons, 100, density=True)
# bins_centers = (bins1[1:] + bins1[:-1]) / 2
# plt.step(bins_centers, h1, label='60 keV $\gamma$ , infinite detector')
# ymax1 = 1.1 * np.max(h1)
# xmax1 = bins1[np.where(h1 == np.max(h1))[0][0]]
# plt.plot([xmax1, xmax1], [0, ymax1], '--', color='blue')
#
# dat = pd.read_csv("../Scintillation_Data/Calibration_for_coverage/scintDist_e.csv")
# num_photons = dat.numPhotons
# num_photons = num_photons[num_photons > 0]
# h, bins = np.histogram(num_photons, 100, density=True)
# bins_centers = (bins[1:] + bins[:-1]) / 2
# plt.step(bins_centers, h, label='60 keV electron , infinite detector')
# ymax = 1.1 * np.max(h)
# xmax = bins[np.where(h == np.max(h))[0][0]]
# plt.plot([xmax, xmax], [0, ymax], '--', color='orange')
# xticks = list(np.linspace(0, 2000, 11))
# xticks.extend([xmax, xmax1])
# xticks.remove(2000)
# xticks.append(2150)
# plt.xticks(xticks)
# plt.ylim([0, 0.0031])
# plt.legend()
# plt.show()

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


# dat = pd.read_csv("../singleElasticAngles.csv")
# angles = dat.angle
# plt.hist(angles, 200, density=True)
# plt.xlabel('angle [rad]')
# plt.ylabel('frequency')
# plt.show()
#
# dat2 = pd.read_csv("../Scintillation_Data/Calibration_for_coverage/elasticPhot.csv")
# num_photons = dat2.numPhotonsCreated
#
# h1, bins1 = np.histogram(num_photons, np.linspace(200, 600, 60))
# mean_N_ph = np.mean(num_photons)
# std_N_ph = np.std(num_photons)
# plt.plot([mean_N_ph, mean_N_ph], [0, np.max(h1)], '--r')
# plt.text(mean_N_ph + 2 * std_N_ph, np.max(h1) * 0.8, "mean =%.1f, std = %.1f" % (mean_N_ph, std_N_ph), fontsize=10,
#          color='red')
# bins1_centers = (bins1[1:] + bins1[:-1]) / 2

# parameters = curve_fit(gauss, bins1_centers, h1,
#                        p0=[0, 350, 420, 15])
# miu = parameters[0][2]
# d_miu = parameters[1][2][2]
# sigma = parameters[0][3]
# d_sigma = parameters[1][3][3]
#
# bins1_centers = bins1_centers
# h1 = h1
# plt.step(bins1_centers, h1, where='mid')
# plt.plot(bins1_centers, gauss(bins1_centers, parameters[0][0], parameters[0][1], miu, sigma),
#          '--',
#          label='Fit result', )
# plt.xlabel('num photons created')
# plt.ylabel('count')
# plt.title('num photons created by 45-55 keV Ar recoil nucleus')
# plt.show()

#
# print("gaussian H + A * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))\n with params H = ", parameters[0][0], ", A = ",
#       parameters[0][1], ", $\mu$ = ", miu, ", $\sigma$ = ", sigma)

# dat = pd.read_csv("../dEdx_43keV.csv")
# # in keV
# init_energy = 1000 * dat.initialEnergy
# dE_dx = dat.dEdx
#
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
#
# energies = energies[avg_dE_dx < 10]
# avg_dE_dx = avg_dE_dx[avg_dE_dx < 10]
# plt.plot(energies, avg_dE_dx, '.')
# plt.xlabel('initial energy[MeV]')
# plt.ylabel('dE/dx[MeV/cm]')
# plt.yscale('log')
# plt.show()
#
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

# dat = pd.read_csv("../photLoc.csv")
# dat2 = pd.read_csv("../eventTypes2.csv")
#
# types = dat2.type
# xPos = dat.xPos
# yPos = dat.yPos
# zPos = dat.zPos
# event_ids = dat.eventID
# #
# num_pes_single_elastic = []
# num_pes_multiple_elastic = []
# num_pes_elastics_nnprime = []
# num_pes_multiple_nnprime = []
# num_pes_nnPrime = []
#
# last_eventID = -1
# same = 0
#
# event_dict = {"single_elastic": num_pes_single_elastic, "multiple_elastic": num_pes_multiple_elastic,
#               "elastics_plus_nnprime": num_pes_elastics_nnprime, "multiple_nnprime": num_pes_multiple_nnprime,
#               "nnPrime": num_pes_nnPrime}
#
# for i in range(len(event_ids)):
#     if i % 1000 == 0:
#         print(event_ids[i])
#
#     event_id = event_ids[i]
#     event_type = types[event_id]
#
#     if event_id == last_eventID:
#         same = 1
#     else:
#         same = 0
#         last_eventID = event_id
#     for event in ["single_elastic", "multiple_elastic", "elastics_plus_nnprime", "multiple_nnprime", "nnPrime"]:
#         if event_type == event:
#             pes_list = event_dict[event]
#             if same:
#                 if is_in_arapucas(yPos[i], zPos[i]):
#                     pes_list[len(num_pes_single_elastic) - 1] += 1
#             else:
#                 if is_in_arapucas(yPos[i], zPos[i]):
#                     pes_list.append(1)
#                 else:
#                     pes_list.append(0)
#
# np.save("single_elastic_pes.npy", num_pes_single_elastic)
# np.save("multiple_elastic_pes.npy", num_pes_multiple_elastic)
# np.save("elastics_nnprime_pes.npy", num_pes_elastics_nnprime)
# np.save("multiple_nnprime_pes.npy", num_pes_multiple_nnprime)
# np.save("nnPrime_pes", num_pes_nnPrime)

# num_pes_single_elastic = np.load("../ColdBox_Creus_pe_data/ColdBox/number_of_pe/single_elastic_pes.npy")
# num_pes_multiple_elastic = np.load("../ColdBox_Creus_pe_data/ColdBox/number_of_pe/multiple_elastic_pes.npy")
# num_pes_elastics_nnprime = np.load("../ColdBox_Creus_pe_data/ColdBox/number_of_pe/elastics_nnprime_pes.npy")
# num_pes_multiple_nnprime = np.load("../ColdBox_Creus_pe_data/ColdBox/number_of_pe/multiple_nnprime_pres.npy")

# dat = pd.read_csv("../Scintillation_Data/ColdBox Simulation/photLoc.csv")
# tAbsorbed = dat.tAbsorbed
# eventIDs = dat.eventID
#
# # print(num_pes_elastics_nnprime)
# ids1 = num_pes_elastics_nnprime > 1
# ids2 = num_pes_single_elastic > 1
# ids3 = num_pes_multiple_elastic > 1
# times1 = []
# times2 = []
# times3 = []
# last_found_idx = 0
# for j in range(100000):
#
#     if j % 100 == 0:
#         print(j)
#
#     if ids1[j] or ids2[j] or ids3[j]:
#         found = False
#         temp = []
#         for i in range(last_found_idx, len(tAbsorbed)):
#             if eventIDs[i] == j:
#                 found = True
#                 temp.append(tAbsorbed[i])
#             if found and eventIDs[i] != j:
#                 last_found_idx = i
#                 break
#         temp = np.array(temp)
#         temp = temp - np.min(temp)
#         temp = list(temp)
#         # print(temp)
#         if ids1[j]:
#             times1.extend(temp)
#         elif ids2[j]:
#             times2.extend(temp)
#         elif ids3[j]:
#             times3.extend(temp)

# print(times)
# np.save("Coldbox_elastics_nnprime_photon_timings", times1)
# np.save("Coldbox_single_elastic_photon_timings", times2)
# np.save("Coldbox_multiple_elastic_photon_timings", times3)
# times1.extend(times2)
# times1.extend(times3)
# np.save("ColdBox_photon_timings", times1)
# times1 = list(np.load("../ColdBox_Creus_pe_data/ColdBox/Timing/Coldbox_elastics_nnprime_photon_timings.npy"))
# times2 = list(np.load("../ColdBox_Creus_pe_data/ColdBox/Timing/Coldbox_single_elastic_photon_timings.npy"))
# times3 = list(np.load("../ColdBox_Creus_pe_data/ColdBox/Timing/Coldbox_multiple_elastic_photon_timings.npy"))
# # times1.extend(times2)
# # times1.extend(times3)
# bins = list(np.linspace(0, 6000, 600))
# # bins1 = list(np.linspace(200, 6000, 200))
# # bins.extend(bins1)
# # bins = np.array(bins)
# # ticks1 = list(np.linspace(0, 1000, 6))
# # ticks2 = list(np.linspace(1000, 6000, 15))
# # ticks1.extend(ticks2)
# # plt.xticks(ticks1)
# h1, bins1 = np.histogram(times1, bins, density=False)
# h2, bins2 = np.histogram(times2, bins, density=False)
# h3, bins3 = np.histogram(times3, bins, density=False)
# #
# plt.step((bins1[1:] + bins1[:-1]) / 2, h1, label="elastics+ (n,n')", color='red')
# plt.step((bins2[1:] + bins2[:-1]) / 2, h2, label="single elastic", color='blue')
# plt.step((bins3[1:] + bins3[:-1]) / 2, h3, label="multiple elastic", color='purple')
# # plt.hist(times1, bins)
#
# # print("multiple elastic prompt fraction in ColdBox = ", sum(h[0][0:300]) / sum(h[0]))
#
# elastics_nnprime_pf = sum(h1[0:10]) / sum(h1)
# single_elastic_pf = sum(h2[0:10]) / sum(h2)
# multiple_elastic_pf = sum(h3[0:10]) / sum(h3)
#
# print("single elastic prompt fraction = ", single_elastic_pf)
# print("multiple elastic prompt fraction = ", multiple_elastic_pf)
# print("elastics + nnprime prompt fraction = ", elastics_nnprime_pf)
#
# plt.xlabel('detection time[ns]')
# plt.xlim([0, 6000])
# plt.ylabel('number of photons')
# plt.title('photon detection time in ColdBox detector')
# plt.yscale('log')
# plt.legend()
# plt.show()

# dat = pd.read_csv("../ColdBox_Creus_pe_data/Creus/DecayTime.csv")
# dat2 = pd.read_csv("../ColdBox_Creus_pe_data/Creus/eventType.csv")
# times = dat.time
# event_ids = dat.eventID
# event_types = dat2.Type
# event_ids_types = dat2.eventID
# my_dict = {}
# for i in range(len(event_ids_types)):
#     if event_ids_types[i] % 1000 == 0:
#         print(event_ids_types[i])
#     my_dict[event_ids_types[i]] = event_types[i]
# times_single_elastic = []
# times_multiple_elastic = []
# times_nnprime = []
# times_elastics_nnprime = []
# event_ids_types_nums = set(event_ids_types)
# for i in range(len(times)):
#     if event_ids[i] % 1000 == 0:
#         print(event_ids[i])
#     if event_ids[i] in event_ids_types_nums:
#         if my_dict[event_ids[i]] == "single_elastic":
#             times_single_elastic.append(times[i])
#         elif my_dict[event_ids[i]] == "multiple_elastic":
#             times_multiple_elastic.append(times[i])
#         elif my_dict[event_ids[i]] == "nnprime":
#             times_nnprime.append(times[i])
#         elif my_dict[event_ids[i]] == "elastics_plus_nnprime":
#             times_elastics_nnprime.append(times[i])
# print(np.min(times_single_elastic))
#
# times_single_elastic = times_single_elastic - np.min(times_single_elastic)
# times_multiple_elastic = times_multiple_elastic - np.min(times_multiple_elastic)
# times_nnprime = times_nnprime - np.min(times_nnprime)
# times_elastics_nnprime = times_elastics_nnprime - np.min(times_elastics_nnprime)

# np.save("single_elastic_times.npy", times_single_elastic)
# np.save("multiple_elastic_times.npy", times_multiple_elastic)
# np.save("elastics_nnprime_times.npy", times_elastics_nnprime)
# np.save("nnprime_times.npy", times_nnprime)

# times_single_elastic = np.load("../ColdBox_Creus_pe_data/Creus/Timing/single_elastic_times.npy")
# times_multiple_elastic = np.load("../ColdBox_Creus_pe_data/Creus/Timing/multiple_elastic_times.npy")
# times_elastics_nnprime = np.load("../ColdBox_Creus_pe_data/Creus/Timing/elastics_nnprime_times.npy")
# times_nnprime = np.load("../ColdBox_Creus_pe_data/Creus/Timing/nnprime_times.npy")

# num_pes_single_elastic = np.load("single_elastic_pes.npy")
# num_pes_multiple_elastic = np.load("multiple_elastic_pes.npy")
# num_pes_elastics_nnprime = np.load("elastics_nnprime_pes.npy")
# num_pes_multiple_nnprime = np.load("multiple_nnprime_pes.npy")
# num_pes_nnPrime = np.load("nnPrime_pes.npy")

# bins = list(np.linspace(0, 6000, 600))
# # bins1 = list(np.linspace(200, 6000, 100))
# # bins.extend(bins1)
# bins = np.array(bins)
# #
# h1, bins1 = np.histogram(times_single_elastic, bins, density=False)
# h2, bins2 = np.histogram(times_multiple_elastic, bins, density=False)
# h3, bins3 = np.histogram(times_nnprime, bins, density=False)
# h4, bins4 = np.histogram(times_elastics_nnprime, bins, density=False)
#
# # for h in [h1, h2, h3, h4]:
# #     print(h[199])
#
# single_elastic_pf = sum(h1[0:10]) / sum(h1)
# multiple_elastic_pf = sum(h2[0:10]) / sum(h2)
# nnprime_pf = sum(h3[0:10]) / sum(h3)
# elastics_nnprime_pf = sum(h4[0:10]) / sum(h4)
#
# print("single elastic prompt fraction = ", single_elastic_pf)
# print("multiple elastic prompt fraction = ", multiple_elastic_pf)
# print("nnprime prompt fraction = ", nnprime_pf)
# print("elastics nnprime prompt fraction = ", elastics_nnprime_pf)
#
# plt.step((bins1[1:] + bins1[:-1]) / 2, h1, label="single elastic", color='red')
# plt.step((bins2[1:] + bins2[:-1]) / 2, h2, label="multiple elastic", color='blue')
# plt.step((bins3[1:] + bins3[:-1]) / 2, h3, label="(n,n')", color='purple')
# plt.step((bins4[1:] + bins4[:-1]) / 2, h4, label="elastics + (n,n')", color='orange')
# plt.xlabel('detection time[ns]')
# # plt.xlim([0, 200])
# # plt.xticks(ticks1)
# plt.ylabel('number of photons')
# plt.title('photon detection time in Creus. et al detector')
# plt.legend()
# plt.yscale('log')
# plt.show()

# for i in range(len(num_pes_multiple_nnprime)):
#     print(num_pes_multiple_nnprime[i])

# plt.hist(num_pes_single_elastic, np.linspace(0, 100, 100), density=True, align='left')
# plt.xlabel('# p.e')
# plt.ylabel('frequency')
# plt.title("number of photoelectrons detected in ColdBox in single elastic events")
# plt.yscale('log')
# plt.xticks(np.linspace(0, 100, 21))
# plt.show()
#
# plt.hist(num_pes_multiple_elastic, np.linspace(0, 100, 100), density=True, align='left')
# plt.xlabel('# p.e')
# plt.ylabel('frequency')
# plt.title("number of photoelectrons detected in ColdBox in multiple elastic events")
# plt.yscale('log')
# plt.xticks(np.linspace(0, 100, 21))
# plt.show()
#
# plt.hist(num_pes_elastics_nnprime, np.linspace(0, 100, 100), density=True)
# plt.xlabel('# p.e')
# plt.ylabel('frequency')
# plt.title("number of photoelectrons detected in ColdBox in elastics + (n,n') events")
# plt.yscale('log')
# plt.xticks(np.linspace(0, 100, 21))
# plt.show()

# zeros = num_pes_single_elastic[num_pes_single_elastic == 0]
# print(len(zeros))

# h1, bins1 = np.histogram(num_pes_single_elastic, np.linspace(0, np.max(num_pes_single_elastic), 1000), density=True)
# h2, bins2 = np.histogram(num_pes_multiple_elastic, np.linspace(0, np.max(num_pes_multiple_elastic), 1000),
#                          density=True)
# h3, bins3 = np.histogram(num_pes_elastics_nnprime, np.linspace(0, np.max(num_pes_elastics_nnprime), 1000),
#                          density=True)
# h4, bins4 = np.histogram(num_pes_nnPrime, np.linspace(0, np.max(num_pes_nnPrime), 1000), density=True)
# # plt.step((bins1[1:] + bins1[:-1]) / 2, h1, label="single elastic", color='red')
# # plt.step((bins2[1:] + bins2[:-1]) / 2, h2, label="multiple elastic", color='blue')
# # plt.step((bins3[1:] + bins3[:-1]) / 2, h3, label="elastics + (n,n')", color='purple')
# plt.step((bins4[1:] + bins4[:-1]) / 2, h4, label="nnPrime", color='pink')
# plt.xlabel('# p.e')
# plt.ylabel('frequency')
# plt.legend()
# plt.title("number of photoelectrons detected in ColdBox by event type")
# plt.yscale('log')
# plt.xlim([0, 250])
# # plt.xscale('log')
# # plt.xticks(np.linspace(0, 100, 21))
# plt.show()

# plt.hist(num_pes_multiple_nnprime, np.linspace(0, 100, 100), density=True)
# plt.xlabel('# p.e')
# plt.ylabel('frequency')
# plt.title("number of photoelectrons detected in ColdBox in multiple (n,n') events")
# plt.yscale('log')
# plt.xticks(np.linspace(0, 100, 21))
# plt.show()

# arr2 = [4,5]
# print(type(arr2) == list)
# arr1 = np.array([0, 3, 5, 6, 0, 3])
# arr2 = np.array([1, 2, 3, 4, 5, 6])
# print(arr2[arr1 != 0])