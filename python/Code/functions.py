from collections import Counter

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.stats import poisson
from scipy.stats import chisquare
from scipy.integrate import simpson

NUM_PARTICLES = 10 ** 4
AVOGADRO_NUMBER = 6.022 * 10 ** 23


def cross_section(N_events, N_incident, thickness, rho=1.4, A=40, string=""):
    print(string, "cross section = ", 10 ** 24 * (N_events * A) / (N_incident * rho * thickness * AVOGADRO_NUMBER))
    print(string, "std = ", sig_cross_section(N_events, N_incident, thickness, rho, A))


def sig_cross_section(N_events, N_incident, thickness, rho=1.4, A=40):
    return ((10 ** 24 * A) / (AVOGADRO_NUMBER * rho * thickness)) * (N_events / N_incident) * np.sqrt(
        1 / N_events + 1 / N_incident)


def plot_and_print(filename, energy):
    distribution = pd.read_csv(filename)
    distances = distribution.distance
    # because particles start 10cm behind the absorber
    for i in range(len(distances)):
        distances[i] -= 10

    plt.hist(distances, 40)
    title = str(energy) + " MeV photons fp through liquid argon"
    plt.title(title)
    plt.ylabel("# of events")
    plt.xlabel("distance traveled [cm]")
    plt.show()
    print("mean value of mean free path of ", energy, "MeV photons is :",
          np.mean(distances), " cm")
    print("std of ", energy, "MeV photons mean free path distribution is: ", np.std(distances), " cm")


def model_func(x, miu):
    return np.e ** (-x * miu)


def sig(ys):
    return [ys[i] * np.sqrt((1 / (ys[i] * NUM_PARTICLES)) + (1 / NUM_PARTICLES)) for i in range(len(ys))]


def estimate_miu_and_std(model, xs, ys, miu0):
    fit = curve_fit(model, xs, ys, miu0, sigma=sig(ys))
    print("estimated miu = ", fit[0][0])
    print("estimated std of miu = ", np.sqrt(fit[1][0][0]))


def two_d_plot(xDist, yDist, title, xlabel, ylabel, xtext, ytext, font, xlim=None, ylim=None, integer_xDist=False):
    if integer_xDist:

        bins_x = np.arange(min(xDist), max(xDist) + 2) - 0.5
        bins_y = np.linspace(0, max(yDist), len(bins_x))
        hist, xbins, ybins, im = plt.hist2d(xDist, yDist, bins=[bins_x, bins_y],
                                            cmap=plt.cm.jet)
        summ = sum(hist[i, j] for i in range(len(hist)) for j in range(len(hist[0])))
        plt.xticks(np.arange(min(xDist), max(xDist) + 1, 1), fontsize=8)

        for i in range(len(bins_x) - 1):
            for j in range(len(bins_y) - 1):
                plt.text(xbins[i] + xtext, ybins[j] + ytext, s="{:.1f}%".format(100 * hist[i, j] / summ)
                         , color="w", ha="center",
                         va="center", fontweight="bold", fontsize=font)

    else:
        hist, xbins, ybins, im = plt.hist2d(xDist, yDist, bins=[np.linspace(min(xDist), max(xDist), 50),
                                                                np.arange(0, 50, 1)],
                                            cmap=plt.cm.jet)
        summ = sum(hist[i, j] for i in range(len(hist)) for j in range(len(hist[0])))
        plt.xticks(np.arange(min(xDist), max(xDist) + 1, (max(xDist) - min(xDist)) / 50), fontsize=7)
        for i in range(len(xbins) - 1):
            for j in range(len(ybins) - 1):
                plt.text(xbins[i] + xtext, ybins[j] + ytext, s="{:.1f}%".format(100 * hist[i, j] / summ)
                         , color="w", fontweight="bold", fontsize=font)

    print("sum of hist[i,j] = ", summ)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.show()
    return hist, xbins, ybins, im


def plot_photon_dist(photon_dist, events_ids, title):
    num_photons_per_event = [0 for _ in range(10 ** 4)]
    num_photons = 0

    for i in range(photon_dist.size):

        if i == 0:
            num_photons = 1
        elif events_ids[i] == events_ids[i - 1]:
            num_photons += 1
        elif events_ids[i] != events_ids[i - 1]:
            num_photons_per_event[events_ids[i - 1]] = num_photons
            num_photons = 1

    print(min(num_photons_per_event))
    bins_x = np.arange(min(num_photons_per_event), max(num_photons_per_event) + 2) - 0.5
    plt.xticks(np.arange(min(num_photons_per_event), max(num_photons_per_event) + 1, 1))
    plt.xlabel("number of photons created")
    plt.ylabel("number of events")
    plt.yscale("log")
    plt.title(title)
    plt.hist(num_photons_per_event, bins=bins_x, rwidth=0.3)
    plt.show()


def find_levels(filename, energy, epsilon=0.5):
    data = pd.read_csv(filename)
    energies = data.energy_keV

    for i in range(len(energies)):
        for j in range(i + 1, len(energies)):

            if abs((float(energies[j]) - float(energies[i])) - energy * 1000) < epsilon:
                print(energies[i])
                print(energies[j])


def poisson_f(k, lamb):
    '''poisson function, parameter lamb is the fit parameter'''
    return poisson.pmf(k, lamb)


def gauss(x, H, A, x0, sigma):
    return H + A * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))


def L_eff(light_yield_n, err_n, light_yield_el, err_el, energy, err_energy):
    err = np.sqrt(
        (err_n / light_yield_el) ** 2 + (light_yield_n * err_el / (light_yield_el ** 2)) ** 2 + (
                (light_yield_n * err_energy) / (light_yield_el * energy ** 2)) ** 2)
    print("L_eff at E = ", energy, " keV is equal to ", (light_yield_n * 60) / (light_yield_el * energy), " +- ", err)


def plot_pmts(filename, energy, pp, cut, initial_params, valley=0):
    dat = pd.read_csv(filename)
    num_up = np.array(dat.numUp)
    num_down = np.array(dat.numDown)
    tot_num_photons = np.array(dat.numPhotons)
    eDep = np.array(dat.eDep)

    num_up = num_up[num_up > 0]
    num_down = num_down[num_down > 0]
    tot_num_photons = tot_num_photons[tot_num_photons > 0]
    eDep = eDep[eDep > 0]

    h_up, bins_up = np.histogram(num_up, np.linspace(min(num_up), max(num_up), 60))
    plt.step((bins_up[1:] + bins_up[:-1]) / 2, h_up, label='number of photons reached top PMT')

    h_down, bins_down = np.histogram(num_down, np.linspace(min(num_up), max(num_up), 60))
    plt.step((bins_down[1:] + bins_down[:-1]) / 2, h_down, label='number of photons reached bottom PMT')

    plt.xlabel("number of scintillation photons")
    plt.ylabel('number of events')
    plt.title("number of scintillation photons absorbed in PMTs following " + str(energy) + "keV $\\gamma$")
    plt.legend(loc='upper right', fontsize=5)
    plt.show()

    # print("mean number of photons absorbed in top PMT = ", np.mean(h_up), " +- ", np.std(h_up))
    # print("mean number of photons absorbed in bottom PMT = ", np.mean(h_down), " +- ", np.std(h_down))

    bins_up = bins_up[1:]
    bins_down = bins_down[1:]
    h = (h_up + h_down) / 2
    h = h[bins_up > cut]
    h_up = h_up[bins_up > cut]
    h_down = h_down[bins_down > cut]
    bins_up = bins_up[bins_up > cut]
    bins_down = bins_down[bins_down > cut]

    if pp:
        # print("number of scintillation photons absorbed in top PMT at ", energy, " keV photopeak = ",
        #       bins_up[h_up == max(h_up)], " +- ",
        #       bins_up[1] - bins_up[0])
        # print("number of scintillation photons absorbed in bottom PMT ", energy, " keV photopeak = ",
        #       bins_down[h_down == max(h_down)],
        #       " +- ",
        #       bins_down[1] - bins_down[0])
        bins_up = bins_up[h >= max(h) / 7]
        h = h[h >= max(h) / 7]
        bins_up_centers = (bins_up[1:] + bins_up[:-1]) / 2
        parameters = curve_fit(gauss, bins_up_centers, h[1:],
                               p0=initial_params)
        miu = parameters[0][2]
        d_miu = parameters[1][2][2]
        sigma = parameters[0][3]
        d_sigma = parameters[1][3][3]
        print("chisquare = ",
              chisquare(h[1:],
                        gauss(bins_up_centers, parameters[0][0], parameters[0][1], miu, sigma))[
                  0] / len(h[1:]))
        print("pvalue = ", chisquare(h[1:], gauss(bins_up_centers, parameters[0][0], parameters[0][1], miu,
                                                  sigma))[1])
        print("gaussian H + A * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))\n with params H = ", parameters[0][0],
              ", A = ",
              parameters[0][1], ", $\\mu$ = ", miu, ", $\\sigma$ = ", sigma)
        print("miu = ", miu, " +- ", np.sqrt(sigma ** 2 + d_sigma ** 2 + d_miu ** 2))
        plt.plot(bins_up_centers,
                 gauss(bins_up_centers, parameters[0][0], parameters[0][1], parameters[0][2], parameters[0][3]),
                 linestyle='--', label='Fit result', )
        plt.step(bins_up_centers, h[1:], label='Geant4 simulation')
        plt.xlabel('number of photons')
        plt.ylabel('number of events')
        plt.title('number of photons absorbed in PMT histogram')
        plt.legend()
        plt.show()

        # print("number of scintillation photons absorbed in PMTs at ", energy, " keV photopeak = ",
        #       bins_up[h == max(h)], " +- ",
        #       bins_up[1] - bins_up[0])
    else:
        possible_bins = []
        bin_size = bins_up[1] - bins_up[0]
        print("theoretical compton edge height = ", (valley + max(h)) / 2)
        for i in range(len(bins_up)):
            if bins_up[i] > bins_up[h == max(h)] and h[i] <= (valley + max(h)) / 2 + bin_size and \
                    h[i] >= (valley + max(h)) / 2 - bin_size:
                possible_bins.append(bins_up[i])
        for i in range(len(possible_bins)):
            print("height of possible bins = ", h[bins_up == possible_bins[i]])
        print("number of scintillation photons absorbed in PMTs at ", energy, " keV compton edge = ",
              possible_bins,
              " +- ",
              bin_size)

    # h, bins = np.histogram(tot_num_photons, np.linspace(min(tot_num_photons), max(tot_num_photons), 80))
    # plt.step((bins[1:] + bins[:-1]) / 2, h, label='number of photons created')
    #
    h, bins = np.histogram(eDep, np.linspace(min(eDep), max(eDep) + 50, 80))
    plt.step((bins[1:] + bins[:-1]) / 2, h, label='normalized energy deposit in sensitive volume')
    #
    plt.xlabel("energy deposited [keV]")
    # plt.xlabel("number of photons")
    plt.ylabel("number of events")
    plt.title("energy deposited by " + str(energy) + "keV $\\gamma$ [keV]")
    # plt.title("total number of scintillation photons created")
    plt.legend()
    plt.show()


def line(x, a, b):
    return a * x + b


def lin(x, a):
    return a * x


def percentage_detected(filename):
    data = pd.read_csv("../Scintillation_Data/liquid_scintillator/" + filename)
    detecteds = data.eDep[data.eDep == 1]
    return len(detecteds) / len(data.eDep)


def calc_chi_square_sigma(chi_squares, L_effs, Rs, dof, min_chi_square, min_L, min_R, sigmas):
    min_L_idx = list(L_effs).index(min_L)
    min_R_idx = list(Rs).index(min_R)
    right_sigma = 0
    up_sigma = 0
    left_sigma = 0
    down_sigma = 0
    for i in range(min_L_idx, len(chi_squares)):
        for j in range(min_R_idx, len(chi_squares[0])):
            if j == min_R_idx and right_sigma == 0 and chi_squares[i, j] - min_chi_square > sigmas * dof:
                right_sigma = L_effs[i] - min_L
            if i == min_L_idx and up_sigma == 0 and chi_squares[i, j] - min_chi_square > sigmas * dof:
                up_sigma = Rs[j] - min_R

    for i in range(min_L_idx, -1, -1):
        for j in range(min_R_idx, -1, -1):
            if j == min_R_idx and left_sigma == 0 and chi_squares[i, j] - min_chi_square > sigmas * dof:
                left_sigma = min_L - L_effs[i]
            if i == min_L_idx and down_sigma == 0 and chi_squares[i, j] - min_chi_square > sigmas * dof:
                down_sigma = min_R - Rs[j]

    return left_sigma, right_sigma, down_sigma, up_sigma


def L_eff_chi_square(light_yield_filename, data_filename, err_down_filename, err_up_filename, fit_start, fit_end, L_y,
                     L_effs, Rs, scatter_angle, bin_size, vmax=300):
    dat = pd.read_csv(light_yield_filename)

    tot_num_photons = dat.numPhotons
    coincidence_time = dat.coincidenceTime
    num_elastic_sensitive = dat.numElasticSensitive
    num_up = dat.numUp
    num_down = dat.numDown
    scattered_inelastic = dat.scatteredInelastically
    external_scatter = dat.scatteredNotSensitive
    recoil_energies = dat.nucleusRecoilEnergy

    dat = pd.read_csv(data_filename)
    dat2 = pd.read_csv(err_down_filename)
    dat3 = pd.read_csv(err_up_filename)
    num_photons = dat.num_photons
    counts = dat.counts
    errs_down = counts - dat2.counts
    errs_up = dat3.counts - counts
    errs = np.sqrt(np.array(errs_up) ** 2 + np.array(errs_down) ** 2)
    # convolution with a gaussian to account for the gain fluctuation of the PMTs, R_p.e is rms of the gaussian
    R_p_e = 0.4

    min_chi_square = 10 ** 6
    chi_squares = np.array([[0 for _ in range(len(Rs))] for _ in range(len(L_effs))])
    min_bins = []
    min_h = []
    best_L_eff = 0
    best_R = 0
    i = 0
    j = 0
    for R in Rs:
        for L_eff in L_effs:
            # creating new arrays which are smeared according to gaussian dist with sigma around original array
            new_energies = []
            for k in range(len(recoil_energies)):
                new_energy = recoil_energies[k] * L_eff * L_y
                sigma1 = R * np.sqrt(new_energy)
                sigma2 = R_p_e * np.sqrt(new_energy)
                new_sigma = np.sqrt(sigma1 ** 2 + sigma2 ** 2)
                new_energies.append(np.random.normal(new_energy, new_sigma))

                # calculating maximal distance L_eff and R travel over 1 sigma
                # new_energies.append(new_energy)
            new_energies = np.array(new_energies)

            # energies in KeVee
            # new_energies = np.array(recoil_energies) * L_eff
            # gaussian_energies = []
            # for k in range(len(recoil_energies)):
            #     gaussian_energies.append(np.random.normal(new_energies[k], R * np.sqrt(new_energies[k])))
            # new_energies = np.convolve(new_energies, gaussian_energies) * L_y

            h, bins = np.histogram(new_energies, np.arange(fit_start, max(new_energies), bin_size))
            # change h and bins to end on the last bin in the fit
            h = h[:(fit_end // bin_size)]
            bins = bins[:(fit_end // bin_size)]
            # plt.step((bins1[1:] + bins1[:-1]) / 2, h1, where='mid', label='number of photons absorbed in top PMT')

            #
            # plt.xlabel('number of photons')
            # plt.ylabel('number of events')
            # plt.title('number of photons absorbed in bottom and top PMT per event')
            # plt.legend(loc='upper right', fontsize=8)
            # plt.show()

            # match integrals over data points and monte carlo
            I_data = sum(counts[i] * (num_photons[i + 1] - num_photons[i]) for i in range(len(counts) - 1))
            h = h / (bin_size * sum(h)) * I_data
            # print(sum(h))
            # print(sum(counts))
            # print("L_eff = ", L_eff, " , R = ", R)
            chi_square = sum(((counts[i] - h[i]) ** 2) / (errs[i] ** 2 + h[i]) for i in range(len(counts)))
            chi_squares[i, j] = chi_square
            if chi_square < min_chi_square:
                min_chi_square = chi_square
                min_h = h
                min_bins = bins
                best_L_eff = L_eff
                best_R = R
            # counts = counts / sum(counts)
            i += 1
        j += 1
        i = 0

    chi_squares = np.array(chi_squares)
    left_sigma, right_sigma, down_sigma, up_sigma = calc_chi_square_sigma(chi_squares, L_effs, Rs, len(counts),
                                                                          min_chi_square, best_L_eff, best_R, 1)
    if scatter_angle == 40:
        left = down = 0.081
        right = up = 0.018
    elif scatter_angle == 50:
        left = down = 0.031
        right = up = 0.052
    elif scatter_angle == 60:
        left = down = 0.055
        right = up = 0.052
    elif scatter_angle == 90:
        left = down = 0.037
        right = up = 0.068

    left_sigma = np.sqrt(left_sigma ** 2 + left ** 2)
    right_sigma = np.sqrt(right_sigma ** 2 + right ** 2)
    up_sigma = np.sqrt(up_sigma ** 2 + up ** 2)
    down_sigma = np.sqrt(down_sigma ** 2 + down ** 2)

    sigma_R = np.array([[down_sigma, up_sigma]]).T
    sigma_L_eff = np.array([[left_sigma, right_sigma]]).T

    left_2sigma, right_2sigma, down_2sigma, up_2sigma = calc_chi_square_sigma(chi_squares, L_effs, Rs, len(counts),
                                                                              min_chi_square, best_L_eff, best_R, 2)
    sigma2_R = np.array([[down_2sigma, up_2sigma]]).T
    sigma2_L_eff = np.array([[left_2sigma, right_2sigma]]).T

    plt.imshow(chi_squares.T, cmap='rainbow', origin='lower',
               extent=[min(L_effs), max(L_effs), min(Rs), max(Rs)],
               aspect='auto', vmax=vmax)

    plt.plot(best_L_eff, best_R, '.', color='black', markersize=10, label='min $\\chi^2$')
    plt.errorbar(best_L_eff, best_R, xerr=sigma_L_eff, yerr=sigma_R, color='black',
                 linestyle='')
    # plt.errorbar(best_L_eff, best_R, xerr=sigma2_L_eff, yerr=sigma2_R, color='yellow',
    #              linestyle='')
    plt.colorbar(label='$\\chi^2$')
    plt.ylabel('R', fontsize=20)
    plt.xlabel('$L_{eff}$', fontsize=20)
    plt.title("chi square map as a function of R and $L_{eff}$ " + " ," + str(
        scatter_angle) + " degree scatter angle",
              fontsize=20)
    plt.legend()
    plt.show()

    print('best $L_{eff}$ = ', best_L_eff, " + ", right_sigma, " - ", left_sigma)
    print("best R = ", best_R, " + ", up_sigma, " - ", down_sigma)
    print('min chi square = ', min_chi_square)

    plt.errorbar(num_photons, counts, yerr=errs, linestyle='', ecolor='black')
    plt.plot(num_photons, counts, '.', color='black', label='Creus et al. , data', )
    min_bins = np.array(min_bins)
    plt.errorbar(min_bins, min_h, yerr=np.sqrt(min_h), label='Aviv')
    plt.legend()
    plt.xlabel('integrated pulse height [p.e]', fontsize=20)
    plt.ylabel("entries/" + str(bin_size) + " p.e", fontsize=20)
    plt.title("number of photons absorbed in PMTs following " + str(
        scatter_angle) + " degree neutron scatter off LAr",
              fontsize=20)
    plt.show()


def is_in_arapuca(ypos, zpos, arapuca):
    if (ypos > arapuca[0][0] and ypos < arapuca[0][1]) and (zpos > arapuca[1][0] and zpos < arapuca[1][1]):
        return True
    else:
        return False


def is_in_arapucas(ypos, zpos):
    arapuca0 = [(-38, 37), (110, 194)]
    arapuca1 = [(118, 195), (38.7, 110)]
    arapuca2 = [(-118, -195), (-36.3, 38.7)]
    arapuca3 = [(37, 118), (-107, -36.3)]

    if is_in_arapuca(ypos, zpos, arapuca0) or is_in_arapuca(ypos, zpos, arapuca1) \
            or is_in_arapuca(ypos, zpos, arapuca2) or is_in_arapuca(ypos, zpos, arapuca3):
        return True
    else:
        return False


def step_from_dist(dists, bins, xlabel, ylabel, title, labels=None, colors=None, yscale=0):
    if type(dists) == list:
        for i in range(len(dists)):
            h, bins2 = np.histogram(dists[i], bins[i])
            plt.step((bins2[1:] + bins2[:-1]) / 2, h, color=colors[i], label=labels[i])
    else:
        h, bins2 = np.histogram(dists, bins)
        plt.step((bins2[1:] + bins2[:-1]) / 2, h)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if yscale:
        plt.yscale('log')
    plt.title(title)
    plt.legend()
    plt.show()


def plot_numPE(sim_pe, filename_data, filenames_errors=None, xs=None, xslabel=None, ys=None, xlim=None, ylim=None,
               label=None,
               gaussian_init_guess=None):
    if xlim is None:
        xlim = []
    exp_data = pd.read_csv(filename_data)
    pes = exp_data.numPE
    counts = exp_data.counts
    binsize = np.average([pes[i + 1] - pes[i] for i in range(len(pes) - 1)])
    I_tot = simpson(counts, dx=binsize)

    # sim_pe = sim_pe[sim_pe != 0]
    height, bin = np.histogram(sim_pe, np.arange(np.min(sim_pe), np.max(sim_pe), binsize))
    avg_bin = np.average([bin[i + 1] - bin[i] for i in range(len(bin) - 1)])
    area = simpson(height, dx=avg_bin)
    if area == 0 or len(height < 5):
        height, bin = np.histogram(sim_pe, np.arange(np.min(sim_pe), np.max(sim_pe), binsize / 2))
        avg_bin = np.average([bin[i + 1] - bin[i] for i in range(len(bin) - 1)])
        area = simpson(height, dx=avg_bin)
    yerror = np.sqrt(height) * I_tot / area
    height = height * I_tot / area
    norm = np.max(height) / np.max(counts)
    yerror = yerror / norm
    height = height / norm
    height = np.asarray(height)
    # print("chi square = ", chisquare(height, ys))

    plt.errorbar((bin[1:] + bin[:-1]) / 2, height, yerr=yerror, fmt='o', markersize=3,
                 label="Geant4 simulation")
    if filenames_errors:
        errors_down = pd.read_csv(filenames_errors[0])
        errors_up = pd.read_csv(filenames_errors[1])
        err_up = errors_up.counts - counts
        err_down = counts - errors_down.counts
        plt.errorbar(pes, counts, yerr=[err_down, err_up], fmt='o', markersize=3, )
    else:
        plt.plot(pes, counts, '.', label="data")
    plt.xlabel("S1 [pe]")
    plt.ylabel("counts")
    if label:
        plt.title(label)
    else:
        plt.title(filename_data.split("/")[-1].split("_")[0])
    if xlim:
        plt.xlim(xlim)
    else:
        dataMax = np.max(pes)
        simMax = np.max(sim_pe)
        plt.xlim([0, 1.2 * np.max([dataMax, simMax])])
    if ylim:
        plt.ylim(ylim)
    else:
        dataMax = np.max(counts)
        simMax = np.max(height)
        plt.ylim([0, np.max([dataMax, simMax]) * 1.2])
    if xs:
        for i in range(len(xs)):
            plt.plot(xs[i], ys, '--', label=xslabel[i])
    if gaussian_init_guess:
        bins_center = (bin[1:] + bin[:-1]) / 2
        parameters = curve_fit(gauss, bins_center, height,
                               p0=gaussian_init_guess)
        miu = parameters[0][2]
        H = parameters[0][0]
        A = parameters[0][1]
        d_miu = parameters[1][2][2]
        sigma = parameters[0][3]
        d_sigma = parameters[1][3][3]
        # print("miu = " + str(miu)+" +- " + str(d_miu) +", sigma = " + str(sigma))
        plt.plot(bins_center, gauss(bins_center, H, A, miu, sigma),
                 '--',
                 label='Fit result : $\\mu = $ ' + "{0:.2f}".format(miu) + ", $\\sigma = $ " + "{0:.2f}".format(sigma))

        parameters = curve_fit(gauss, pes, counts,
                               p0=gaussian_init_guess)
        miu = parameters[0][2]
        H = parameters[0][0]
        A = parameters[0][1]
        d_miu = parameters[1][2][2]
        sigma = parameters[0][3]
        d_sigma = parameters[1][3][3]
        # print("miu = " + str(miu) + " +- " + str(d_miu) + ", sigma = " + str(sigma))
        plt.plot(bins_center, gauss(bins_center, H, A, miu, sigma),
                 '--',
                 label='Data fit : $\\mu = $ ' + "{0:.2f}".format(miu) + ", $\\sigma = $ " + "{0:.2f}".format(sigma))
    plt.legend()
    plt.show()
    # plt.step((bin[1:] + bin[:-1]) / 2, height,
    #          label=filename.split("_")[2].split("_")[0] + " $\\frac{mm}{MeV}$ birks, "
    #                + filename.split("_")[3].split(".")[0] + " cm absLength")


def compare_gaussian(filenames, bounds, guesses, num_bins, parameter):
    mius = []
    sigmas = []
    for i in range(len(filenames)):
        dat = pd.read_csv(filenames[i])
        param = dat[parameter]
        param = param[param < bounds[i][1]]
        param = param[param > bounds[i][0]]
        h, bins = np.histogram(param, bins=num_bins)
        bins_center = (bins[1:] + bins[:-1]) / 2
        parameters = curve_fit(gauss, bins_center, h,
                               p0=guesses[i])
        miu = parameters[0][2]
        H = parameters[0][0]
        A = parameters[0][1]
        d_miu = parameters[1][2][2]
        sigma = parameters[0][3]
        mius.append(miu)
        sigmas.append(sigma)
        d_sigma = parameters[1][3][3]
        print("miu = " + str(miu) + " +- " + str(d_miu) + ", sigma = " + str(sigma))
        plt.step(bins_center, h, label='histogram ' + str(i))
        plt.plot(bins_center, gauss(bins_center, H, A, miu, sigma),
                 '--',
                 label='Fit result : $\\mu = $ ' + "{0:.2f}".format(miu) + ", $\\sigma = $ " + "{0:.2f}".format(sigma))
    plt.xlabel('$N_{PH}$', fontsize=16)
    plt.ylabel("count", fontsize=16)
    plt.title("number of photons produced by $\\gamma$ rays with different energies", fontsize=16)
    plt.legend(loc='best')
    plt.show()
    # np.save("../current/mius.npy", mius)
    # np.save("../current/sigmas.npy", sigmas)


def expo(a, miu, x):
    return a * np.exp(-x / miu)


def plot_prompt_fraction_graph(filename, particle, use_name=False):
    df = pd.read_csv(filename)
    #
    # times = df.tAbsorbed
    # events = df.eventID
    # print(times)

    # Parameters
    bin_width = 1  # nanoseconds
    prompt_window = (-20, 30)

    # Prepare result
    prompt_fractions = {}

    time_cap = 10000000
    # Group by event
    for event_id, group in df.groupby('eventID'):
        times = group['tAbsorbed'].values
        times = np.array(times)
        times = times[times < time_cap]

        # Create histogram
        min_time = times.min()
        max_time = times.max()
        if type(max_time) == str:
            prompt_fractions[event_id] = (np.nan, 0)
            continue
        bins = np.arange(min_time, max_time + bin_width, bin_width)
        hist, bin_edges = np.histogram(times, bins=bins)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

        # Find peak
        try:
            peak_index = np.argmax(hist)
        except ValueError:
            prompt_fractions[event_id] = (np.nan, 0)
            continue
        peak_time = bin_centers[peak_index]

        # Integration windows
        prompt_start = peak_time + prompt_window[0]
        prompt_end = peak_time + prompt_window[1]

        # Compute integrals
        in_prompt = (bin_centers >= prompt_start) & (bin_centers <= prompt_end)
        in_total = (bin_centers >= prompt_start)

        # Simpson integration
        if np.sum(in_prompt) >= 2 and np.sum(in_total) >= 2:  # Simpson requires at least 2 intervals
            numerator = simpson(hist[in_prompt], x=bin_centers[in_prompt])
            denominator = simpson(hist[in_total], x=bin_centers[in_total])
            pf = numerator / denominator if denominator > 0 else np.nan
            prompt_fractions[event_id] = (pf, denominator)
        else:
            pf = np.nan
            prompt_fractions[event_id] = (pf, 0)

    event_ids = []
    total_integral = []
    prompt_fraction_vals = []

    for event_id, group in df.groupby('eventID'):
        total = len(group)
        pf = prompt_fractions[event_id]
        if not np.isnan(pf[0]):
            event_ids.append(event_id)
            total_integral.append(pf[1])
            prompt_fraction_vals.append(pf[0])

    # 2D scatter plot
    fig = plt.figure(figsize=(8, 6))
    plt.scatter(prompt_fraction_vals, total_integral, alpha=0.6, color='navy', edgecolor='k')
    plt.ylabel("Total integral [a.u]",fontsize=20)
    plt.xlabel("Prompt Fraction",fontsize=20)
    plt.grid(True)
    plt.tight_layout()
    if use_name:
        yields = filename.split("/")[2].split('\\')[1]
        s = filename.split("/")[2].split('\\')[3].split("_")
        decay_times = filename.split("/")[2].split('\\')[2]
        name = "prompt_fraction_vs_tot_integral_" + yields + "_" + s[1] + "_" + s[2].split(".")[
            0] + "_" + decay_times + ".png"
        plt.title(
            "Prompt Fraction Histogram: " + particle + ", " + yields + ", " + decay_times + ", " + s[2].split(".")[
                0], fontsize=10)
    else:
        name = particle
        plt.title(name, fontsize=20)
    # print(name)
    fig.savefig(f"pf_bytype{particle}.png", dpi=300,bbox_inches='tight')
    plt.close(fig)
    return prompt_fraction_vals,total_integral
    # plt.show()


def plot_prompt_fraction_grid(file_list, particle_list, use_name=False):
    # Parameters
    bin_width = 10  # nanoseconds
    prompt_window = (-300, 500)

    import os
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from scipy.integrate import simpson

    n_files = len(file_list)
    n_cols = 2
    n_rows = (n_files + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 4 * n_rows))
    axes = axes.flatten()

    for i, filename in enumerate(file_list):
        df = pd.read_csv(filename)
        prompt_fractions = {}

        for event_id, group in df.groupby('eventID'):
            times = group['tAbsorbed'].values
            if isinstance(times.max(), str) or times.max() > 1e7:
                prompt_fractions[event_id] = (np.nan, 0)
                continue

            bins = np.arange(times.min(), times.max() + bin_width, bin_width)
            hist, bin_edges = np.histogram(times, bins=bins)
            bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

            try:
                peak_time = bin_centers[np.argmax(hist)]
            except ValueError:
                prompt_fractions[event_id] = (np.nan, 0)
                continue

            prompt_start = peak_time + prompt_window[0]
            prompt_end = peak_time + prompt_window[1]
            in_prompt = (bin_centers >= prompt_start) & (bin_centers <= prompt_end)
            in_total = (bin_centers >= prompt_start)

            if np.sum(in_prompt) >= 2 and np.sum(in_total) >= 2:
                numerator = simpson(hist[in_prompt], x=bin_centers[in_prompt])
                denominator = simpson(hist[in_total], x=bin_centers[in_total])
                pf = numerator / denominator if denominator > 0 else np.nan
                prompt_fractions[event_id] = (pf, denominator)
            else:
                prompt_fractions[event_id] = (np.nan, 0)

        prompt_vals = []
        total_integrals = []
        for eid, (pf, total) in prompt_fractions.items():
            if not np.isnan(pf):
                prompt_vals.append(pf)
                total_integrals.append(total)

        ax = axes[i]
        ax.scatter(prompt_vals, total_integrals, alpha=0.6, color='navy', edgecolor='k')
        if use_name:
            yields = filename.split("/")[2].split('\\')[1]
            decay_times = filename.split("/")[2].split('\\')[2]
            ax.set_title(
                "particle = " +
                os.path.basename(filename).split("_")[1] + ", det size = " + os.path.basename(filename).split("_")[
                    2].split(".")[0] + ", yield ratio = " + yields + ", decay times = " + decay_times,
                fontsize=8)
        else:
            side = len(filename.split("/")[5].split("_")) == 4
            # ax.set_title("particle = " + particle_list[i] + ", yield ratio = " + filename.split("/")[3] + ", decay times = " +
            #              filename.split("/")[4] + ", det size = " + filename.split("/")[5].split("_")[
            #                  2] + ", 1side = " + str(side), fontsize=7)
            ax.set_title(f"$\\{particle_list[i]}$" if particle_list[i]=="mu" else particle_list[i],fontsize=35)
        ax.set_xlabel("Prompt Fraction",fontsize=15)
        ax.set_ylabel("Total integral [a.u]",fontsize=15)
        # ax.set_xlim([0.4,0.7])
        # ax.set_ylim([0,sorted(total_integrals)[-7]])
        ax.grid(True)

    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])

    plt.suptitle("Prompt Fraction vs Total Integral for different particles", fontsize=14)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig("prompt_fraction_grid.png", dpi=300)
    plt.close(fig)

def discrete_count(distribution):
    freq = Counter(distribution)
    total = sum(freq.values())
    normalized_freq = {k: v / total for k, v in freq.items()}
    x = sorted(normalized_freq.keys())
    y = [normalized_freq[k] for k in x]
    return x,y

