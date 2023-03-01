import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.stats import poisson

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


def L_eff(light_yield_n, err_n, light_yield_el, err_el, energy):
    err = np.sqrt((err_n / light_yield_el) ** 2 + (light_yield_n * err_el / light_yield_el ** 2) ** 2)
    print("L_eff at E = ", energy, " keV is equal to ", light_yield_n / light_yield_el, " +- ", err)


def plot_pmts(filename, energy):
    dat = pd.read_csv(filename)
    num_up = np.array(dat.numUp)
    num_down = np.array(dat.numDown)
    tot_num_photons = np.array(dat.numPhotons)
    eDep = np.array(dat.eDep)

    num_up = num_up[num_up > 0]
    num_down = num_down[num_down > 0]
    tot_num_photons = tot_num_photons[tot_num_photons > 0]
    eDep = eDep[eDep > 0]

    h_up, bins_up = np.histogram(num_up, np.linspace(min(num_up), max(num_up), 100))
    plt.step((bins_up[1:] + bins_up[:-1]) / 2, h_up, label='number of photons reached top PMT')

    h_down, bins_down = np.histogram(num_down, np.linspace(min(num_down), max(num_down), 100))
    plt.step((bins_down[1:] + bins_down[:-1]) / 2, h_down, label='number of photons reached bottom PMT')

    plt.xlabel("number of scintillation photons")
    plt.ylabel('number of events')
    plt.title("number of scintillation photons absorbed in PMTs following " + str(energy) + "keV $\gamma$")
    plt.legend(loc='upper right', fontsize=5)
    plt.show()

    # print("mean number of photons absorbed in top PMT = ", np.mean(h_up), " +- ", np.std(h_up))
    # print("mean number of photons absorbed in bottom PMT = ", np.mean(h_down), " +- ", np.std(h_down))

    bins_up = bins_up[1:]
    bins_down = bins_down[1:]
    h_up = h_up[bins_up > 120]
    h_down = h_down[bins_down > 120]
    bins_up = bins_up[bins_up > 120]
    bins_down = bins_down[bins_down > 120]

    # fit with curve_fit
    # parameters = curve_fit(poisson_f, bins_up, h_up,
    #                        p0=180)
    # print(parameters)
    #
    # # plot poisson-deviation with fitted parameter
    # xs = 0.5 * (bins_up[1:] + bins_up[:-1])
    #
    # plt.plot(xs, poisson_f(xs, parameters[0][0]), marker='o', linestyle='', label='Fit result', )
    # plt.step((bins_up[1:] + bins_up[:-1]) / 2, h_up, label='number of photons reached top PMT')
    # plt.legend()
    # plt.show()

    print("number of scintillation photons absorbed in top PMT at photopeak / compton edge = ",
          bins_up[list(h_up).index(max(h_up))] * 1.11,
          "+-", bins_up[1] - bins_up[0])
    print("number of scintillation photons absorbed in bottom PMT at photopeak / compton edge = ",
          bins_down[list(h_down).index(max(h_down))] * 1.11,
          "+-", bins_down[1] - bins_down[0])
    #

    # h, bins = np.histogram(tot_num_photons, np.linspace(min(tot_num_photons), max(tot_num_photons), 80))
    # plt.step((bins[1:] + bins[:-1]) / 2, h, label='number of photons created')
    #
    h, bins = np.histogram(eDep, np.linspace(min(eDep), max(eDep), 80))
    plt.step((bins[1:] + bins[:-1]) / 2, h, label='normalized energy deposit in sensitive volume')
    #
    plt.xlabel("energy deposited [keV]")
    # plt.xlabel("number of photons")
    plt.ylabel("number of events")
    plt.title("energy deposited by " + str(energy) + "keV $\gamma$ [keV]")
    # plt.title("total number of scintillation photons created")
    plt.legend()
    plt.show()


def line(x, a, b):
    return a * x + b
