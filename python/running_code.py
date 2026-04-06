from cProfile import label
from collections import defaultdict
from fileinput import filename
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LogNorm
import os
import glob

from Code.functions import L_eff
from functions import *
import scipy as sc
from collections import Counter

import csv
import math
from matplotlib.patches import Patch

# -----------------------------
# Input files
# -----------------------------
dat1 = pd.read_csv("../current/scatterStats_35cm.csv")
dat2 = pd.read_csv("../current/scatterStats_25cm.csv")

datasets = {
    "35 cm": dat1,
    "25 cm": dat2
}

detectors = ["A0", "A1", "A2", "A3", "A4"]
categories = [
    "Single elastic only",
    "Single elastic + external scatter",
    "Multiple elastic"
]

colors = {
    "Single elastic only": "tab:blue",
    "Single elastic + external scatter": "tab:orange",
    "Multiple elastic": "tab:green"
}

hatches = {
    "35 cm": "",
    "25 cm": "//"
}

# -----------------------------
# Helper function
# -----------------------------
def compute_category_fractions(dat, detector_name=None):
    numElastic = dat.numElasticSensitive
    detector = dat.detector
    nonSensitiveScatter = dat[["ExtScatter", "CryoScatter", "innerLayerScatter"]].any(axis=1)

    if detector_name is None:
        mask = np.ones(len(dat), dtype=bool)
    else:
        mask = detector == detector_name

    total = mask.sum()
    if total == 0:
        return {
            "Single elastic only": 0.0,
            "Single elastic + external scatter": 0.0,
            "Multiple elastic": 0.0
        }

    single_elastic_only = ((numElastic == 1) & (~nonSensitiveScatter) & mask).sum() / total
    single_elastic_plus_nonSensitive = ((numElastic == 1) & (nonSensitiveScatter) & mask).sum() / total
    multiple_elastic = ((numElastic > 1) & mask).sum() / total

    return {
        "Single elastic only": single_elastic_only,
        "Single elastic + external scatter": single_elastic_plus_nonSensitive,
        "Multiple elastic": multiple_elastic
    }

# -----------------------------
# Plot 1: overall category fractions
# -----------------------------
def plot_overall_category_fractions(datasets):
    x = np.arange(len(categories))
    width = 0.32

    plt.figure(figsize=(8, 5))

    for geom_label in ["35 cm", "25 cm"]:
        vals = compute_category_fractions(datasets[geom_label])
        shift = -width / 2 if geom_label == "35 cm" else width / 2

        plt.bar(
            x + shift,
            [vals[c] for c in categories],
            width=width,
            color="tab:blue",
            hatch=hatches[geom_label],
            edgecolor="black"
        )

    plt.xticks(
        x,
        ["Single elastic", "Single elastic\n+ external scatter", "Multiple elastic"]
    )
    plt.ylabel("Fraction of events", fontsize=14)
    plt.title("Event category fractions", fontsize=16)
    plt.ylim(0, 1)
    plt.tick_params(axis="both", labelsize=12)

    legend_handles = [
        Patch(facecolor="tab:blue", edgecolor="black", hatch="", label="35 cm outer diameter"),
        Patch(facecolor="tab:blue", edgecolor="black", hatch="//", label="25 cm outer diameter")
    ]
    plt.legend(handles=legend_handles, title="Geometry")
    plt.tight_layout()
    plt.savefig("../current/event_types_comparison.png", dpi=300)
    plt.show()

# -----------------------------
# Plot 2: category fractions by detector
# -----------------------------
def plot_category_fractions_by_detector(datasets, detectors):
    x = np.arange(len(detectors))
    cat_offsets = [-0.25, 0.0, 0.25]
    geom_shift = 0.05
    bar_width = 0.09

    plt.figure(figsize=(11, 5))

    for j, category in enumerate(categories):
        base_pos = x + cat_offsets[j]

        for geom_label, shift_sign in [("35 cm", -1), ("25 cm", 1)]:
            vals = [
                compute_category_fractions(datasets[geom_label], detector_name=d)[category]
                for d in detectors
            ]

            plt.bar(
                base_pos + shift_sign * geom_shift,
                vals,
                width=bar_width,
                color=colors[category],
                hatch=hatches[geom_label],
                edgecolor="black"
            )

    plt.xticks(x, detectors)
    plt.ylabel("Fraction of events", fontsize=14)
    plt.xlabel("Detector", fontsize=14)
    plt.title("Elastic scattering categories by detector", fontsize=16)
    plt.ylim(0, 1)
    plt.tick_params(axis="both", labelsize=12)

    category_handles = [
        Patch(facecolor=colors[c], edgecolor="black", label=c)
        for c in categories
    ]
    geometry_handles = [
        Patch(facecolor="white", edgecolor="black", hatch=hatches[g], label=g)
        for g in ["35 cm", "25 cm"]
    ]

    leg1 = plt.legend(handles=category_handles, title="Event category", loc="upper left")
    plt.gca().add_artist(leg1)
    plt.legend(handles=geometry_handles, title="Outer diameter", loc="upper right")
    plt.savefig("../current/event_types_comparison_by_detector.png", dpi=300)
    plt.tight_layout()
    plt.show()

# -----------------------------
# Run both plots
# -----------------------------
plot_overall_category_fractions(datasets)
plot_category_fractions_by_detector(datasets, detectors)


################# ------------------------------------------------------------#######################


dat = pd.read_csv("../current/scatterStats_35cm.csv")
tof = dat.tof
# tof = tof[outerScatter > 3]
# tof = tof[cryoScatter == 0]
# tof = tof[innerLayerScatter == 0]
# tof = tof[nonSensitiveScatter > 0]

x = np.sort(tof)
y = np.arange(1, len(x) + 1) / len(x)   # fraction in [0,1]

plt.figure()
plt.plot(x, y)
plt.xlabel("TOF [ns]")  # change units if needed
plt.ylabel("CDF (fraction of events with TOF ≤ x)")
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()

n1 = 1.24  # LAr
n2 = 1.5  # glass

# thetas_i = np.linspace(0, np.pi / 2, 1000)
# thetas_t = np.arcsin((n1 / n2) * np.sin(thetas_i))
# cos_thetas_i = np.cos(thetas_i)
# cos_thetas_t = np.cos(thetas_t)
# R_s = ((n1 * cos_thetas_i - n2 * cos_thetas_t) / (n1 * cos_thetas_i + n2 * cos_thetas_t)) ** 2
# R_p = ((n2 * cos_thetas_i - n1 * cos_thetas_t) / (n2 * cos_thetas_i + n1 * cos_thetas_t)) ** 2
# R = (R_p + R_s) / 2
#
# thetas_i_deg = thetas_i * 180 / np.pi
#
# plt.plot(thetas_i_deg, R, '.')
# plt.xlabel('incidence angle [degrees]')
# plt.ylabel('Reflection coefficient (unpolarized light)')
# plt.title('Reflection vs incidence angle for LAr->glass')
# plt.show()


# geometry -> coverage
R = 5
h = 10
tile_size = 0.6
SiPM_A = tile_size ** 2
nphi = int((2 * np.pi * R) // tile_size)
A_row = nphi * SiPM_A
A_PMT = 2.54 ** 2
A_total = 2 * np.pi * R ** 2 + 2 * np.pi * R * h


def coverage_percent(top, bot, rows):
    return 100.0 * (A_PMT * (top + bot) + A_row * rows) / A_total


# extract (top, bot, rows) from filename parts like your help[1], help[3], help[5]
def parse_config(path):
    help = os.path.basename(path).split('_')
    top = int(help[1])
    bot = int(help[3])
    rows = int(help[5])
    return top, bot, rows


dir = '../current/money_plot/fullGeant4Sim/outerCellDiameter_35cm/'
num_subdirs = len(glob.glob(os.path.join(dir, "*/")))
labels = []
mean_pes = np.zeros((num_subdirs, 19))
std_pes = np.zeros((num_subdirs, 19))

coverages = []
means_1TPB = []
stds_1TPB = []
means_07TPB = []
stds_07TPB = []

i = -1
for dir2 in glob.glob(os.path.join(dir, "R=0.*5*/")):
    if len(glob.glob(os.path.join(dir2, '*'))) == 0:
        continue
    i += 1
    helper = os.path.basename(os.path.normpath(dir2))
    label = helper.split('PMT')[0]
    labels.append(label)

    # map config -> (mean, sem)
    cfg_to_stats = {}

    for file in glob.glob(os.path.join(dir2, "*.csv")):
        top, bot, rows = parse_config(file)

        data = pd.read_csv(file)
        npe = data.numPE.to_numpy()

        mean = np.mean(npe)
        s = np.std(npe, ddof=1)
        sem = s / np.sqrt(len(npe))

        cfg_to_stats[(top, bot, rows)] = (mean, sem)

    # sort configs by coverage (ascending)
    cfgs_sorted = sorted(cfg_to_stats.keys(),
                         key=lambda cfg: coverage_percent(*cfg))

    coverages = [coverage_percent(*cfg) for cfg in cfgs_sorted]
    means = [cfg_to_stats[cfg][0] for cfg in cfgs_sorted]
    sems = [cfg_to_stats[cfg][1] for cfg in cfgs_sorted]

    mean_pes[i] = means
    std_pes[i] = sems
    # if label == "R=0.3PMTR=0.4SiPM":
    #     means_1TPB,stds_1TPB = means,sems
    #     cvgs = coverages
    # if label == "R=0.3PMTR=0.4SiPM_0.7TPB":
    #     means_07TPB,stds_07TPB = means,sems
    # plt.errorbar(coverages, means, yerr=sems, fmt='o', label=labels[i])

# avg = (mean_pes[0] + mean_pes[1]) / 2
# print(avg)
#
# means_1TPB = np.array(means_1TPB)
# stds_1TPB = np.array(stds_1TPB)
# means_07TPB = np.array(means_07TPB)
# stds_07TPB = np.array(stds_07TPB)
# cvgs = np.array(cvgs)

# plt.xlabel('Coverage [%]', fontsize=20)
# plt.ylabel('Mean number of photoelectrons', fontsize=20)
# plt.title('Mean p.e. vs coverage for different PMT reflection coefficients', fontsize=20)
# plt.legend()
# plt.show()

avgs = (mean_pes[0] + mean_pes[1]) / 2
systematic_error = abs((mean_pes[0] - mean_pes[1])) / 2
statistical_error = 0.5 * np.sqrt(std_pes[0] ** 2 + std_pes[1] ** 2)
tot_error = np.sqrt(statistical_error ** 2 + systematic_error ** 2)
plt.errorbar(coverages[:8], avgs[:8], yerr=tot_error[:8], fmt='o', capthick=2,capsize=2, label='PMTs only')
plt.errorbar(coverages[8:-1], avgs[8:-1], yerr=tot_error[8:-1], fmt='s', capthick=2,capsize=2, label='PMTs + SiPM bands')
plt.errorbar(coverages[:-1], avgs[:-1], yerr=statistical_error[:-1], fmt='none', color='black',elinewidth=2,zorder=10,capsize=2)
plt.axhspan(2, 4, alpha=0.3,color='red',label='electronic noise region')
plt.axhspan(20, 30, alpha=0.3,color='purple',label='realistic trigger threshold region for $\\nu$ detectors')
plt.axhspan(40, 80, alpha=0.3,color='green',label='neutron identification threshold region')
plt.xticks(np.arange(0, 80, 10), fontsize=15)
plt.yticks(np.arange(0, 160, 20), fontsize=15)
plt.xlabel('Coverage [%]', fontsize=24)
plt.ylabel('Mean number of photoelectrons', fontsize=24)
# plt.title('Mean p.e. vs coverage ', fontsize=20)
plt.grid()
plt.legend(fontsize=11,loc='upper left')
plt.show()

# plt.errorbar(cvgs,means_1TPB,yerr=stds_1TPB, fmt='o', label="TPB QE = 1")
# plt.errorbar(cvgs,means_07TPB,yerr=stds_07TPB, fmt='o', label="TPB QE = 0.7")
# plt.xlabel('Coverage [%]', fontsize=20)
# plt.ylabel('Mean number of photoelectrons', fontsize=20)
# plt.title('Mean p.e. vs coverage for different TPB QE', fontsize=20)
# plt.legend()
# plt.show()
#
# ratios = means_07TPB / means_1TPB
# errors = np.sqrt((stds_07TPB/means_1TPB)**2 + (means_07TPB * stds_1TPB /(means_1TPB**2))**2)
# plt.errorbar(cvgs, ratios, yerr=errors,fmt = 'o',label='ratios')
# a,b = np.polyfit(cvgs, ratios, 1)
# xs = np.linspace(0,max(cvgs),100)
# ys = [0.7 for _ in range(100)]
# print(a,b)
# plt.plot(cvgs,a*cvgs+b,'--',label='linear fit')
# plt.plot(xs,ys,'-',label='expected value')
# plt.xlabel('Coverage [%]', fontsize=20)
# plt.ylabel('ratios', fontsize=20)
# plt.title('linearity of p.e. ratio between TPB efficiency = 1 and 0.7', fontsize=20)
# plt.legend()
# plt.show()


dat = pd.read_csv("../current/money_plot/InstrumentData/SiPM_QE.csv")
wavelength = dat.wavelength
energy = 1239.8 / wavelength
efficiency = dat.efficiency / 100
print(list(energy))
print(list(efficiency))
plt.plot(energy, efficiency, '.')
plt.xlabel('energy[eV]')
plt.ylabel('efficiency')
plt.title('SiPM QE vs photon energy')
plt.show()

dat = pd.read_csv("../current/money_plot/InstrumentData/PMT_QE.csv")
wavelength = dat.wavelength
energy = 1239.8 / wavelength
efficiency = dat.efficiency / 100
print(list(energy))
print(list(efficiency))
plt.plot(energy, efficiency, '.')
plt.xlabel('energy[eV]')
plt.ylabel('efficiency')
plt.title('R8520 PMT QE vs photon energy')
plt.show()

dat1 = pd.read_csv("../current/money_plot/InstrumentData/tpb_emission_spectrum.csv")
waveLength1 = dat1.waveLength
energy1 = 1239.8 / waveLength1
intensity = dat1.intensity
intensity = intensity / max(intensity)

dat2 = pd.read_csv("../current/money_plot/InstrumentData/ESR_reflectance.csv")
waveLength2 = dat2.waveLength
energy2 = 1239.8 / waveLength2
reflectance = dat2.reflectance
plt.plot(energy2, reflectance, '.', label='reflectance')

plt.plot(energy1, intensity, '.', label='webplotdigitizer extracted')
# print(list(energy2),list(reflectance))
# print(list(energy1), list(intensity))


eV = 1
x1 = [2.0664 * eV, 2.1014 * eV, 2.1377 * eV, 2.1752 * eV, 2.2140 * eV,
      2.2543 * eV, 2.2960 * eV, 2.3393 * eV, 2.3843 * eV, 2.4311 * eV,
      2.4797 * eV, 2.5303 * eV, 2.5830 * eV, 2.6380 * eV, 2.6953 * eV,
      2.7552 * eV, 2.8178 * eV, 2.8834 * eV, 2.9520 * eV, 3.0240 * eV,
      3.0996 * eV, 3.1791 * eV, 3.2627 * eV, 3.3509 * eV, 3.4440 * eV]
y1 = [0.2690, 0.3074, 0.2970, 0.2833, 0.3046,
      0.3235, 0.3047, 0.3183, 0.3429, 0.3931,
      0.4450, 0.5075, 0.5886, 0.6881, 0.8016,
      0.9149, 1.0000, 0.9872, 0.9217, 0.8213,
      0.5133, 0.1813, 0.1154, 0.1350, 0.1030]
x2 = [1.80 * eV, 2.00 * eV, 2.20 * eV, 2.40 * eV, 2.60 * eV, 2.80 * eV,
      2.90 * eV, 3.00 * eV, 3.05 * eV, 3.10 * eV, 3.20 * eV, 3.30 * eV,
      3.40 * eV, 3.50 * eV, 3.70 * eV, 4.00 * eV, 10.00 * eV]
y2 = [0.00, 0.02, 0.08, 0.20, 0.45, 0.80,
      0.92, 1.00, 0.98, 0.92, 0.80, 0.60,
      0.40, 0.22, 0.02, 0.00, 0.00]

plt.plot(x1, y1, '.', label='accurate')
plt.plot(x2, y2, '.', label='approximate')
plt.xlabel('energy [eV]')
plt.ylabel('emission intensity (relative)')
plt.legend()
plt.show()


def sat_exp(x, y0, A, k):
    return y0 + A * (1 - np.exp(-k * x))


def sat_mm(x, y0, A, x0):
    return y0 + A * (x / (x + x0))


def log_model(x, y0, a, x0):
    return y0 + a * np.log1p(x / x0)


def red_chi2(y, yfit, sigma, npar):
    r = (y - yfit) / sigma
    return np.sum(r * r) / max(len(y) - npar, 1)


xs, ys, labels, yerrs = np.load("../current/money_plot/Python_MC/mu_coverage_data.npy")
xs = np.asarray(xs, dtype=float)
ys = np.asarray(ys, dtype=float)
yerrs = np.asarray(yerrs, dtype=float)
plt.errorbar(xs, ys, yerr=yerrs, fmt='o-', label='tpb =1 run')

xs2, ys2, labels2, yerrs2 = np.load("../current/mu_coverage_data.npy")
xs2 = np.asarray(xs2, dtype=float)
ys2 = np.asarray(ys2, dtype=float)
yerrs2 = np.asarray(yerrs2, dtype=float)
plt.errorbar(xs2, ys2, yerr=yerrs2, fmt='o-', label='no tpb run (effectively tpb_efic = 1)')

# xmin, xmax = np.min(xs), np.max(xs)
# ymin, ymax = np.min(ys), np.max(ys)
# y0_guess = ys[np.argmin(xs)]
# A_guess = max(ymax - y0_guess, 1e-6)
# k_guess = 2.0 / max(xmax - xmin, 1e-6)  # “bends” over the x-range
# p0_exp = (y0_guess, A_guess, k_guess)
# p0_mm = (y0_guess, A_guess, 0.5 * (xmax - xmin))
# x0_guess = 0.1 * (xmax - xmin) if (xmax - xmin) > 0 else 1.0
# a_guess = (ymax - y0_guess) / max(np.log1p((xmax - xmin) / x0_guess), 1e-6)
# p0_log = (y0_guess, a_guess, x0_guess)
# models_parameters = []
# for func, p0 in zip([sat_exp, sat_mm, log_model], [p0_exp, p0_mm, p0_log]):
#     popt, pcov = curve_fit(
#         func, xs, ys,
#         p0=p0,
#         sigma=yerrs,  # optional
#         absolute_sigma=True if yerrs is not None else False,
#         maxfev=20000
#     )
#
#     perr = np.sqrt(np.diag(pcov))
#     models_parameters.append((popt, perr))
#     # print("params:", popt)
#     # print("errors:", perr)
#     yfit = func(xs, *popt)
#     xx = np.linspace(xmin, xmax, 300)
#     plt.plot(xx, func(xx, *popt),label=str(func).split(' ')[1])

# for x, y, lab in zip(xs, ys, labels):
#     if x < 18:
#         plt.annotate(lab, (x, y), textcoords="offset points", xytext=(15, -10), fontsize=8)
#     else:
#         plt.annotate(lab, (x, y), textcoords="offset points", xytext=(-45, -30), fontsize=8)

plt.xlabel("Instrumented surface coverage (%)", fontsize=15)
plt.ylabel("Mean number of photoelectrons per event", fontsize=15)
plt.legend(loc="best")
plt.grid(True, alpha=0.3)
plt.title("$N_{p.e}$ vs Optical Coverage for 2.5 MeV n single elastic scatter", fontsize=15)
plt.tight_layout()
plt.show()


# dat = pd.read_csv("TPB_lum_Temp.csv")
# temp = dat.temp
# lum = dat.lum
# plt.plot(temp,lum,'.')
# plt.xlabel('Temperature (K)')
# plt.ylabel('Photoluminesence')
# plt.show()

def calc_absLen(k, lamda):
    return lamda * 10 ** (-9) / (4 * np.pi * k)


ks = [8.2248, 6.0014, 4.2626, 3.2540, 2.5991, 2.1420, 1.7793, 1.5015, 1.1557]
absLens = []
lamdas = [826, 496, 354, 275, 225, 191, 165, 146, 124]
for i in range(len(ks)):
    absLens.append(calc_absLen(ks[i], lamdas[i]))

print(absLens)


def pmt_centers_in_circle(n, R_cm=5.0, s_cm=2.05, gap_cm=0.0, prefer="radial"):
    """
    Place up to n square PMTs (side s_cm) on the circular base of a cylinder (radius R_cm),
    ensuring each square lies fully inside the circle.

    PMT centers lie on a square grid with pitch = s_cm + gap_cm.

    prefer:
      - "radial": choose centers closest to (0,0) first (usually nicest)
      - "row": row-major fill (bottom->top, left->right)
    """
    if n < 0:
        raise ValueError("n must be >= 0")
    if n == 0:
        return []

    pitch = s_cm + gap_cm
    half_diag = (s_cm / 2.0) * math.sqrt(2.0)  # farthest corner distance from center

    # Conservative: require center to be within (R - half_diag)
    # so all corners remain inside circle.
    r_allow = R_cm - half_diag
    if r_allow <= 0:
        return []  # PMT is too big to fit at all

    # Build a grid that covers the allowable disk.
    # Max center coordinate magnitude that can still fit:
    max_c = r_allow
    m = int(math.floor(max_c / pitch))
    coords = [k * pitch for k in range(-m, m + 1)]

    candidates = [(x, y) for y in coords for x in coords
                  if (x * x + y * y) <= (r_allow * r_allow + 1e-12)]

    if prefer == "radial":
        candidates.sort(key=lambda p: p[0] * p[0] + p[1] * p[1])
    elif prefer == "row":
        # bottom->top, left->right (y then x)
        candidates.sort(key=lambda p: (p[1], p[0]))
    else:
        raise ValueError("prefer must be 'radial' or 'row'")

    return candidates[:n]


def max_pmts_in_circle(R_cm=5.0, s_cm=2.05, gap_cm=0.0, prefer="radial"):
    """
    Return (n_max, centers) for square PMTs fully contained in a circle.
    """
    pitch = s_cm + gap_cm
    half_diag = (s_cm / 2.0) * math.sqrt(2.0)
    r_allow = R_cm - half_diag
    if r_allow <= 0:
        return 0, []

    max_c = r_allow
    m = int(math.floor(max_c / pitch))
    coords = [k * pitch for k in range(-m, m + 1)]

    centers = [(x, y) for y in coords for x in coords
               if (x * x + y * y) <= (r_allow * r_allow + 1e-12)]

    if prefer == "radial":
        centers.sort(key=lambda p: p[0] * p[0] + p[1] * p[1])
    elif prefer == "row":
        centers.sort(key=lambda p: (p[1], p[0]))

    return len(centers), centers


# print(pmt_centers_in_circle(4,prefer="row"))
# n_max, centers = max_pmts_in_circle(prefer="row")
# print(n_max)
# print(centers[:10])  # first few (closest to center)

# ============================ GEOMETRY HELPERS ============================= #

def pmt_centers(n, s_cm=2.05):
    """
    Return square PMT centers (side = s_cm) packed around (0,0).
    Accepts n in 0..4 (0 => no PMTs).
    """
    if n == 0:  return []
    if n == 1:  return [(0.0, 0.0)]
    if n == 2:  return [(-s_cm / 2, 0.0), (s_cm / 2, 0.0)]
    if n == 3:  return [(-s_cm / 2, -s_cm / 2), (s_cm / 2, -s_cm / 2), (0.0, s_cm / 2)]
    if n == 4:  return [(-s_cm / 2, -s_cm / 2), (s_cm / 2, -s_cm / 2),
                        (-s_cm / 2, s_cm / 2), (s_cm / 2, s_cm / 2)]
    raise ValueError("pmts must be 0..4")


def in_pmt_square(u_cm, v_cm, cx_cm, cz_cm, s_cm=2.05):
    return (abs(u_cm - cx_cm) <= s_cm / 2) and (abs(v_cm - cz_cm) <= s_cm / 2)


def sipm_rows(R_cm=5.0, H_cm=10.0, tile_cm=0.6, rows=0):
    """
    Returns (row_y_list, nphi, dphi) for SiPM rows on the SIDE surface.
    Each row is a full ring centered at y = row_y_list[i].

    rows > 0           → number of rows, capped at floor(H_cm / tile_cm)
    """

    # Angular segmentation
    nphi = int((2 * np.pi * R_cm) // tile_cm)
    dphi = 2 * np.pi / max(nphi, 1)

    # Maximum number of rows that fit in height
    max_rows = int(np.floor(H_cm / tile_cm))

    # Cap rows safely
    rows = max(0, min(rows, max_rows))

    if rows == 0:
        row_y = []
    else:
        # Symmetric placement about y = 0
        offsets = np.arange(rows) - (rows - 1) / 2.0
        row_y = (offsets * tile_cm).tolist()

    return row_y, nphi, dphi


# for num in np.arange(0,17):
#     row,phi,dphi = sipm_rows(rows=num)
#     print(row)
# ========================= DETECTION / BRANCHING ========================== #

_rng = np.random.default_rng()


def consume_photon(hit_seq, cfg, tpb_eff=1.0,
                   R_cm=5.0, tile_cm=0.6,
                   R_side=0.95, R_pmt=0.15, R_sipm=0.40):
    """
    hit_seq: list of dicts sorted by time, each with:
      TOP/BOTTOM: {'surf':'TOP'|'BOTTOM','u':x_cm,'v':z_cm,'t':ns,'w':float}
      SIDE (legacy): {'surf':'SIDE','u':phi_rad,'v':y_cm,'t':ns,'w':float}
      SIDE (compact):{'surf':'SIDE','k':tile_index,'v':y_cm,'t':ns,'w':float}
    cfg: {'pmts_top':0..4, 'pmts_bot':0..4, 'sipm_rows':0..3}
    Returns (detected:bool, t_detect:float|None, weight:float)
    """
    pmts_top = pmt_centers_in_circle(cfg['pmts_top'])
    pmts_bot = pmt_centers_in_circle(cfg['pmts_bot'])
    row_y, nphi, dphi = sipm_rows(R_cm=R_cm, tile_cm=tile_cm, rows=cfg['sipm_rows'])
    half_y = tile_cm / 2.0

    for h in hit_seq:
        surf, t, w = h['surf'], h['t'], h.get('w', 1.0)
        PDE, R = 0.0, R_side

        if surf in ('TOP', 'BOTTOM'):
            # PDE = 0.266 * tpb_eff
            # R = R_pmt
            if pmts_top or pmts_bot:
                u, v = h.get('u'), h.get('v')  # x_cm, z_cm
                centers = pmts_top if surf == 'TOP' else pmts_bot
                for (cx, cz) in centers:
                    if in_pmt_square(u, v, cx, cz):
                        PDE = 0.266 * tpb_eff
                        R = R_pmt
                        break

        elif surf == 'SIDE':
            # PDE = 0.117 * tpb_eff
            # R = R_sipm
            # Full ring coverage per row: only y-band matters
            y_cm = h['v']
            in_band = any(abs(y_cm - y0) <= half_y for y0 in row_y)
            if in_band:
                PDE = 0.117 * tpb_eff
                R = R_sipm
            else:
                PDE = 0.0
                R = R_side

        # branch
        r = _rng.random()
        if r < PDE:
            return True, t, w
        elif r < PDE + R:
            continue
        else:
            return False, None, 0.0

    return False, None, 0.0


# ============================== CSV LOADERS =============================== #

def _header_indices(header):
    """Legacy photons CSV header → indices (case-insensitive)."""
    h = [c.strip() for c in header]

    def find(name, default=None):
        name = name.lower()
        for i, c in enumerate(h):
            if c.lower() == name: return i
        return default

    return {
        'event': find('eventID', find('event', 0)),
        'track': find('trackID', find('track', 1)),
        'surface': find('surface', 2),
        'u': find('u', 3),
        'v': find('v', 4),
        't': find('t_ns', find('t', None)),
        'w': find('weight', None),
    }


def load_hits_legacy(csv_path):
    """
    Legacy photons CSV with header:
      eventID,trackID,surface,u,v[,t_ns][,weight]
    Returns (photons_dict, sorted_event_ids).
    """
    photons, events = {}, set()
    with open(csv_path, 'r', newline='') as f:
        reader = csv.reader(f)
        header = next(reader)
        idx = _header_indices(header)
        for row in reader:
            if not row or row[0].startswith('#'): continue
            try:
                evt = int(row[idx['event']]);
                trk = int(row[idx['track']])
                surf = row[idx['surface']].strip().upper()
                u = float(row[idx['u']]);
                v = float(row[idx['v']])
                t = float(row[idx['t']]) if idx['t'] is not None else math.inf
                w = float(row[idx['w']]) if idx['w'] is not None else 1.0
            except Exception:
                continue
            photons.setdefault((evt, trk), []).append(
                {'surf': surf, 'u': u, 'v': v, 't': t, 'w': w}
            )
            events.add(evt)
    for hits in photons.values():
        hits.sort(key=lambda h: h['t'])
    return photons, sorted(events)


def load_hits_compact(csv_path, pos_scale=10, time_scale=100, nphi=52):
    """
    Compact photons CSV (no header): event,track,surf,u,v,t_10ps
      surf ∈ {'T','B','S'}; for S: u = tile index (0..nphi-1), v = y_mm
      ns = t_10ps / 100  (10 ps ticks) ;  cm = mm_ticks / 10  (pos_scale=10)
    """
    photons, events = {}, set()
    with open(csv_path, 'r') as f:
        for line in f:
            s = line.strip()
            if not s or s[0] == '#': continue
            parts = s.split(',')
            if len(parts) < 6: continue
            try:
                evt = int(parts[0])
                trk = int(parts[1])
                surfC = parts[2].strip().upper()
                surf = {'T': 'TOP', 'B': 'BOTTOM', 'S': 'SIDE'}.get(surfC, 'SIDE')
                u_i = int(parts[3])
                v_i = int(parts[4])
                t_i = int(parts[5])
            except Exception:
                continue

            if surf in ('TOP', 'BOTTOM'):
                u_cm = u_i / pos_scale
                v_cm = v_i / pos_scale
                hit = {'surf': surf, 'u': u_cm, 'v': v_cm,
                       't': t_i / time_scale, 'w': 1.0}
            else:  # SIDE
                y_cm = v_i / pos_scale
                hit = {'surf': 'SIDE', 'k': u_i, 'v': y_cm,
                       't': t_i / time_scale, 'w': 1.0}

            photons.setdefault((evt, trk), []).append(hit)
            # if evt not in events:
            #     print(evt)
            events.add(evt)

    for hits in photons.values():
        hits.sort(key=lambda h: h['t'])
    return photons, sorted(events)


# --------- Event meta CSV (eventID, detector, tOne, counters...) ---------- #

def _to_bool(s):
    if s is None: return False
    s = str(s).strip().lower()
    return s in ('1', 'true', 't', 'yes', 'y')


def load_event_meta(meta_csv_path):
    """
    Reads per-event CSV with header like:
      eventID,detector,tOne,numElasticSensitive,numInelasticSensitive,ExtScatter,numInelastic,good
    Returns dict: eventID -> {
        't1': float (ns),
        'ls': str (upper case detector name),
        'elastic_sens': int,
        'inelastic_sens': int,
        'ext_scatter': int,
        'inelastic': int,
        'good': bool
    }
    """
    out = {}
    with open(meta_csv_path, 'r', newline='') as f:
        reader = csv.reader(f)
        header = next(reader, None)
        if header is None:
            return out
        cols = [h.strip().lower() for h in header]

        def ix(name, default=None):
            try:
                return cols.index(name.lower())
            except ValueError:
                return default

        i_evt = ix('eventid', 0)
        i_det = ix('detector', 1)
        i_t1 = ix('tone', 2)  # supports 'tOne'
        i_el_s = ix('numelasticsensitive', None)
        i_in_s = ix('numinelasticsensitive', None)
        i_ext = ix('extscatter', None)
        i_in = ix('numinelastic', None)
        i_good = ix('good', None)

        # print("meta header cols:", cols)
        # print("i_evt,i_det,i_t1:", i_evt, i_det, i_t1)
        for row in reader:
            if not row or (row[0] and row[0].startswith('#')): continue
            try:
                evt = int(row[i_evt])
                ls = (row[i_det] if i_det is not None and i_det < len(row) else '').strip().upper()
                t1 = float(row[i_t1]) if (i_t1 is not None and i_t1 < len(row) and row[i_t1] != '') else math.inf
                # print(t1)
                elastic_sens = int(row[i_el_s]) if (
                        i_el_s is not None and i_el_s < len(row) and row[i_el_s] != '') else 0
                inelastic_sens = int(row[i_in_s]) if (
                        i_in_s is not None and i_in_s < len(row) and row[i_in_s] != '') else 0
                ext_scatter = int(row[i_ext]) if (i_ext is not None and i_ext < len(row) and row[i_ext] != '') else 0
                inelastic = int(row[i_in]) if (i_in is not None and i_in < len(row) and row[i_in] != '') else 0
                good = _to_bool(row[i_good]) if (i_good is not None and i_good < len(row)) else False
            except Exception:
                continue
            out[evt] = {
                't1': t1, 'ls': ls,
                'elastic_sens': elastic_sens,
                'inelastic_sens': inelastic_sens,
                'ext_scatter': ext_scatter,
                'inelastic': inelastic,
                'good': good
            }
            # k = next(iter(out))
            # print("SANITY meta_map sample:", k, out[k], "keys=", out[k].keys())
    return out


# ========================== SIM + HISTOGRAM CORE ========================== #

def npe_histogram(npe_arr, bins=50, range=None):
    counts, edges = np.histogram(npe_arr, bins=bins, range=range)
    centers = 0.5 * (edges[:-1] + edges[1:])
    return counts, edges, centers


def simulate_event_counts_with_tof(
        photons, event_ids, cfg, *,
        meta_map=None, tof_min=None, tof_max=None,
        ls_allow=None,
        require_good=None,  # None=ignore, True=only good, False=only not good
        min_elastic_sens=None,
        min_inelastic_sens=None,
        max_ext_scatter=None,
        max_inelastic=None,
        tpb_eff=1.0,
        R_cm=5.0, tile_cm=0.6,
        R_side=0.95, R_pmt=0.15, R_sipm=0.4
):
    """
    Per-config simulation with event meta filtering.
    Returns:
      kept_event_ids, npe_arr (aligned), t0_map, tof_map, ls_used
    """

    # print("simulate: meta_map is None?", meta_map is None)
    # if meta_map is not None:
    #     k = next(iter(meta_map))
    #     print("simulate: meta sample:", k, meta_map[k], "keys=", meta_map[k].keys())
    # normalize LS filter
    ls_set = None
    if ls_allow:
        ls_set = {s.strip().upper() for s in (ls_allow if isinstance(ls_allow, (list, tuple, set)) else [ls_allow])}

    npe = {}
    t0_det = {}
    tof = {}
    ls_used = {}
    for evt in event_ids:
        npe[evt] = 0.0
        t0_det[evt] = math.inf

    # Simulate inner-cell detection times and counts
    printed_events = set()
    for (evt, trk), hits in photons.items():
        if evt not in npe: continue

        if evt not in printed_events:
            print(cfg, evt)
            printed_events.add(evt)

        detected, tdet, w = consume_photon(
            hits, cfg, tpb_eff=tpb_eff,
            R_cm=R_cm, tile_cm=tile_cm,
            R_side=R_side, R_pmt=R_pmt, R_sipm=R_sipm
        )
        if detected:
            npe[evt] += w
            if tdet is not None and tdet < t0_det[evt]:
                t0_det[evt] = tdet

    kept = []
    for evt in event_ids:
        # must have inner detection
        if not np.isfinite(t0_det[evt]):
            # print('t0 isnt finite')
            continue

        rec = meta_map.get(evt) if meta_map is not None else None
        # print(rec)
        if rec is None:
            # print('rec is none')
            continue  # no LS info

        # LS name filter
        ls = rec.get('ls', '').upper()
        # print(ls)
        if ls_set is not None and ls not in ls_set:
            # print('ls not in ls set')
            continue

        # counters/flags filters
        if require_good is True and not rec.get('good', False): continue
        if require_good is False and rec.get('good', False):     continue

        if (min_elastic_sens is not None) and (rec.get('elastic_sens', 0) < min_elastic_sens): continue
        if (min_inelastic_sens is not None) and (rec.get('inelastic_sens', 0) < min_inelastic_sens): continue
        if (max_ext_scatter is not None) and (rec.get('ext_scatter', 0) > max_ext_scatter): continue
        if (max_inelastic is not None) and (rec.get('inelastic', 0) > max_inelastic): continue

        # TOF
        # print("evt:", evt, "rec type:", type(rec), "rec keys:", getattr(rec, "keys", lambda: "NO_KEYS")())
        # print("rec:", rec)
        t1 = rec.get('t1', math.inf)
        # print("t1 from rec.get:", t1)
        if not np.isfinite(t1):
            continue
        tof_evt = t1 - t0_det[evt]
        # print(tof_evt)
        if (tof_min is not None and tof_evt < tof_min) or (tof_max is not None and tof_evt > tof_max):
            continue

        kept.append(evt)
        tof[evt] = tof_evt
        ls_used[evt] = ls

    kept.sort()
    npe_arr = np.array([npe[e] for e in kept], dtype=float)
    return kept, npe_arr, t0_det, tof, ls_used


def _all_configs():
    """
    Ordered configs:
      1T, 2T, 3T, 4T,
      4T1B, 4T2B, 4T3B, 4T4B,
      4T4B+1row, 4T4B+2rows, 4T4B+3rows
    """
    out = []
    for t in (1, 2, 3, 4):
        out.append({'pmts_top': t, 'pmts_bot': 0, 'sipm_rows': 0})
    for b in (1, 2, 3, 4):
        out.append({'pmts_top': 4, 'pmts_bot': b, 'sipm_rows': 0})
    for rows in (1, 2, 4, 8, 16):
        out.append({'pmts_top': 4, 'pmts_bot': 4, 'sipm_rows': rows})
    return out


def _cfg_label(cfg):
    t = cfg['pmts_top']
    b = cfg['pmts_bot']
    r = cfg['sipm_rows']
    base = f"{t}T-{b}B"
    if r > 0: base += f"-{r}R"
    return base


def run_all_configs(
        csv_path, *,
        loader='compact',
        bins=50,  # (unused here; hist done later)
        seed=123,
        # geometry/encoding
        R_cm=5.0, H_cm=10.0, tile_cm=0.6, pos_scale=10, time_scale=100, nphi=52,
        # physics knobs
        tpb_eff=1.0, R_side=0.95, R_pmt=0.15, R_sipm=0.4,
        # Event meta & filters
        meta_csv_path=None,  # event meta CSV path (header as specified)
        ls_allow=None,  # e.g., "A0" or ["A0","A1"]
        require_good=None,  # None ignore, True require, False require NOT good
        min_elastic_sens=None,
        min_inelastic_sens=None,
        max_ext_scatter=None,
        max_inelastic=None,
        tof_min=None, tof_max=None  # ns
):
    """
    Loads hits once, then for each configuration:
      - simulates inner-cell detection (t0_det, npe)
      - reads event meta (detector, tOne, counters, good)
      - filters by LS name, counters, 'good', and TOF window
    Returns list of per-config dicts.
    """
    global _rng
    _rng = np.random.default_rng(seed)

    # Load photons
    if loader == 'legacy':
        photons, all_event_ids = load_hits_legacy(csv_path)
    else:
        photons, all_event_ids = load_hits_compact(csv_path, pos_scale, time_scale, nphi)

    # Load event meta
    meta_map = load_event_meta(meta_csv_path) if meta_csv_path is not None else None

    # coverage helpers
    A_bases = 2.0 * np.pi * (R_cm ** 2)
    A_side = 2.0 * np.pi * R_cm * H_cm
    A_total = A_bases + A_side
    A_pmt_1 = 2.54 * 2.54
    nphi_side = int((2 * np.pi * R_cm) // tile_cm)
    A_sipm_row = nphi_side * (tile_cm * tile_cm)

    def coverage(cfg):
        A_pmt = (cfg['pmts_top'] + cfg['pmts_bot']) * A_pmt_1
        A_sipm = cfg['sipm_rows'] * A_sipm_row
        return 100.0 * (A_pmt + A_sipm) / A_total

    results = []
    for i, cfg in enumerate(_all_configs()):
        _rng = np.random.default_rng(seed + i + 1)

        kept_ids, npe_arr, t0_map, tof_map, ls_used = simulate_event_counts_with_tof(
            photons, all_event_ids, cfg,
            meta_map=meta_map,
            ls_allow=ls_allow,
            require_good=require_good,
            min_elastic_sens=min_elastic_sens,
            min_inelastic_sens=min_inelastic_sens,
            max_ext_scatter=max_ext_scatter,
            max_inelastic=max_inelastic,
            tof_min=tof_min, tof_max=tof_max,
            tpb_eff=tpb_eff,
            R_cm=R_cm, tile_cm=tile_cm,
            R_side=R_side, R_pmt=R_pmt, R_sipm=R_sipm
        )

        results.append({
            'cfg': dict(cfg),
            'label': _cfg_label(cfg),
            'event_ids': kept_ids,  # events that PASSED all filters
            'npe': npe_arr,  # raw PE/event for kept events
            't0_map': t0_map,  # earliest inner-cell detected time (ns), all events
            'tof_map': tof_map,  # computed TOF (ns), only kept events
            'ls_used': ls_used,  # which LS was used, only kept events
            'coverage_pct': coverage(cfg)
        })
    return results


# ========================= FITTING / SUMMARY HELPERS ====================== #

def fit_gaussian_from_hist(x_centers, y_counts, min_count=1):
    """
    Fit y ≈ H + A * exp(-0.5*((x-μ)/σ)^2) using only bins with y>=min_count.
    Returns fit dict on success, or None on failure.
    """
    mask = (y_counts >= min_count)
    if mask.sum() < 3:
        return None

    xf = x_centers[mask]
    yf = y_counts[mask]

    A0 = float(yf.max())
    mu0 = float(np.average(xf, weights=yf))
    sigma0 = float(np.sqrt(np.average((xf - mu0) ** 2, weights=yf)))
    H0 = 0.0

    try:
        from scipy.optimize import curve_fit
        def gauss(x, H, A, mu, sigma):
            return H + A * np.exp(-0.5 * ((x - mu) / sigma) ** 2)

        popt, pcov = curve_fit(
            gauss, xf, yf,
            p0=[H0, A0, mu0, max(sigma0, 1e-6)],
            maxfev=20000
        )
        perr = np.sqrt(np.diag(pcov)) if pcov is not None else [np.nan] * 4
        H, A, mu, sigma = popt
        dH, dA, dmu, dsigma = perr
    except Exception:
        return None

    return {
        'H': float(H), 'A': float(A), 'mu': float(mu), 'sigma': abs(float(sigma)),
        'dH': float(dH), 'dA': float(dA), 'dmu': float(dmu), 'dsigma': float(dsigma),
        'x': xf, 'y': yf
    }


def _auto_integer_bins(all_results, bin_width=1):
    """
    Build integer-centered bins covering ALL configs: [-0.5, max+0.5] with given width.
    """
    max_npe = 0
    for r in all_results:
        if isinstance(r.get('npe'), np.ndarray) and r['npe'].size:
            m = int(np.nanmax(r['npe']))
            if m > max_npe: max_npe = m
    if max_npe < 1: max_npe = 1
    edges = np.arange(-0.5, max_npe + 0.5 + bin_width, bin_width, dtype=float)
    return edges


def summarize_configs(results, *, bins=None, min_count=1, trim_low=0.0, trim_high=0.99):
    """
    Re-bins each configuration, tries a Gaussian fit on trimmed data, and ALWAYS returns:
      mu_est, sigma_est (sample stats on trimmed data),
      sem_est (sigma_est / sqrt(n)),
      median_est, q16, q84 (central 68% interval).
    Trimming keeps only quantiles [trim_low, trim_high] of per-event NPE before stats.
    """

    def _auto_integer_bins(all_results, bin_width=1):
        max_npe = 0
        for r in all_results:
            if isinstance(r.get('npe'), np.ndarray) and r['npe'].size:
                m = int(np.nanmax(r['npe']))
                if m > max_npe: max_npe = m
        if max_npe < 1: max_npe = 1
        return np.arange(-0.5, max_npe + 0.5 + bin_width, bin_width, dtype=float)

    if bins is None or (isinstance(bins, (int, float)) and bins == 'auto'):
        bins = _auto_integer_bins(results, bin_width=1)

    out = []
    for res in results:
        npe = res['npe']
        if npe is None or len(npe) == 0:
            out.append({
                'label': res['label'], 'cfg': res['cfg'], 'coverage_pct': res['coverage_pct'],
                'counts': np.array([]), 'edges': np.array([]), 'centers': np.array([]),
                'fit': None,
                'mu_est': np.nan, 'sigma_est': np.nan, 'sem_est': np.nan,
                'median_est': np.nan, 'q16': np.nan, 'q84': np.nan,
                'n_events': 0, 'event_ids': res.get('event_ids', [])
            })
            continue

        # Trim outliers
        if 0.0 < trim_low or trim_high < 1.0:
            lo = np.quantile(npe, trim_low)
            hi = np.quantile(npe, trim_high)
            npe_use = npe[(npe >= lo) & (npe <= hi)]
            if len(npe_use) == 0:
                npe_use = npe
        else:
            npe_use = npe

        # Histogram for visualization/fit (optional)
        counts, edges = np.histogram(npe_use, bins=bins)
        centers = 0.5 * (edges[:-1] + edges[1:])

        # Try Gaussian fit on bins ≥ min_count (still optional)
        fit = None
        try:
            mask = (counts >= min_count)
            if mask.sum() >= 3:
                from scipy.optimize import curve_fit
                def gauss(x, H, A, mu, sigma):
                    return H + A * np.exp(-0.5 * ((x - mu) / sigma) ** 2)

                xf, yf = centers[mask], counts[mask]
                A0 = float(yf.max())
                mu0 = float(np.average(xf, weights=yf))
                sigma0 = float(np.sqrt(np.average((xf - mu0) ** 2, weights=yf)))
                H0 = 0.0
                popt, pcov = curve_fit(
                    gauss, xf, yf, p0=[H0, A0, mu0, max(sigma0, 1e-6)], maxfev=20000
                )
                perr = np.sqrt(np.diag(pcov)) if pcov is not None else [np.nan] * 4
                H, A, mu, sigma = popt
                dH, dA, dmu, dsigma = perr
                fit = {'H': float(H), 'A': float(A), 'mu': float(mu), 'sigma': abs(float(sigma)),
                       'dH': float(dH), 'dA': float(dA), 'dmu': float(dmu), 'dsigma': float(dsigma),
                       'x': xf, 'y': yf}
        except Exception:
            fit = None

        # Robust summaries on trimmed data
        mu_est = float(np.mean(npe_use))
        sigma_est = float(np.std(npe_use, ddof=1)) if len(npe_use) > 1 else 0.0
        sem_est = sigma_est / np.sqrt(len(npe_use)) if len(npe_use) > 0 else np.nan
        median_est = float(np.median(npe_use))
        q16, q84 = (float(np.quantile(npe_use, 0.16)),
                    float(np.quantile(npe_use, 0.84))) if len(npe_use) > 0 else (np.nan, np.nan)

        out.append({
            'label': res['label'], 'cfg': res['cfg'], 'coverage_pct': res['coverage_pct'],
            'counts': counts, 'edges': edges, 'centers': centers,
            'fit': fit,
            'mu_est': mu_est, 'sigma_est': sigma_est, 'sem_est': sem_est,
            'median_est': median_est, 'q16': q16, 'q84': q84,
            'n_events': int(len(npe)),
            'event_ids': res.get('event_ids', [])
        })
    return out


# =============================== PLOTTING ================================= #

def plot_mu_vs_coverage(summaries, yerr_mode='sem', title="Mean PE vs Coverage"):
    """
    Plot μ vs coverage with robust error bars.
    yerr_mode ∈ {'none','sigma','sem','central68','dmu'}
      - 'sigma'     : sample std (spread) on trimmed data
      - 'sem'       : std / sqrt(n)  (uncertainty on the mean)  ← recommended
      - 'central68' : (q84 - q16)/2  (half-width of central 68%)
      - 'dmu'       : fit uncertainty on μ (if fit available; else NaN)
      - 'none'      : no error bars
    """
    import matplotlib.pyplot as plt

    xs, ys, yerrs, labels = [], [], [], []
    for s in summaries:
        mu = s.get('mu_est', np.nan)
        if not np.isfinite(mu):
            continue
        xs.append(s['coverage_pct'])
        ys.append(mu)
        if yerr_mode == 'sigma':
            yerrs.append(s.get('sigma_est', np.nan))
        elif yerr_mode == 'sem':
            yerrs.append(s.get('sem_est', np.nan))
        elif yerr_mode == 'central68':
            q16 = s.get('q16', np.nan)
            q84 = s.get('q84', np.nan)
            yerrs.append((q84 - q16) / 2.0 if (np.isfinite(q16) and np.isfinite(q84)) else np.nan)
        elif yerr_mode == 'dmu':
            yerrs.append(s['fit'].get('dmu', np.nan) if s.get('fit') else np.nan)
        else:
            yerrs.append(None)
        labels.append(f"{s['label']} (n={s.get('n_events', 0)})")

    if not xs:
        print("[plot_mu_vs_coverage] Nothing to plot (no valid μ).")
        return

    plt.figure()
    if yerr_mode == 'none':
        plt.plot(xs, ys, 'o-')
    else:
        yerr_arr = np.array([ye if (ye is not None and np.isfinite(ye)) else np.nan for ye in yerrs])
        plt.errorbar(xs, ys, yerr=yerr_arr, fmt='.')

    for x, y, lab in zip(xs, ys, labels):
        plt.annotate(lab, (x, y), textcoords="offset points", xytext=(6, 6), fontsize=8)

    plt.xlabel("Instrumented surface coverage (%)")
    plt.ylabel("Mean photoelectrons per event (μ)")
    plt.grid(True, alpha=0.3)
    plt.title(title)
    plt.tight_layout()
    plt.show()

    return xs, ys, labels, yerr_arr


# ============================ EVENT-ID HELPERS ============================ #

def passed_event_ids(results, by_config=False, mode='union'):
    """
    by_config=False  → return a single sorted list
      mode='union'        : all events that passed in ANY config
      mode='intersection' : events that passed in EVERY config
    by_config=True   → return dict {label: [eventIDs...]} per config
    """
    if by_config:
        return {r['label']: list(r['event_ids']) for r in results}

    sets = [set(r['event_ids']) for r in results]
    if not sets:
        return []
    if mode == 'intersection':
        keep = set.intersection(*sets)
    else:
        keep = set.union(*sets)
    return sorted(keep)


def write_passed_event_ids(results, out_csv, by_config=True, mode='union'):
    """
    Save event IDs that passed the filters.
      by_config=True  → CSV columns: label,eventID (one row per (config,event))
      by_config=False → CSV column : eventID (union or intersection across configs)
    """
    with open(out_csv, 'w', newline='') as f:
        w = csv.writer(f)
        if by_config:
            w.writerow(['label', 'eventID'])
            for r in results:
                for e in r['event_ids']:
                    w.writerow([r['label'], e])
        else:
            ids = passed_event_ids(results, by_config=False, mode=mode)
            w.writerow(['eventID'])
            for e in ids:
                w.writerow([e])


def debug_configs(results, top_k=5):
    for r in results:
        npe = r['npe']
        print(f"{r['label']:>8} | kept={len(r['event_ids']):5d}", end='')
        if len(npe):
            print(
                f"  npe: mean={np.mean(npe):.1f}  std={np.std(npe, ddof=1):.1f}  min={np.min(npe):.0f}  max={np.max(npe):.0f}")
            if top_k and len(npe) > 0:
                idx = np.argsort(npe)[-top_k:]
                print("   top events:", [(int(r['event_ids'][i]), int(npe[i])) for i in idx])
        else:
            print("  (no events)")


res = run_all_configs(
    "../current/money_plot/Python_MC/mergedHits_tpb=1.csv",
    loader='compact',
    meta_csv_path="../current/money_plot/Python_MC/mergedStats_tpb=1.csv",
    ls_allow=["A0", "A1", "A2", "A3", "A4"],
    tof_min=40, tof_max=50,
    seed=42,
    R_cm=5.0, H_cm=10.0, tile_cm=0.6,
    tpb_eff=1.0
)

# Optional: quick debug to see tails/outliers
# debug_configs(res, top_k=3)

# Robust summary: trim top 1% to tame huge outliers
summaries = summarize_configs(res, bins=None, min_count=1, trim_low=0.0, trim_high=1)

# Plot with SEM (uncertainty on the mean) instead of full σ (spread)
xs, ys, labels, yerr_arr = plot_mu_vs_coverage(summaries, yerr_mode='sem')
arrs = [xs, ys, labels, yerr_arr]
# np.save("../current/money_plot/mu_coverage_data.npy", arrs)

# Event IDs that passed
per_cfg_ids = passed_event_ids(res, by_config=True)

#########################################################################################################

# angles = [25, 40, 50, 60, 90]
# dat = pd.read_csv("../current/numPE_150_cm_absLen_2.5_MeV_all.csv")
# pes = dat.numPE
# det = dat.detector
#
# A0_pes = pes[det == "A0"]
# A1_pes = pes[det == "A1"]
# A2_pes = pes[det == "A2"]
# A3_pes = pes[det == "A3"]
# A4_pes = pes[det == "A4"]
#
# all_pes = [A0_pes, A1_pes, A2_pes, A3_pes, A4_pes]
#
# for i in range(len(all_pes)):
#     h, bins = np.histogram(all_pes[i], 50)
#     avg_bins = (bins[1:] + bins[:-1]) / 2
#     avg_bins = avg_bins[h > 5]
#     h = h[h > 5]
#     guess = [2, 100, 20 * (i ** (2 / 3)), 10]
#     parameters = curve_fit(gauss, avg_bins[2:], h[2:],
#                            p0=guess)
#     miu = parameters[0][2]
#     H = parameters[0][0]
#     A = parameters[0][1]
#     d_miu = parameters[1][2][2]
#     sigma = parameters[0][3]
#     d_sigma = parameters[1][3][3]
#     print("miu = " + str(miu) + " +- " + str(d_miu) + ", sigma = " + str(sigma))
#     plt.plot(avg_bins, gauss(avg_bins, H, A, miu, sigma),
#              '--',
#              label='Fit result : $\\mu = $ ' + "{0:.2f}".format(miu) + ", $\\sigma = $ " + "{0:.2f}".format(sigma))
#     plt.step(avg_bins, h)
#     plt.xlabel('# p.e.')
#     plt.ylabel('count')
#     plt.title("# p.e., " + str(angles[i]) + " degree scattering angle")
#     plt.show()

# dat = pd.read_csv("../current/scatterStats_35cm.csv")
# tof = dat.tof
# det = dat.detector[tof > 0]
# tof = tof[tof > 0]
# A0_tof = tof[det == "A0"]
# A1_tof = tof[det == "A1"]
# A2_tof = tof[det == "A2"]
# A3_tof = tof[det == "A3"]
# A4_tof = tof[det == "A4"]
#
# tofs = [A0_tof,A1_tof,A2_tof,A3_tof,A4_tof]

#
# for i in range(len(tofs)):
#     tofs[i] = tofs[i][tofs[i] < 100]
#     h,bins = np.histogram(tofs[i],100)
#     avg_bins = (bins[1:] + bins[:-1])/2
#     guess = [10,50,50,10]
#     parameters = curve_fit(gauss, avg_bins, h,
#                            p0=guess)
#     miu = parameters[0][2]
#     H = parameters[0][0]
#     A = parameters[0][1]
#     d_miu = parameters[1][2][2]
#     sigma = parameters[0][3]
#     d_sigma = parameters[1][3][3]
#     print("miu = " + str(miu) + " +- " + str(d_miu) + ", sigma = " + str(sigma))
#     plt.plot(avg_bins, gauss(avg_bins, H, A, miu, sigma),
#              '--',
#              label='Fit result : $\\mu = $ ' + "{0:.2f}".format(miu) + ", $\\sigma = $ " + "{0:.2f}".format(sigma))
#     plt.step(avg_bins,h)
#     plt.xlabel('time of flight [ns]')
#     plt.ylabel('count')
#     plt.title("TOF elastic fiducial events, " + str(angles[i]) + " degree scattering angle")
#     plt.show()


# energies = [2.5, 14, 30]
# lengths = [33, 34, 35, 36, 37]

# with open("../current/batch_runs.mac", "w") as f:
#     f.write("/control/verbose 0\n/run/verbose 0\n/event/verbose 0\n/tracking/verbose 0\n/run/initialize\n\n")
#     for E in energies:
#         f.write(f"/gun/energy {E} MeV\n")
#         f.write("/run/beamOn 100000\n")
