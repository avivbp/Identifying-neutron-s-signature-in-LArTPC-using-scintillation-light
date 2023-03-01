import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from functions import *
from scipy.stats import lognorm
from scipy.stats import chisquare
from scipy.optimize import curve_fit

dat = pd.read_csv(
    "../Scintillation_Data\\nucleus_recoil_scint\\lightYield_30.csv")

tot_num_photons = dat.numPhotons
coincidence_time = dat.coincidenceTime
num_elastic_sensitive = dat.numElasticSensitive
num_up = dat.numUp
num_down = dat.numDown
scattered_inelastic = dat.scatteredInelastically
external_scatter = dat.scatteredNotSensitive

plt.hist(tot_num_photons, 100)
plt.xlabel('number of scintillation photons created')
plt.ylabel('number of events')
plt.title('number of scintillation photons created per event')
plt.show()

h1, bins1 = np.histogram(num_up, np.linspace(0, max(num_up), 100))
plt.step((bins1[1:] + bins1[:-1]) / 2, h1, label='number of photons reached top PMT')

h2, bins2 = np.histogram(num_down, np.linspace(0, max(num_down), 100))
plt.step((bins2[1:] + bins2[:-1]) / 2, h2, label='number of photons reached bottom PMT')

plt.xlabel('number of photons')
plt.ylabel('number of events')
plt.title('number of photons that reached bottom and top PMT per event')
plt.legend(loc='upper right', fontsize=8)
plt.show()

print("mean number of photons absorbed in top PMT = ", np.mean(h1), " +- ", np.std(h1))
print("mean number of photons absorbed in bottom PMT = ", np.mean(h2), " +- ", np.std(h2))

# popt = curve_fit(lognorm, bins1, h1, [0.9, 10])
# print(popt)
#
# fig, ax = plt.subplots(1, 1)
# s = 0.954
# mean, var, skew, kurt = lognorm.stats(s, moments='mvsk')
# x = np.linspace(lognorm.ppf(0.01, s),
#                 lognorm.ppf(0.99, s), 100)
# ax.plot(x, lognorm.pdf(x, s),
#         'r-', lw=5, alpha=0.6, label='lognorm pdf')
# plt.xlabel('x')
# plt.ylabel('f')
# plt.show()

# fit1 = lognorm((bins1[1:] + bins1[:-1]) / 2)
# fit2 = lognorm((bins2[1:] + bins2[:-1]) / 2)
#
# print(chisquare(h1, fit1))
# print(chisquare(h2, fit2))

# h, bins = np.histogram(np.array(num_up) / np.array(tot_num_photons) * 100, np.linspace(0, 15, 20))
# plt.step((bins[1:] + bins[:-1]) / 2, h, label='percentage of photons reached top PMT')
#
# h, bins = np.histogram(np.array(num_down) / np.array(tot_num_photons) * 100, np.linspace(0, 15, 20))
# plt.step((bins[1:] + bins[:-1]) / 2, h, label='percentage of photons reached bottom PMT')
#
# plt.xlabel('number of photons[%]')
# plt.ylabel('number of events')
# plt.title('percentage of photons that reached bottom and top PMT per event')
# plt.legend(loc='upper right', fontsize=8)
# plt.show()
#
# light_yield = (np.array(num_up) + np.array(num_down)) / np.array(tot_num_photons)
# plt.hist(light_yield * 100, 50)
# plt.xlabel('light yield[%]')
# plt.ylabel('number of events')
# plt.title('light yield per event at $25^{\circ}$ neutron scatter angle')
# plt.show()

# plt.hist(coincidence_time, 100)
# plt.xlabel('event coincidence time between PMTs[ns]')
# plt.ylabel('number of events')
# plt.title('coincidence time between PMTs per event')
# plt.show()

plot_pmts("../Scintillation_Data/electron_recoil_scint\\PMTS_60.csv", 60)


# #
# plot_pmts("C:\\Users\\Aviv\\PycharmProjects\\pythonProject3\\PMTS_122.csv", 122)
# #
# plot_pmts("C:\\Users\\Aviv\\PycharmProjects\\pythonProject3\\PMTS_511.csv", 511)
#
# plot_pmts("C:\\Users\\Aviv\\PycharmProjects\\pythonProject3\\PMTS_662.csv", 662)
#
# plot_pmts("C:\\Users\\Aviv\\PycharmProjects\\pythonProject3\\PMTS_1274.csv", 1274)

def avg_yield(a, da, b, db):
    print("avg light yield = ", (a + b) / 2, " +- ", np.sqrt(da ** 2 + db ** 2) / 2)


avg_yield(177.0057575757576, 2.9898989898989896, 182.74636363636364, 3.030303030303031)
avg_yield(389.8118181818182, 6.484848484848484, 386.739696969697, 5.98989898989899)
avg_yield(1035.7645454545457, 23.303030303030283, 1078.079090909091, 23.10101010101009)
avg_yield(1525.1512121212124, 26.4040404040404, 1503.1530303030304, 27.61616161616162)
avg_yield(3452.200909090909, 54.54545454545456, 3341.425151515152, 47.02020202020202)

electron_energies = np.array([60, 122, 341, 478, 1061])
light_yields = np.array(
    [179.8760606060606, 388.27575757575755, 1056.9218181818183, 1514.1521212121215, 3396.8130303030302])
# light_yields_err = [2.128510771951996, 4.413959384176372, 16.406461295713022, 19.103832940460904, 36.00731179076987]
light_yields_err = np.array([10, 20, 30, 40, 50])

plt.errorbar(electron_energies, light_yields, yerr=light_yields_err, fmt='.')
plt.plot(electron_energies, line(electron_energies, 3.198, -11.59321776))
plt.xlabel("Electron Energies [keV]")
plt.ylabel("Light Yield [a.u.]")
plt.title("Light yield as a function of electron recoil energy")
plt.show()

popt = curve_fit(line, electron_energies, light_yields, sigma=light_yields_err, p0=[1, 0])
print(popt)

a = popt[0][0]
da = np.sqrt(popt[1][0][0])
b = popt[0][1]
db = np.sqrt(popt[1][1][1])

print("line equation = ", a, "x", " + ", b)
L_eff(14, 1, line(16.4, a, b), da * 16.4 + db, 16.4)
