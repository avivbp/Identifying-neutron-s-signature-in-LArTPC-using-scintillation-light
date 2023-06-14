import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from functions import *

dat = pd.read_csv('../neutronEscapeEnergy.csv')
escape_energy_keV = dat.neutronEscapeEnergy
print("number of events = ", len(escape_energy_keV))
h, bins = np.histogram(escape_energy_keV, 100)
plt.step((bins[1:] + bins[:-1]) / 2, h)
# plt.hist(escape_energy_keV, 100)
plt.xlabel('energy[keV]')
plt.ylabel('count')
plt.yscale('log')
plt.xticks(np.arange(0, 500, 20))
plt.title('5.3MeV neutron escape energy starting from center of DUNE')
plt.show()

escape_energy_keV = escape_energy_keV[escape_energy_keV < 70]
escape_energy_keV = escape_energy_keV[escape_energy_keV > 40]
print("number of escaped neutrons around 60keV = ", len(escape_energy_keV))


