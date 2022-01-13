import mesa_reader as mr
import numpy as np
from pylab import *
from math import log10, pi

from scipy import interpolate
from scipy.interpolate import griddata

from numpy import loadtxt
import pandas as pd 
 
from constants import *
from plot_settings import *
from functions import *

def find_core(p):
	for i in range(len(p.mass)-1,1,-1):
		if p.gradR_core_composition[i] < p.gradA_core_composition[i]:
			return i
	return 1

mods = [1.5,2.0,3.0,5.0,10,12,14,20,30,40,50,60]

# Set up color cycle
color = plt.cm.viridis(np.linspace(0, 1, len(mods)))
mpl.rcParams['axes.prop_cycle'] = cycler('color', color)


# Data Location
FIGURES='./figures/' # Place to save plots
prefix = '/Users/ajermyn/Dropbox/Active_Projects/CBM_trends/output/runs/'
DIR = prefix + 'main_Z_MW_time_2022_01_07_15_39_00_sha_5262' + '/runs/' # The directory where you unpacked the data

fig = plt.figure(figsize=(6,5))
ax = plt.subplot(111)

for j in mods:
	p = mr.MesaData(DIR + str(j) + '/LOGS/profile_mid_MS.data')
	ind = find_core(p)
	m_core = p.mass[ind]

	sel = (p.mass < 1.2 * m_core) & (p.mass > 0.8 * m_core)

	ax.plot(p.mass[sel]/m_core, p.gradR_core_composition[sel]/p.gradA_core_composition[sel], label=r'$'+str(j)+'\,M_\odot$')

ax.set_ylim([0.9,1.2])
ax.set_xlim([0.8,1.2])
ax.legend(ncol=2)
ax.set_ylabel(r'$\nabla_{\rm rad}/\nabla_{\rm ad}$')
ax.set_xlabel(r'$m/m_{\rm boundary}$')
plt.tight_layout()
plt.savefig('figures/grad_profiles.pdf')