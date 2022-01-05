import mesa_reader as mr
import numpy as np
from pylab import *
from math import log10, pi

from scipy import interpolate
from scipy.interpolate import griddata
from scipy.signal import argrelextrema

from numpy import loadtxt
import pandas as pd 
 
from constants import *
from plot_settings import *
from functions import *
from getters import *


p1 = mr.MesaData('work_2.4/LOGS/profile68.data') # Last profile, mid-MS
p2 = mr.MesaData('work_9/LOGS/profile91.data') # Last profile, mid-MS

# Find CZ boundaries
a = p1.gradr-p1.grada
bounds1 = np.where(np.sign(a[:-1]) != np.sign(a[1:]))[0] + 1
bounds1 = bounds1.astype(int)
if a[bounds1[0]-1] > 0:
	bounds1 = [0] + list(bounds1[:-1])
a = p2.gradr-p2.grada
bounds2 = np.where(np.sign(a[:-1]) != np.sign(a[1:]))[0] + 1
bounds2 = bounds2.astype(int)
if a[bounds2[0]-1] > 0:
	bounds2 = [0] + list(bounds2[:-1])

fig = plt.figure(figsize=(11,7))
spec = gridspec.GridSpec(ncols=2,nrows=3, figure=fig, width_ratios=(1,1), height_ratios=(1,1,1), wspace=0.05, hspace=0)
axes = list(list(fig.add_subplot(spec[i,j]) for j in range(2)) for i in range(3))

axes[0][0].plot(np.log10(p1.temperature), p1.gradr, label=r'$\nabla_{\rm rad}$')
axes[0][0].plot(np.log10(p1.temperature), p1.grada, label=r'$\nabla_{\rm ad}$')
[axes[0][0].axvspan(np.log10(p1.temperature[bounds1[2*i]]), np.log10(p1.temperature[bounds1[2*i+1]]), alpha=0.3) for i in range(len(bounds1)//2)]
axes[0][0].legend()

axes[1][0].plot(np.log10(p1.temperature), np.log10(p1.pressure_scale_height*rsun)-5, label=r'$\log h/\mathrm{km}$')
axes[1][0].plot(np.log10(p1.temperature), np.log10(p1.conv_vel), label=r'$\log v_{\rm c}/\mathrm{cm\,s^{-1}}$')
[axes[1][0].axvspan(np.log10(p1.temperature[bounds1[2*i]]), np.log10(p1.temperature[bounds1[2*i+1]]), alpha=0.3) for i in range(len(bounds1)//2)]
axes[1][0].text(4.72,6.5,'HeII')
axes[1][0].text(4.3,6.5,'HeI')
axes[1][0].text(4.1,6.5,'HI')
axes[1][0].legend(loc='upper left')

axes[2][0].plot(np.log10(p1.temperature), np.log10(p1.nu), label=r'$\log \nu/\mathrm{cm^2\,s^{-1}}$')
axes[2][0].plot(np.log10(p1.temperature), np.log10(p1.alpha)-10, label=r'$\log \alpha/10^{10}\mathrm{cm^2\,s^{-1}}$')
[axes[2][0].axvspan(np.log10(p1.temperature[bounds1[2*i]]), np.log10(p1.temperature[bounds1[2*i+1]]), alpha=0.3) for i in range(len(bounds1)//2)]
axes[2][0].legend(loc='upper left')

axes[0][0].set_title(r'$2.4M_\odot,\, X=0.36$')
axes[2][0].set_xlabel(r'$\log\, T/\mathrm{K}$')
axes[0][0].set_xlim([5.45,round(min(np.log10(p1.temperature)),2)])
axes[0][0].set_ylim([0.05,0.45])
axes[1][0].set_xlim([5.45,round(min(np.log10(p1.temperature)),2)])
axes[1][0].set_ylim([0,7.5])
axes[2][0].set_xlim([5.45,round(min(np.log10(p1.temperature)),2)])
axes[2][0].set_ylim([3,11])
axes[0][0].set_xticklabels([])
axes[1][0].set_xticklabels([])
print(min(np.log10(p1.temperature)))

axes[0][1].plot(np.log10(p2.temperature), p2.gradr, label=r'$\nabla_{\rm rad}$')
axes[0][1].plot(np.log10(p2.temperature), p2.grada, label=r'$\nabla_{\rm ad}$')
[axes[0][1].axvspan(np.log10(p2.temperature[bounds2[2*i]]), np.log10(p2.temperature[bounds2[2*i+1]]), alpha=0.3) for i in range(len(bounds2)//2)]

axes[1][1].plot(np.log10(p2.temperature), np.log10(p2.pressure_scale_height*rsun)-5, label=r'$\log h/\mathrm{km}$')
axes[1][1].plot(np.log10(p2.temperature), np.log10(p2.conv_vel), label=r'$\log v_{\rm c}/\mathrm{cm\,s^{-1}}$')
[axes[1][1].axvspan(np.log10(p2.temperature[bounds2[2*i]]), np.log10(p2.temperature[bounds2[2*i+1]]), alpha=0.3) for i in range(len(bounds2)//2)]
axes[1][1].text(4.65,6.5,'HeII')
axes[1][1].text(5.26,6.5,'Fe')

axes[2][1].plot(np.log10(p2.temperature), np.log10(p2.nu), label=r'$\log \nu/\mathrm{cm^2\,s^{-1}}$')
axes[2][1].plot(np.log10(p2.temperature), np.log10(p2.alpha)-10, label=r'$\log \alpha/10^{10}\mathrm{cm^2\,s^{-1}}$')
[axes[2][1].axvspan(np.log10(p2.temperature[bounds2[2*i]]), np.log10(p2.temperature[bounds2[2*i+1]]), alpha=0.3) for i in range(len(bounds2)//2)]

axes[0][1].set_title(r'$9M_\odot,\, X=0.36$')
axes[2][1].set_xlabel(r'$\log\, T/\mathrm{K}$')
axes[0][1].set_xlim([5.45,min(np.log10(p2.temperature))])
axes[0][1].set_ylim([0.05,0.45])
axes[1][1].set_xlim([5.45,min(np.log10(p2.temperature))])
axes[1][1].set_ylim([0,7.5])
axes[2][1].set_xlim([5.45,min(np.log10(p2.temperature))])
axes[2][1].set_ylim([3,11])
axes[0][1].set_xticklabels([])
axes[1][1].set_xticklabels([])
axes[0][1].yaxis.tick_right()
axes[1][1].yaxis.tick_right()
axes[2][1].yaxis.tick_right()

plt.savefig('figures/profile.pdf')