from lowPr_highTa_onset import Ra_crit
import numpy as np

import matplotlib.pyplot as plt

Pr = 1e-8
Tas = 10**np.linspace(-5,5,num=100,endpoint=True)
Ras = list(Ra_crit(Ta,Pr) for Ta in Tas)

plt.plot(Tas, Ras)
plt.xscale('log')
plt.show()