import numpy as np
import matplotlib.pyplot as plt

#Pr = 1e-4 - 1e-10
#Ta = 1e5 - 1e14

def ra_crit_overstability(x, pr, ta):
    """ eqn 222 in Chapter 3 of Chandrasekhar 1961 """
    return 2 * (1 + pr) * (1 / x) * ( (1 + x)**3 + pr**2 * ta / (1 + pr)**2 )

def ra_crit_direct(x, ta):
    """ eqn 227 in Chapter 3 of Chandrasekhar 1961 """
    return (1 / x) * ( (1 + x)**3 + ta )

def x_star(ta, pr):
    """ eqn 225 in Chapter 3 of Chandrasekhar 1961 """
    return (ta * (1 - pr) / (1 + pr))**(1/3) - 1


Npr = 100
Nta = 100
Nx  = 100
pr_grid = np.logspace(-10, -4, Npr)
ta_grid = np.logspace(4, 14, Nta)
x_grid  = np.logspace(-1, 1, Nx)

current_ra_crits = np.zeros(Nx)

ra_crit = np.zeros((Npr, Nta))
for i, Pr in enumerate(pr_grid):
    for j, Ta in enumerate(ta_grid):
        x_s = x_star(Ta, Pr)
        for k, x in enumerate(x_grid):
            if x < x_s:
                current_ra_crits[k] = ra_crit_overstability(x, Pr, Ta)
            else:
                current_ra_crits[k] = ra_crit_direct(x, Ta)
        ra_crit[i,j] = (np.pi)**4 * np.min(current_ra_crits)

yy, xx = np.meshgrid(ta_grid, pr_grid)
pmesh = plt.pcolormesh(np.log10(xx), np.log10(yy), np.log10(ra_crit))
plt.xlabel('log10(Pr)')
plt.ylabel('log10(Ta)')
cbar = plt.colorbar(pmesh)
cbar.set_label('log10(Ra_crit)')
plt.show()
