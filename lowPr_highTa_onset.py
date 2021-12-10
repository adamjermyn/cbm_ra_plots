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

def Ra_crit(Ta, Pr):
    Nx  = 200
    x_grid  = np.logspace(-2, 2, Nx)
    current_ra_crits = np.zeros(Nx)

    ra_crit = 0
    x_s = x_star(Ta, Pr)
    for k, x in enumerate(x_grid):
        if x < x_s:
            current_ra_crits[k] = ra_crit_overstability(x, Pr, Ta)
        else:
            current_ra_crits[k] = ra_crit_direct(x, Ta)
    return (np.pi)**4 * np.min(current_ra_crits)
