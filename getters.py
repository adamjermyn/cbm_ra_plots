import numpy as np

secyer = 3.14e7
Msun = 2e33

def original_m_core_over_m(h):
	return (h.m_core/Msun)/h.star_mass

def Ra_HI_getter(h):
	return np.log10(h.HI_Ra)
def Ra_HeI_getter(h):
	return np.log10(h.HeI_Ra)
def Ra_HeII_getter(h):
	return np.log10(h.HeII_Ra)
def Ra_FeCZ_getter(h):
	return np.log10(h.FeCZ_Ra)

def Re_HI_getter(h):
	return np.log10(h.HI_Re)
def Re_HeI_getter(h):
	return np.log10(h.HeI_Re)
def Re_HeII_getter(h):
	return np.log10(h.HeII_Re)
def Re_FeCZ_getter(h):
	return np.log10(h.FeCZ_Re)

def Pr_HI_getter(h):
	return np.log10(h.HI_Pr)
def Pr_HeI_getter(h):
	return np.log10(h.HeI_Pr)
def Pr_HeII_getter(h):
	return np.log10(h.HeII_Pr)
def Pr_FeCZ_getter(h):
	return np.log10(h.FeCZ_Pr)

def viscous_times_getter(h):
	return np.array([h.HI_dz**2/h.HI_nu,h.HeI_dz**2/h.HeI_nu,h.HeII_dz**2/h.HeII_nu,h.FeCZ_dz**2/h.FeCZ_nu]).T

def Ra_getter(h):
	Ras = np.array([h.HI_Ra,h.HeI_Ra,h.HeII_Ra,h.FeCZ_Ra])
	Ras[Ras == 0] = np.inf
	return np.log10(np.min(Ras,axis=0))
def Re_getter(h):
	Res = np.array([h.HI_Re,h.HeI_Re,h.HeII_Re,h.FeCZ_Re])
	Res[Res == 0] = np.inf
	return np.log10(np.min(Res,axis=0))
def Pr_getter(h):
	Prs = np.array([h.HI_Pr,h.HeI_Pr,h.HeII_Pr,h.FeCZ_Pr])
	Prs[Prs == 0] = np.inf
	return np.log10(np.min(Prs,axis=0))

def mdot_getter(h):
	return h.Mdot

def num_core_scale_heights(h):
	return np.log(10**h.log_center_Rho/h.rho_core_top)

def dlnm_core_getter(h):
	return h.dm_core*Msun/h.m_core

def dr_core_div_h_getter(h):
	return h.dr_core_div_h

def dlnr_core_getter(h):
	return h.dr_core/h.r_core

def rcore_div_h_getter(h):
	return h.r_core/h.hp_core_top

def m_core_over_m(h):
	return (h.m_core/Msun + h.dm_core)/h.star_mass

def tau_b_getter(h):
	Bs = np.array([h.b_HeII_max,h.b_HeI_max,h.b_FeCZ_max])
	ts = np.array([h.HeII_tau_eta, h.HeI_tau_eta, h.FeCZ_tau_eta]) / secyer

	tau = np.zeros(ts[0].shape)
	for i in range(len(tau)):
		j = np.argmax(Bs[:,i])
		tau[i] = ts[j,i]

	return tau

def tau_m_getter(h):
	Bs = np.array([h.b_HeII_max,h.b_HeI_max,h.b_FeCZ_max])
	ts = np.array([h.HeII_xm / h.Mdot, h.HeI_xm / h.Mdot, h.FeCZ_xm / h.Mdot])

	tau = np.zeros(ts[0].shape)
	for i in range(len(tau)):
		j = np.argmax(Bs[:,i])
		tau[i] = ts[j,i]

	return tau

def tau_m_div_tau_b_getter(h):
	Bs = np.array([h.b_HeII_max,h.b_HeI_max,h.b_FeCZ_max])
	tsm = np.array([h.HeII_xm / h.Mdot, h.HeI_xm / h.Mdot, h.FeCZ_xm / h.Mdot])
	tsb = np.array([h.HeII_tau_eta, h.HeI_tau_eta, h.FeCZ_tau_eta]) / secyer

	tau = np.zeros(tsm[0].shape)
	for i in range(len(tau)):
		j = np.argmax(Bs[:,i])
		tau[i] = tsm[j,i]/tsb[j,i]

	return tau

def tau_c_getter(h):
	Bs = np.array([h.b_HeII_max,h.b_HeI_max,h.b_FeCZ_max])
	Ts = np.array([h.turnover_HeII, h.turnover_HeI, h.turnover_FeCZ]) / secyer

	tau = np.zeros(Ts[0].shape)
	for i in range(len(tau)):
		j = np.argmax(Bs[:,i])
		tau[i] = Ts[j,i]

	return tau


def tau_buoy_getter(h):
	Bs = np.array([h.b_HeII_max,h.b_HeI_max,h.b_FeCZ_max])
	Ts = np.array([h.HeII_buoyant_time, h.HeI_buoyant_time, h.FeCZ_buoyant_time]) / secyer

	tau = np.zeros(Ts[0].shape)
	for i in range(len(tau)):
		j = np.argmax(Bs[:,i])
		tau[i] = Ts[j,i]

	return tau