##Calculates comoving volume of space
##Also Friedmann equation and various cosmological distances
##Needed for universe.py and mass_functioin.py

##Arkadiy Ukolov (2022)


import math
import numpy as np
from scipy.integrate import quad

pc = 3.086 * 10 ** (16) #m
G = 6.674 * 10 ** (-11) #m^3 kg^-1 s^-2
c = 299792458 #m s^-1
year = 365 * 24 * 60 * 60 #seconds

h = 0.6893 #0.7
hubble_0 = 100 * h #km s^-1 Mpc^-1
D_H = c / (hubble_0 * 1000) #Mpc

lower_limit_for_integrating_z = 0.00001

w_params = [
[-1.05,-0.04], #linear #-0.913,-0.2
[-0.878,-0.6], #CPL
[-1.11,0.43], #Barboza
[-0.912,-1.121], #LC
[-0.918,0.05], #Jassal
]

quint_model = [
'Linear-redshift',
'Chevallier-Polarski-Linder',
'Barboza-Alcaniz',
'Low Correlation',
'Jassal-Bagla-Padmanabhan',
]

def quintessense_w(a,w):
	integ = 0
	if w == 'lin':
		z = 1 / a - 1
		w = w_params[0][0] + w_params[0][1] * z
	
	return (1 + w) * a ** (-1) #careful here - w depends on a, not on z

def hubble_function(z_lim,w,omegas):
	omega_m = omegas[0]
	omega_k = omegas[1]
	omega_l = omegas[2]
	a = 1 / (1 + z_lim)
	a_lim = a

	dark_part = 0
	if omega_l != 0:
		dark_part = omega_l * np.exp(-3 * quad(quintessense_w, 1, a_lim, args = (w))[0])
	
	E = (omega_m * a ** (-3) + omega_k * a ** (-2) + dark_part) ** 0.5

	return E

def integrand(z_lim,w,omegas):
	return (1/hubble_function(z_lim,w,omegas))

def D_C(limiting_z, w, omegas): #line-of-sight co-moving distance == chi
	D_C = D_H * quad(integrand, lower_limit_for_integrating_z, limiting_z, args = (w,omegas))[0]
	return D_C

def D_M(limiting_z, w, omegas): #transverse co-moving distance == proper motion distance
	omega_k = omegas[1]
	if omega_k == 0:
		D_M = D_C(limiting_z, w, omegas)
	elif omega_k > 0:
		D_M = D_H * omega_k ** (-1/2) * np.sinh(omega_k ** (1/2) * D_C(limiting_z, w, omegas) / D_H)
	elif omega_k < 0:
		D_M = D_H * abs(omega_k) ** (-1/2) * np.sin((abs(omega_k)) ** (1/2) * D_C(limiting_z, w, omegas) / D_H)

	return D_M

def D_A(limiting_z, w, omegas): #angular diametre distance == apparent size distance
	D_A = D_M(limiting_z, w, omegas) / (1 + limiting_z)
	return D_A

def D_L(z,w,omegas): #luminosity distance at z
	D_L = (1+z) * D_M(z,w,omegas)
	return D_L

def dV(z, w, omegas): #comoving volume element at z
	dV = D_H * (1 + z) ** 2 * D_A(z, w, omegas) ** 2 / hubble_function(z,w,omegas)
	return dV

def comov_V(z, w, omegas): #calculate total comoving volume at z
	return quad(dV, lower_limit_for_integrating_z, z, args = (w, omegas))[0] * 4 * np.pi

def shell_V(V): #generate shells of comoving volumes per range of z
	shells = np.empty(len(V))
	j = len(V)-1
	while j > 0:
		shells[j] = V[j] - V[j-1] #V[j-1] - V[j]
		j -= 1

	return shells

