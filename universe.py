##Written for MPhys project at the University of Manchester
##Models number density of galaxies per redshift per square degree as a function of cosmology
##Supports cosmologies with evolving dark energy

##Arkadiy Ukolov (2022)


import sys
import math
import json
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.optimize import fmin
from astropy.io import fits
from dictionaries import *
from comov_vol import *


pc = 3.086 * 10 ** (16) #m
G = 6.674 * 10 ** (-11) #m^3 kg^-1 s^-2
c = 299792458 #m s^-1
year = 365 * 24 * 60 * 60 #seconds

h = 0.6893 #0.7
hubble_0 = 100 * h #km s^-1 Mpc^-1
D_H = c / (hubble_0 * 1000) #Mpc

deg_in_sky = 41257

lower_limit_for_integrating_z = 0.00001

crit_density = 3 * hubble_0 ** 2 / (8 * np.pi * G)

items = 0
z_upper_limit = 0.0
z = np.zeros(0)
z_incr = 0 #z[1] - z[0] #z[0] - z[1]


def N_func(n,shell): #turn volume into number of galaxies in a given shell
	return (n * shell) / (deg_in_sky * z_incr)

def schechter(M,z): #mass density function (galaxy number density)
	M_st = 10.88
	l_fi_1 = -3.17
	alpha_1 = -1.08

	diff = M - M_st
	fi_1 = 10 ** l_fi_1
	func = np.log(10) * np.exp(-10 ** diff) * 10 ** diff * (fi_1 * 10 ** (diff * alpha_1))

	return func

def integrate_schechter(z_step, w, omegas, lower_lim, upper_lim):
	r = (D_L(z_step, -1, MODELS['$\Lambda$CDM']['omegas']) / D_L(z_step, w, omegas)) ** 2
	M_1 = np.log10((10 ** lower_lim) / r) #lower_lim 

	M_2 = upper_lim
	n = quad(schechter, M_1, M_2, args = (z_step))[0]

	return n

def integral_of_age(z,w,omegas):
	function = 1 / ((1+z) * hubble_function(z,w,omegas) * hubble_0)
	return function

def age_of_epoch(z,w,omegas): #time that has passed since z
	integral = quad(integral_of_age, lower_limit_for_integrating_z, z, args = (w,omegas))[0] #in km^(-1) s Mpc
	integral = integral / 1000 * pc * 10 ** 6 / (year * 10 ** 9) #in Gyr

	return integral

def merger_function(z,time_difference): #works for >10^10 solar masses
	f_initial = 0.01
	m = 3
	f = f_initial * (1+z) ** m


	mu = 0.5
	a = 0.65 * mu
	b = 1.6
	tau = a * (1+z) ** b
	
	R = 0.7 * 10 ** (-3) * (1+z) ** (3.77)#merger rate per galaxy per Gyr

	return R

def subtract_mergers(N, N_mergers_add_stars, z, w, omegas, shells_V,lower_lim,upper_lim): #subtract merged galaxies as z --> 0
	mass_gain = 0.25
	ignore_galaxies = 0.5

	i = len(N) - 2#1
	size = len(N) - 1
	#N[i+1] = 0 #uncomment this for accumulated effect
	while i >= 0: #z[i+1]
		time_difference = age_of_epoch(z[size], w, omegas)-age_of_epoch(z[i], w, omegas) #in Gyr
		N_m = time_difference * quad(merger_function, z[i], z[size], args = (time_difference))[0]

		above_the_cut = np.log10(10 ** lower_lim + 10 ** lower_lim * mass_gain)

		#above the cut
		number_did_transition_above_cut = integrate_schechter(z[i], w, omegas, above_the_cut, upper_lim)
		number_above_cut = -N_m * N_func(number_did_transition_above_cut,shells_V[i])

		#inside the cut, ignore some galaxies
		number_did_transition_within_cut = integrate_schechter(z[i], w, omegas, lower_lim, above_the_cut)
		number_within_cut = -N_m * N_func(number_did_transition_within_cut,shells_V[i]) * ignore_galaxies

		#below the cut, added due to mass gain
		new_lower_lim = np.log10(10 ** lower_lim - 10 ** lower_lim * mass_gain)
		number_did_transition_below_cut = integrate_schechter(z[i], w, omegas, new_lower_lim, lower_lim)
		number_below_cut = N_m * N_func(number_did_transition_below_cut,shells_V[i])

		N_mergers_add_stars[i] = number_below_cut + N[i] #+ N[i] use N_mergers_add_stars[i+1] to see accumulated change
		N[i] = (number_above_cut + number_within_cut + number_below_cut) + N[i] #+ N[i] #use + N[i+1] for accumulated effect

		i -= 1

	return N, N_mergers_add_stars

def stellar_formation_function(z): #in solar masses per year PER GALAXY
	#f_0 = 10 #solar masses per year
	#m = 3
	#s = f_0 * (1+z) ** m

	#for 10.5 -- 10.75    #for the 10.75 -- 11 : should implement in bits

	#sd, 0, bin shift + 1; serr, 1; serr, 0;
	A = 0.34 #0.24 #0.66 #0.51
	alpha = 4.7 #5.73 #7.61 #4.27
	beta = -1.05 #-1.33 #-2.32 #-0.99
	k = -2.01 #-1.84 #-1.57 #-1.98
	log_s = (A * (1 + z) ** alpha) * np.exp(beta * (1 + z)) + k
	s = 10 ** log_s
	return s

def add_stars(shells_V, N, z, w, omegas, lower_lim, upper_lim):
	i = len(N) - 2
	size = len(N) - 1
	#N[i+1] = 0 #use this if want to see accumulated effect
	while i >= 0: #z[i+1] #FORM STARS FROM HIGH REDSHIFT DOWN TO LOW REDSHIFT
		time_difference = (age_of_epoch(z[size], w, omegas)-age_of_epoch(z[i], w, omegas)) * 10 ** 9 # in years
		stars_formed = quad(stellar_formation_function, z[i], z[size])[0] * time_difference  #z[i+1])[0] * time_difference #mass formed per galaxy

		if stars_formed > 10 ** lower_lim:
			new_lower_lim = 0
		else:
			new_lower_lim = np.log10(10 ** lower_lim - stars_formed)

		number_did_transition = integrate_schechter(z[i], w, omegas, new_lower_lim, lower_lim)
		N_transited = N_func(number_did_transition,shells_V[i])

		new_N = N_transited + N[i]# + N[i] #use N[i+1] if want to see accumulated effect
		N[i] = new_N

		i -= 1

	return N



def plot_schechter(des_z,low_lim,high_lim):
	fig = plt.figure()
	ax = fig.add_subplot(1, 1, 1)

	M_small = np.linspace(low_lim,high_lim,100)
	sch = schechter(M_small,des_z)

	ax.plot(M_small,sch, label = ('Sch at z = ' + str(des_z)))
	ax.set_xlabel('$\log{M/M_{\odot}}$')
	ax.set_yscale('log')
	ax.set_ylabel('$\phi$ $Mpc^{-3}$') 


def plot_mergers(z):
	fig = plt.figure()
	ax = fig.add_subplot(1, 1, 1)

	f_initial = 0.001
	m = 5.4
	f = f_initial * (1+z) ** m

	merg = merger_function(z)
	ax.plot(z,merg, label = ('Mergers'))
	ax.set_xlabel('z')
	ax.set_ylabel('Galaxy Major Merger Rate (Gyr^-1)') 
	ax.legend()


def plot_star_formation(z):
	fig = plt.figure()
	ax = fig.add_subplot(1, 1, 1)

	form = stellar_formation_function(z)
	ax.plot(z,form, label = ('Formed stars'))
	ax.set_xlabel('z')
	ax.set_ylabel('Stellar formation rate (Solar masses per year)') 
	ax.legend()


def plot_comoving_volume_and_number_density(models_to_try, lower_lim, upper_lim, axis_1, axis_2, INSTRUCTIONS):
	mergers = INSTRUCTIONS[0]
	stars = INSTRUCTIONS[1]
	shift = INSTRUCTIONS[2]
	normalise_by_instr = INSTRUCTIONS[3]
	use_data = INSTRUCTIONS[4]
	show_unchanged = INSTRUCTIONS[5]
	legend_on = INSTRUCTIONS[6]
	y_N_label = INSTRUCTIONS[7]
	shift_by_LCDM = INSTRUCTIONS[8]
	shift_by_unchanged = INSTRUCTIONS[9]
	save_to_json = INSTRUCTIONS[10]
	save_models_name = INSTRUCTIONS[11]

	massive_N = {}
	counter_of_models = 0
	V_LCDM = []
	V_LCDM_upper_lim = []
	V_LCDM_lower_lim = []
	N_LCDM_lower_lim = []
	N_LCDM_upper_lim = []
	LCDM = '$\Lambda$CDM'
	for model in models_to_try:
		omegas = MODELS[model]['omegas']
		w_list = MODELS[model]['w']
		linestyle_list = MODELS[model]['line']
		w_position = 0

		for w_name in w_list:
			V = []
			shells_V = []
			N = []
			i = 0

			if w_name in quint_model:
				w = 'lin'
			else:
				w = float(w_name)

			for z_step in z:
				V_step = comov_V(z_step,w,omegas)

				n = integrate_schechter(z_step, w, omegas, lower_lim, upper_lim)

				if model == LCDM:
					shell_up_V = shell_V(V_LCDM_upper_lim)
					shell_down_V = shell_V(V_LCDM_lower_lim)
				
				V.append(V_step)
				shells_V = shell_V(V)
				N.append(N_func(n,shells_V[i]))

				i += 1

			N_unchanged = np.asarray(N)
			N_stars_only = np.asarray(N)
			N_mergers_add_stars = np.zeros(len(N))

			if mergers:
				N,N_mergers_add_stars = subtract_mergers(N,N_mergers_add_stars,z,w,omegas,shells_V,lower_lim,upper_lim)

			N_without_stars = np.asarray(N)

			if stars:
				N = add_stars(shells_V, N, z, w, omegas, lower_lim, upper_lim)
				N_stars_only = add_stars(shells_V, N_stars_only, z, w, omegas, lower_lim, upper_lim)

			#normalise by LCDM
			N_norm = []
			N_unchanged_norm = []
			N_without_stars_norm = []
			N_stars_only_norm = []
			N_mergers_add_stars_norm = []

			if normalise_by_instr == 1:
				normalise_by = 1
			elif normalise_by_instr == 3:
				normalise_by = N[len(N)-1]

			for i in range(len(N)):
				if abs(normalise_by) >= 0:
					N_norm.append(N[i] / normalise_by)
					N_unchanged_norm.append(N_unchanged[i] / normalise_by)
					N_without_stars_norm.append(N_without_stars[i] / normalise_by)
					N_stars_only_norm.append(N_stars_only[i] / normalise_by)
					N_mergers_add_stars_norm.append(N_mergers_add_stars[i] / normalise_by)
				
				else:
					N_norm.append(1)
					N_unchanged_norm.append(1)
					N_without_stars_norm.append(1)
					N_stars_only_norm.append(1)
					N_mergers_add_stars_norm.append(1)

			N_norm = np.asarray(N_norm)
			N_without_stars_norm = np.asarray(N_without_stars_norm)
			N_unchanged_norm = np.asarray(N_unchanged_norm)
			N_stars_only_norm = np.asarray(N_stars_only_norm)
			N_mergers_add_stars_norm = np.asarray(N_mergers_add_stars_norm)
			N_not_shifted = np.asarray(N_norm)

			if shift:
				N_norm -= N_norm[len(N)-1]
				N_unchanged_norm -= N_unchanged_norm[len(N_unchanged_norm)-1]
				N_without_stars_norm -= N_without_stars_norm[len(N_without_stars_norm)-1] #comment this for accumulated change
				N_stars_only_norm -= N_stars_only_norm[len(N_stars_only_norm)-1] #comment this for accumulated change
				N_mergers_add_stars_norm -= N_mergers_add_stars_norm[len(N_mergers_add_stars_norm)-1]

			if shift_by_unchanged:
				N_norm -= N_unchanged_norm
				N_without_stars_norm -= N_unchanged_norm #comment this for accumulated change
				N_stars_only_norm -= N_unchanged_norm #comment this for accumulated change
				N_unchanged_norm -= N_unchanged_norm
				#N_mergers_add_stars_norm -= N_unchanged_norm

			if model == LCDM: #relies on LCDM going first in the list
				N_LCDM = N_norm.copy()

			if shift_by_LCDM:
				N_norm -= N_LCDM

				N_without_stars_norm -= N_LCDM #comment this for accumulated change
				N_unchanged_norm -= N_LCDM #comment this for accumulated change
				N_stars_only_norm -= N_LCDM
				N_mergers_add_stars_norm -= N_LCDM

			if save_to_json:
				counter_of_models += 1
				massive_N[model] = {}
				massive_N[model]['redshift'] = z.tolist()
				#massive_N[model]['w'] = w_name
				massive_N[model]['number_density'] = N_norm.tolist()

				if counter_of_models >= len(models_to_try):
					save_to_json = False
					json_file = open(save_models_name, 'w')
					json.dump(massive_N, json_file, indent = 4)
					json_file.close()
					print('MODELS SAVED')

			model_name = model + ', w = ' + str(w_name)
			if axis_1 == 'none':
				pass
			else:
				label = model
				axis_1.plot(z,V, label = label, linestyle = linestyle_list[w_position])
				axis_1.set_xlabel('Redshift')
				axis_1.set_ylabel('$V_c$ $Mpc^3$')
				axis_1.set_xlim([0, z_upper_limit-z_incr])
				axis_1.legend()

				print('V for ' + model_name + ' is plotted')

			if axis_2 == 'none':
				pass
			else:
				if model == LCDM:
					if not use_data:
						pass
					else:
						data_UDS = extract_useful_data('counts_deg2_binz025_UDS.dat',1)
						data_VISTA = extract_useful_data('counts_deg2_binz025_UltraVISTA.dat',1)
						data_VIDEO = extract_useful_data('counts_deg2_binz025_VIDEO.dat',1)

						data_UDS = normalise_data(data_UDS, shift, normalise_by, shift_by_LCDM, N_LCDM)
						data_VISTA = normalise_data(data_VISTA, shift, normalise_by, shift_by_LCDM, N_LCDM)
						data_VIDEO = normalise_data(data_VIDEO, shift, normalise_by, shift_by_LCDM, N_LCDM)

						average = average_of_datasets([data_UDS,data_VIDEO,data_VISTA])

						print('Chi sqrd for UDS ' + model, chi_sq_red(N_norm, data_UDS))
						print('Chi sqrd for VISTA ' + model, chi_sq_red(N_norm, data_VISTA))
						print('Chi sqrd for VIDEO ' + model, chi_sq_red(N_norm, data_VIDEO))
						print('Chi sqrd for AVERAGE ', chi_sq_red(N_norm, average))#[average,average,sigma]))

						axis_2.errorbar(average[:,0],average[:,1],yerr=average[:,2], xerr=0.125, fmt='ks', label = 'average')
						axis_2.errorbar(data_UDS[:,0],data_UDS[:,1],yerr=data_UDS[:,2], xerr=0.125, fmt='b*', ms = 7, label = 'UDS data')
						axis_2.errorbar(data_VISTA[:,0],data_VISTA[:,1],yerr=data_VISTA[:,2], xerr=0.125, fmt='rh', ms = 5, label = 'UltraVISTA data')
						axis_2.errorbar(data_VIDEO[:,0],data_VIDEO[:,1],yerr=data_VIDEO[:,2], xerr=0.125, fmt='gD', ms = 5, label = 'VIDEO data')

					if show_unchanged:
						axis_2.plot(z,N_unchanged_norm, label = 'Baseline', linewidth = 1.2, linestyle = (0, (5, 5, 1, 3, 1, 3)), c = 'g') #MODELS[model]['line'], c = 'g')
						axis_2.plot(z,N_without_stars_norm, label = 'Negative mergers', linewidth = 1.2, linestyle = (0, (10, 10)), c = 'r')
						axis_2.plot(z,N_stars_only_norm, label = 'Star formation', linewidth = 1.2, linestyle = (0, (7, 3)), c = 'k')
						axis_2.plot(z,N_mergers_add_stars_norm, label = 'Positive mergers', linewidth = 1.2, linestyle = (0, (10, 10)), c = 'b')

				label = model_name

				if model == LCDM:
					linewidth = 2
				elif model == '$\Lambda$CDM, different w':
					linewidth = 1.1
				else:
					linewidth = 1

				axis_2.plot(z,N_norm, label = label, linewidth = linewidth, linestyle = linestyle_list[w_position]) #comment for accumulated effect

				if model == LCDM:
					linestyle = ['dashdot', 'dashed']
				else:
					linestyle = [MODELS[model]['line'], MODELS[model]['line']]

				print(model_name, N_norm[len(N_norm)-1])

				#plot changes in ALL models:
				#axis_2.plot(z,N_not_shifted, label = (model_name + ' only stars'), linestyle = linestyle[0])
				#axis_2.plot(z,N_unchanged_norm, label = (model_name + ' unchanged'), linestyle = linestyle[0])
				#axis_2.plot(z,N_without_stars_norm, label = (model_name + ' mergers only'), linestyle = linestyle[1])

				axis_2.set_xlabel('Redshift', fontsize=15)
				axis_2.set_xlim([0, z_upper_limit])
				#ax2.set_yscale('log')
				axis_2.set_ylabel(y_N_label, fontsize=15)
				axis_2.tick_params(axis='both', which='major', labelsize=12)

				if legend_on:
					axis_2.legend(prop={"size":12})

				print('N for ' + model_name + ' is plotted')


			w_position += 1


def normalise_data(data_UDS, shift, normalise_by, shift_by_LCDM, N_LCDM):
	#data_UDS[:,1] = data_UDS[:,1] / (0.125 * 2) #by redshift bin
	#data_UDS[:,2] = data_UDS[:,2] / (0.125 * 2)

	B = data_UDS[len(data_UDS) - 1,1].copy()
	A_unc = data_UDS[:,2].copy()
	B_unc = data_UDS[len(data_UDS) - 1,2].copy()

	if normalise_by != 1: #divide by N[3]
		A = data_UDS[:,1].copy()
		#A_unc = data_UDS[:,2].copy() #uncomment this if putting normalise after shift
		#B_unc = data_UDS[len(data_UDS) - 1,2].copy()

		data_UDS[:,1] = A / B
		data_UDS[:,2] = np.sqrt((A_unc / A) ** 2 + (B_unc / B) ** 2) * data_UDS[:,1]

	if shift:
	#shift by N[3]
		A_unc = data_UDS[:,2].copy() ##comment this if putting normalise after shift
		B_unc = data_UDS[len(data_UDS) - 1,2].copy() #comment this if putting normalise after shift
		
		data_UDS[:,1] = data_UDS[:,1] - data_UDS[len(data_UDS) - 1,1] #- 1 #-1 here to avoide division by 0 in normalisation further on
		data_UDS[:,2] = np.sqrt(A_unc ** 2 + B_unc ** 2)

	if shift_by_LCDM:
		for i in range(len(data_UDS[:,1])):
			itemise = int(0.125 * 2 / z_incr * (i+1) - 0.125 / z_incr)
			if itemise < (z_upper_limit / z_incr):
				data_UDS[i,1] -= N_LCDM[itemise]
				#data_UDS[:,2] = np.sqrt(A_unc ** 2 + B_unc ** 2)

	return data_UDS

def chi_sq_red(N_norm, data_UDS):
	chi_sqrd = 0
	k = 1
	for data in data_UDS:
		itemise = int(0.125 * 2 / z_incr * k - 0.125 / z_incr) 
		if (itemise < (z_upper_limit / z_incr)):
			chi_sqrd += ((N_norm[itemise] - data[1])/data[2])**2
			k += 1

	red = chi_sqrd / (k - 2)
	return red

def average_of_datasets(datalist, diff_bins = False, red_start = 0, red_end = 0):
	if not diff_bins:
		average = np.empty((len(datalist[0][:,0]),3))
		upper_sum = np.zeros(len(datalist[0][:,0]))
		sigma_sqrd = np.zeros(len(datalist[0][:,0]))
		for j in range(len(datalist[0][:,0])): #so this goes down elements of a dataset
			for i in range(len(datalist)): #and this goes across datasets
				upper_sum[j] += datalist[i][j,1]/(datalist[i][j,2] ** 2)
				sigma_sqrd[j] += 1/(datalist[i][j,2] ** 2)

		for i in range(len(average[:])):
			average[i][0] = datalist[0][i,0]
			average[i][1] = upper_sum[i] / sigma_sqrd[i]
			average[i][2] = np.sqrt(1 / sigma_sqrd[i])

	else:
		maxim = 0
		minim = 1000
		for dataset in datalist:
			if len(dataset[0]) < minim: #because it is an array of ROWS and not columns...
				minim = len(dataset[0])
			if len(dataset[0]) > maxim:
				maxim = len(dataset[0])
		
		average_points = int((maxim + minim) / 2) #number of bins
		bin_size = (red_end - red_start) / average_points 

		average = np.empty((average_points,3))
		upper_sum = np.zeros(average_points)
		sigma_sqrd_inverted = np.zeros(average_points)
		for j in range(average_points): #so this goes down elements of average dataset
			for i in range(len(datalist)): #and this goes across STEPS
				step_data = datalist[i]
				for s in range(len(step_data)): #this goes across DATASETS in a STEP
					dataset_in_step = step_data[s]
					for k in range(len(dataset_in_step)): #this checks EACH point (row) of the ith dataset
						if dataset_in_step[k,0] <= bin_size * (j+1) and dataset_in_step[k,0] >= bin_size * j:
							sigma_sqrd_inverted[j] += 1 / (dataset_in_step[k,2] ** 2)
							upper_sum[j] += dataset_in_step[k,1] / (dataset_in_step[k,2] ** 2)

		for i in range(average_points):
			average[i][0] = i * bin_size + bin_size / 2 + red_start
			average[i][1] = upper_sum[i] / sigma_sqrd_inverted[i]
			average[i][2] = np.sqrt(1 / sigma_sqrd_inverted[i]) #* 2

	return average

def data_parsing(name,skip_lines):
	data_list = []
	data_line = []
	address = "data/" + name
	f = open(address, "r")
	file_text = f.read()
	f.close()
	lines = file_text.split('\n')

	i = 0
	for line in lines:
		if i < skip_lines:
			i += 1
		else:
			data_in_line = line.split(' ')
			data_line = []
			for data in data_in_line:
				if data != '':
					data_line.append(float(data))

			if data_line != []:
				data_list.append(np.asarray(data_line))

	data_list = np.asarray(data_list)
	return data_list


def extract_useful_data(name,skip_lines):
	data_list = data_parsing(name,skip_lines)
	length = int(z_upper_limit / 0.25)
	data = np.empty((length,3))

	for i in range(length):
		data[i][0] = (data_list[i][0] + data_list[i][1]) / 2
		data[i][1] = data_list[i][2]
		data[i][2] = data_list[i][3]

	return data



def main():
	global items, z_upper_limit, z, z_incr

	items = int(sys.argv[1])
	z_upper_limit = float(sys.argv[2])
	z = np.linspace(0,z_upper_limit,items)
	z_incr = z[1] - z[0] #z[0] - z[1]
	print('z step is ', round(z_incr,3))

	fig = plt.figure()
	ax1 = fig.add_subplot(1, 1, 1)

	fig_2 = plt.figure() #figsize=(8, 6)
	ax2 = fig_2.add_subplot(1,1,1)#(4, 2, 1)

	mergers = True 				#0
	stars = True 				#1
	shift = False				#2			False: N; True: N - N(z=3)
	normalise_by = 1			#3			1: nothing; 3:  N_unshifted(z=3)
	use_data = False 			#4
	show_unchanged = True 		#5
	legend_on = True 			#6
	label_y = 'N(z)' 			#7			#'$\Delta$N/N(z=' + str(z_upper_limit) + ')'  
	shift_by_LCDM = False		#8
	shift_by_unchanged = True 	#9
	save_models = False 		#10
	save_models_name = 'saved_models.json' #11

	INSTRUCTIONS = [mergers,stars,shift,normalise_by,use_data,show_unchanged,legend_on,label_y,shift_by_LCDM,shift_by_unchanged,save_models,save_models_name]

	models_to_try = ['$\Lambda$CDM', '$\Lambda$CDM, w = -1.15', '$\Lambda$CDM, w = -0.85']
	plot_comoving_volume_and_number_density(models_to_try, 10.75, 12, ax1, ax2, INSTRUCTIONS)

	plt.show()


def set_global_params(items_set, z_upper_limit_set):
	global items, z_upper_limit, z, z_incr

	items = items_set
	z_upper_limit = z_upper_limit_set
	z = np.linspace(0,z_upper_limit,items) #np.linspace(z_upper_limit,0,items)
	z_incr = z[1] - z[0] #z[0] - z[1]
	print('z step is ', round(z_incr,3))




#main() #comment this if launch code from other files





