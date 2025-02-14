import numpy as np
from mpmath import mp
import auxiliary_functions as aux
import matplotlib.pyplot as plt
import os

mp.dps = 50
classical_clock_rate = 6 * mp.power(10,9) * 3600 * 24 * 365
dimensions = [20*i + 100 for i in range(46)]

resources = {
	'NVSieve_plain': {
		'physical_qubits_baseline': [],
		'physical_qubits_active': [],
		'time_baseline': [],
		'time_active': [],
		'reaction_limit': [],
		'hashing_time': [],
		'classical_time': []
	},
	'NVSieve_angular': {
		'physical_qubits_baseline': [],
		'physical_qubits_active': [],
		'time_baseline': [],
		'time_active': [],
		'reaction_limit': [],
		'hashing_time': [],
		'classical_time': []
	},
	'NVSieve_spherical': {
		'physical_qubits_baseline': [],
		'physical_qubits_active': [],
		'time_baseline': [],
		'time_active': [],
		'reaction_limit': [],
		'hashing_time': [],
		'classical_time': []
	},
	'NVSieve_filter': {
		'physical_qubits_baseline': [],
		'physical_qubits_active': [],
		'time_baseline': [],
		'time_active': [],
		'reaction_limit': [],
		'hashing_time': [],
		'classical_time': []
	},
	'GaussSieve_plain': {
		'physical_qubits_baseline': [],
		'physical_qubits_active': [],
		'time_baseline': [],
		'time_active': [],
		'reaction_limit': [],
		'hashing_time': [],
		'classical_time': []
	},
	'GaussSieve_angular': {
		'physical_qubits_baseline': [],
		'physical_qubits_active': [],
		'time_baseline': [],
		'time_active': [],
		'reaction_limit': [],
		'hashing_time': [],
		'classical_time': []
	},
	'GaussSieve_spherical': {
		'physical_qubits_baseline': [],
		'physical_qubits_active': [],
		'time_baseline': [],
		'time_active': [],
		'reaction_limit': [],
		'hashing_time': [],
		'classical_time': []
	},
	'GaussSieve_filter': {
		'physical_qubits_baseline': [],
		'physical_qubits_active': [],
		'time_baseline': [],
		'time_active': [],
		'reaction_limit': [],
		'hashing_time': [],
		'classical_time': []
	},
	'NVSieve_plain_no_qram': {
		'physical_qubits_baseline': [],
		'physical_qubits_active': [],
		'time_baseline': [],
		'time_active': [],
		'reaction_limit': []
	},
	'NVSieve_angular_no_qram': {
		'physical_qubits_baseline': [],
		'physical_qubits_active': [],
		'time_baseline': [],
		'time_active': [],
		'reaction_limit': []
	},
	'NVSieve_spherical_no_qram': {
		'physical_qubits_baseline': [],
		'physical_qubits_active': [],
		'time_baseline': [],
		'time_active': [],
		'reaction_limit': []
	},
	'NVSieve_filter_no_qram': {
		'physical_qubits_baseline': [],
		'physical_qubits_active': [],
		'time_baseline': [],
		'time_active': [],
		'reaction_limit': []
	},
	'GaussSieve_plain_no_qram': {
		'physical_qubits_baseline': [],
		'physical_qubits_active': [],
		'time_baseline': [],
		'time_active': [],
		'reaction_limit': []
	},
	'GaussSieve_angular_no_qram': {
		'physical_qubits_baseline': [],
		'physical_qubits_active': [],
		'time_baseline': [],
		'time_active': [],
		'reaction_limit': []
	},
	'GaussSieve_spherical_no_qram': {
		'physical_qubits_baseline': [],
		'physical_qubits_active': [],
		'time_baseline': [],
		'time_active': [],
		'reaction_limit': []
	},
	'GaussSieve_filter_no_qram': {
		'physical_qubits_baseline': [],
		'physical_qubits_active': [],
		'time_baseline': [],
		'time_active': [],
		'reaction_limit': []
	},
	'NVSieve_plain_nist': {
		'physical_qubits_baseline': [],
		'physical_qubits_active': [],
		'time_baseline': [],
		'time_active': [],
		'reaction_limit': []
	},
	'NVSieve_angular_nist': {
		'physical_qubits_baseline': [],
		'physical_qubits_active': [],
		'time_baseline': [],
		'time_active': [],
		'reaction_limit': []
	},
	'NVSieve_spherical_nist': {
		'physical_qubits_baseline': [],
		'physical_qubits_active': [],
		'time_baseline': [],
		'time_active': [],
		'reaction_limit': []
	},
	'NVSieve_filter_nist': {
		'physical_qubits_baseline': [],
		'physical_qubits_active': [],
		'time_baseline': [],
		'time_active': [],
		'reaction_limit': []
	},
	'GaussSieve_plain_nist': {
		'physical_qubits_baseline': [],
		'physical_qubits_active': [],
		'time_baseline': [],
		'time_active': [],
		'reaction_limit': []
	},
	'GaussSieve_angular_nist': {
		'physical_qubits_baseline': [],
		'physical_qubits_active': [],
		'time_baseline': [],
		'time_active': [],
		'reaction_limit': []
	},
	'GaussSieve_spherical_nist': {
		'physical_qubits_baseline': [],
		'physical_qubits_active': [],
		'time_baseline': [],
		'time_active': [],
		'reaction_limit': []
	},
	'GaussSieve_filter_nist': {
		'physical_qubits_baseline': [],
		'physical_qubits_active': [],
		'time_baseline': [],
		'time_active': [],
		'reaction_limit': []
	}
}

hashing_parameters = {
	'classical': {
		'NVSieve_angular': [],
		'NVSieve_spherical': [],
		'NVSieve_filter': [],
		'GaussSieve_angular': [],
		'GaussSieve_spherical': [],
		'GaussSieve_filter': []
	},
	'quantum': {
		'NVSieve_angular': [],
		'NVSieve_spherical': [],
		'NVSieve_filter': [],
		'GaussSieve_angular': [],
		'GaussSieve_spherical': [],
		'GaussSieve_filter': []
	}
}

List_size = {
	'NVSieve_plain': [],
	'NVSieve_angular': [],
	'NVSieve_spherical': [],
	'NVSieve_filter': [],
	'GaussSieve_plain': [],
	'GaussSieve_angular': [],
	'GaussSieve_spherical': [],
	'GaussSieve_filter': []
}

for case in ['classical', 'quantum']:
	with open('data/parameters_hashing_'+case+'.txt', "r") as data_file:
		for line in data_file:
			if 'NVSieve_angular' in line:
				string = 'NVSieve_angular'
			elif 'GaussSieve_angular' in line:
				string = 'GaussSieve_angular'
			elif 'NVSieve_spherical' in line:
				string = 'NVSieve_spherical'
			elif 'GaussSieve_spherical' in line:
				string = 'GaussSieve_spherical'
			elif 'NVSieve_filter' in line:
				string = 'NVSieve_filter'
			elif 'GaussSieve_filter' in line:
				string = 'GaussSieve_filter'
			elif (line.strip('\n').split(" "))[0].isnumeric():
				hashing_parameters[case][string].append( float((line.strip('\n').split(" "))[1]) )


for D in dimensions:

	# Computes the classical runtime

	L_NVSieve = D * mp.ceil(mp.exp(0.163*D + 0.102*mp.log(D) + 1.73))	# NVSieve list size
	S_NVSieve = mp.ceil(mp.exp(0.163*D + 0.102*mp.log(D) + 1.73))		# NVSieve list of centers
	L_GaussSieve = mp.ceil(mp.power(2, 0.193*D + 2.325))			# GaussSieve list size
	I_NVSieve = D								# Number of NVSieve iterations of sieving steps
	I_GaussSieve = mp.ceil(mp.power(2, 0.283*D + 0.335))			# Number of GaussSieve iterations of sieving steps
	
	searching_time = 3 * D * L_NVSieve * L_NVSieve / classical_clock_rate
	resources['NVSieve_plain']['classical_time'].append( searching_time )
	
	
	t = hashing_parameters['classical']['NVSieve_angular'][int((D-dimensions[0])/len(dimensions))]
	p2 = aux.p2_angular(t, D)
	
	hashing_time = 9 * aux.k_angular(t) * t * L_NVSieve / classical_clock_rate
	searching_time = 3 * D * L_NVSieve * L_NVSieve * p2 / classical_clock_rate
	resources['NVSieve_angular']['classical_time'].append( hashing_time + searching_time )
	
	
	t = hashing_parameters['classical']['NVSieve_spherical'][int((D-dimensions[0])/len(dimensions))]
	p2 = aux.p2_spherical(t, D)
	
	hashing_time = 5 * D * mp.ceil(mp.power(2, mp.sqrt(D))) * aux.k_spherical(t, D) * t * L_NVSieve / classical_clock_rate
	searching_time = 3 * D * L_NVSieve * L_NVSieve * p2 /classical_clock_rate
	resources['NVSieve_spherical']['classical_time'].append( hashing_time + searching_time  )
	
	
	alpha = hashing_parameters['classical']['NVSieve_filter'][int((D-dimensions[0])/len(dimensions))]
	cap_value = aux.Cap(alpha, D)
	t = aux.t_filter(alpha, D)
	
	hashing_time = 2 * mp.ceil(mp.log(D,2)) * cap_value * L_NVSieve * t / classical_clock_rate
	searching_time = 3 * D * mp.power(L_NVSieve * cap_value, 2) * t / classical_clock_rate
	resources['NVSieve_filter']['classical_time'].append( hashing_time + searching_time )
	
	
	searching_time = (125 * D - 19) * L_GaussSieve * I_GaussSieve / classical_clock_rate
	resources['GaussSieve_plain']['classical_time'].append( searching_time )
	
	
	t = hashing_parameters['classical']['GaussSieve_angular'][int((D-dimensions[0])/len(dimensions))]
	p2 = aux.p2_angular(t, D)
	
	hashing_time = 9 * aux.k_angular(t) * t * L_GaussSieve / classical_clock_rate
	searching_time = (125 * D - 19) * L_GaussSieve * I_GaussSieve * p2 / classical_clock_rate
	resources['GaussSieve_angular']['classical_time'].append( hashing_time + searching_time )
	
	
	t = hashing_parameters['classical']['GaussSieve_spherical'][int((D-dimensions[0])/len(dimensions))]
	p2 = aux.p2_spherical(t, D)
	
	hashing_time = 5 * D * mp.ceil(mp.power(2, mp.sqrt(D))) * aux.k_spherical(t, D) * t * L_GaussSieve / classical_clock_rate
	searching_time = (125 * D - 19) * L_GaussSieve * I_GaussSieve * p2 / classical_clock_rate
	resources['GaussSieve_spherical']['classical_time'].append( hashing_time + searching_time )
	
	
	alpha = hashing_parameters['classical']['GaussSieve_filter'][int((D-dimensions[0])/len(dimensions))]
	cap_value = aux.Cap(alpha, D)
	t = aux.t_filter(alpha, D)
	
	hashing_time = 2 * mp.ceil(mp.log(D,2)) * cap_value * L_GaussSieve * t / classical_clock_rate
	searching_time = (125 * D - 19) * L_GaussSieve * I_GaussSieve * mp.power(cap_value, 2) * t / classical_clock_rate
	resources['GaussSieve_filter']['classical_time'].append( hashing_time + searching_time )

##################################################

for D in dimensions:

	# Computes the classical hashing time for the quantum sieves and the list size C

	L_NVSieve = D * mp.ceil(mp.exp(0.163*D + 0.102*mp.log(D) + 1.73))	# NVSieve list size
	S_NVSieve = mp.ceil(mp.exp(0.163*D + 0.102*mp.log(D) + 1.73))		# NVSieve list of centers
	L_GaussSieve = mp.ceil(mp.power(2, 0.193*D + 2.325))			# GaussSieve list size
	I_NVSieve = D								# Number of NVSieve iterations of sieving steps
	I_GaussSieve = mp.ceil(mp.power(2, 0.283*D + 0.335))			# Number of GaussSieve iterations of sieving steps
	
	
	List_size['NVSieve_plain'].append( S_NVSieve )
	resources['NVSieve_plain']['hashing_time'].append( 0 )
	
	
	t = hashing_parameters['quantum']['NVSieve_angular'][int((D-dimensions[0])/len(dimensions))]
	p2 = aux.p2_angular(t, D)
	
	hashing_time = 9 * aux.k_angular(t) * t * L_NVSieve / classical_clock_rate
	List_size['NVSieve_angular'].append( S_NVSieve * p2 )
	resources['NVSieve_angular']['hashing_time'].append( hashing_time )
	
	
	t = hashing_parameters['quantum']['NVSieve_spherical'][int((D-dimensions[0])/len(dimensions))]
	p2 = aux.p2_spherical(t, D)
	
	hashing_time = 5 * D * mp.ceil(mp.power(2, mp.sqrt(D))) * aux.k_spherical(t, D) * t * L_NVSieve / classical_clock_rate
	List_size['NVSieve_spherical'].append( S_NVSieve * p2 )
	resources['NVSieve_spherical']['hashing_time'].append( hashing_time  )
	
	
	alpha = hashing_parameters['quantum']['NVSieve_filter'][int((D-dimensions[0])/len(dimensions))]
	cap_value = aux.Cap(alpha, D)
	t = aux.t_filter(alpha, D)
	
	hashing_time = 2 * mp.ceil(mp.log(D,2)) * cap_value * L_NVSieve * t / classical_clock_rate
	List_size['NVSieve_filter'].append( S_NVSieve * t * mp.power(cap_value, 2) )
	resources['NVSieve_filter']['hashing_time'].append( hashing_time )
	
	
	List_size['GaussSieve_plain'].append( L_GaussSieve )
	resources['GaussSieve_plain']['hashing_time'].append( 0 )
	
	
	t = hashing_parameters['quantum']['GaussSieve_angular'][int((D-dimensions[0])/len(dimensions))]
	p2 = aux.p2_angular(t, D)
	
	hashing_time = 9 * aux.k_angular(t) * t * L_GaussSieve / classical_clock_rate
	List_size['GaussSieve_angular'].append( L_GaussSieve * p2 )
	resources['GaussSieve_angular']['hashing_time'].append( hashing_time )
	
	
	t = hashing_parameters['quantum']['GaussSieve_spherical'][int((D-dimensions[0])/len(dimensions))]
	p2 = aux.p2_spherical(t, D)
	
	hashing_time = 5 * D * mp.ceil(mp.power(2, mp.sqrt(D))) * aux.k_spherical(t, D) * t * L_GaussSieve / classical_clock_rate
	List_size['GaussSieve_spherical'].append( L_GaussSieve * p2 )
	resources['GaussSieve_spherical']['hashing_time'].append( hashing_time )
	
	
	alpha = hashing_parameters['quantum']['GaussSieve_filter'][int((D-dimensions[0])/len(dimensions))]
	cap_value = aux.Cap(alpha, D)
	t = aux.t_filter(alpha, D)
	
	hashing_time = 2 * mp.ceil(mp.log(D,2)) * cap_value * L_GaussSieve * t / classical_clock_rate
	List_size['GaussSieve_filter'].append( L_GaussSieve * t * mp.power(cap_value, 2) )
	resources['GaussSieve_filter']['hashing_time'].append( hashing_time )

##################################################

for D in dimensions:

	# Given the list size C in List_size[sieve], computes the number of physical qubits and circuit time required by a baseline and an active-volume architectures

	L_NVSieve = D * mp.ceil(mp.exp(0.163*D + 0.102*mp.log(D) + 1.73))	# NVSieve list size
	S_NVSieve = mp.ceil(mp.exp(0.163*D + 0.102*mp.log(D) + 1.73))		# NVSieve list of centers
	L_GaussSieve = mp.ceil(mp.power(2, 0.193*D + 2.325))			# GaussSieve list size
	I_NVSieve = D								# Number of NVSieve iterations of sieving steps
	I_GaussSieve = mp.ceil(mp.power(2, 0.283*D + 0.335))			# Number of GaussSieve iterations of sieving steps
	
	for sieve in ['NVSieve_plain', 'NVSieve_angular', 'NVSieve_spherical', 'NVSieve_filter']:

		physical_qubits_baseline, physical_qubits_active, time_baseline, time_active, reaction_limit = aux.resources(D, List_size[sieve][int((D-dimensions[0])/len(dimensions))], 1, 'NVSieve')
		resources[sieve]['physical_qubits_baseline'].append(physical_qubits_baseline)
		resources[sieve]['physical_qubits_active'].append(physical_qubits_active)
		resources[sieve]['time_baseline'].append(I_NVSieve/2 * L_NVSieve * time_baseline)
		resources[sieve]['time_active'].append(I_NVSieve/2 * L_NVSieve * time_active)
		resources[sieve]['reaction_limit'].append(I_NVSieve/2 * L_NVSieve * reaction_limit)
		
		physical_qubits_baseline, physical_qubits_active, time_baseline, time_active, reaction_limit = aux.resources_no_qram(D, List_size[sieve][int((D-dimensions[0])/len(dimensions))], 1, 'NVSieve')
		resources[sieve+'_no_qram']['physical_qubits_baseline'].append(physical_qubits_baseline)
		resources[sieve+'_no_qram']['physical_qubits_active'].append(physical_qubits_active)
		resources[sieve+'_no_qram']['time_baseline'].append(I_NVSieve/2 * L_NVSieve * time_baseline)
		resources[sieve+'_no_qram']['time_active'].append(I_NVSieve/2 * L_NVSieve * time_active)
		resources[sieve+'_no_qram']['reaction_limit'].append(I_NVSieve/2 * L_NVSieve * reaction_limit)
		
		physical_qubits_baseline, physical_qubits_active, time_baseline, time_active, reaction_limit = aux.resources_nist(D, List_size[sieve][int((D-dimensions[0])/len(dimensions))], 1, mp.power(2, 40), 'NVSieve')
		resources[sieve+'_nist']['physical_qubits_baseline'].append(physical_qubits_baseline)
		resources[sieve+'_nist']['physical_qubits_active'].append(physical_qubits_active)
		resources[sieve+'_nist']['time_baseline'].append(I_NVSieve/2 * L_NVSieve * time_baseline)
		resources[sieve+'_nist']['time_active'].append(I_NVSieve/2 * L_NVSieve * time_active)
		resources[sieve+'_nist']['reaction_limit'].append(I_NVSieve/2 * L_NVSieve * reaction_limit)
		
		
	
	for sieve in ['GaussSieve_plain', 'GaussSieve_angular', 'GaussSieve_spherical', 'GaussSieve_filter']:	
	
		resources[sieve]['physical_qubits_baseline'].append(0)
		resources[sieve]['physical_qubits_active'].append(0)
		resources[sieve]['time_baseline'].append(0)
		resources[sieve]['time_active'].append(0)
		resources[sieve]['reaction_limit'].append(0)
	
		physical_qubits_baseline, physical_qubits_active, time_baseline, time_active, reaction_limit = aux.resources(D, List_size[sieve][int((D-dimensions[0])/len(dimensions))], 1, 'GaussSieve1')
		resources[sieve]['physical_qubits_baseline'][-1] = max(resources[sieve]['physical_qubits_baseline'][-1], physical_qubits_baseline)
		resources[sieve]['physical_qubits_active'][-1] = max(resources[sieve]['physical_qubits_active'][-1], physical_qubits_active)
		resources[sieve]['time_baseline'][-1] += 9 * I_GaussSieve * time_baseline
		resources[sieve]['time_active'][-1] += 9 * I_GaussSieve * time_active
		resources[sieve]['reaction_limit'][-1] += 9 * I_GaussSieve * reaction_limit
		
		physical_qubits_baseline, physical_qubits_active, time_baseline, time_active, reaction_limit = aux.resources(D, List_size[sieve][int((D-dimensions[0])/len(dimensions))], 0, 'GaussSieve1')
		resources[sieve]['physical_qubits_baseline'][-1] = max(resources[sieve]['physical_qubits_baseline'][-1], physical_qubits_baseline)
		resources[sieve]['physical_qubits_active'][-1] = max(resources[sieve]['physical_qubits_active'][-1], physical_qubits_active)
		resources[sieve]['time_baseline'][-1] += I_GaussSieve * time_baseline
		resources[sieve]['time_active'][-1] += I_GaussSieve * time_active
		resources[sieve]['reaction_limit'][-1] += I_GaussSieve * reaction_limit
		
		physical_qubits_baseline, physical_qubits_active, time_baseline, time_active, reaction_limit = aux.resources(D, List_size[sieve][int((D-dimensions[0])/len(dimensions))], 0, 'GaussSieve2')
		resources[sieve]['physical_qubits_baseline'][-1] = max(resources[sieve]['physical_qubits_baseline'][-1], physical_qubits_baseline)
		resources[sieve]['physical_qubits_active'][-1] = max(resources[sieve]['physical_qubits_active'][-1], physical_qubits_active)
		resources[sieve]['time_baseline'][-1] += I_GaussSieve * time_baseline
		resources[sieve]['time_active'][-1] += I_GaussSieve * time_active
		resources[sieve]['reaction_limit'][-1] += I_GaussSieve * reaction_limit
		
		
				
		resources[sieve+'_no_qram']['physical_qubits_baseline'].append(0)
		resources[sieve+'_no_qram']['physical_qubits_active'].append(0)
		resources[sieve+'_no_qram']['time_baseline'].append(0)
		resources[sieve+'_no_qram']['time_active'].append(0)
		resources[sieve+'_no_qram']['reaction_limit'].append(0)
	
		physical_qubits_baseline, physical_qubits_active, time_baseline, time_active, reaction_limit = aux.resources_no_qram(D, List_size[sieve][int((D-dimensions[0])/len(dimensions))], 1, 'GaussSieve1')
		resources[sieve+'_no_qram']['physical_qubits_baseline'][-1] = max(resources[sieve+'_no_qram']['physical_qubits_baseline'][-1], physical_qubits_baseline)
		resources[sieve+'_no_qram']['physical_qubits_active'][-1] = max(resources[sieve+'_no_qram']['physical_qubits_active'][-1], physical_qubits_active)
		resources[sieve+'_no_qram']['time_baseline'][-1] += 9 * I_GaussSieve * time_baseline
		resources[sieve+'_no_qram']['time_active'][-1] += 9 * I_GaussSieve * time_active
		resources[sieve+'_no_qram']['reaction_limit'][-1] += 9 * I_GaussSieve * reaction_limit
		
		physical_qubits_baseline, physical_qubits_active, time_baseline, time_active, reaction_limit = aux.resources_no_qram(D, List_size[sieve][int((D-dimensions[0])/len(dimensions))], 0, 'GaussSieve1')
		resources[sieve+'_no_qram']['physical_qubits_baseline'][-1] = max(resources[sieve+'_no_qram']['physical_qubits_baseline'][-1], physical_qubits_baseline)
		resources[sieve+'_no_qram']['physical_qubits_active'][-1] = max(resources[sieve+'_no_qram']['physical_qubits_active'][-1], physical_qubits_active)
		resources[sieve+'_no_qram']['time_baseline'][-1] += I_GaussSieve * time_baseline
		resources[sieve+'_no_qram']['time_active'][-1] += I_GaussSieve * time_active
		resources[sieve+'_no_qram']['reaction_limit'][-1] += I_GaussSieve * reaction_limit
		
		physical_qubits_baseline, physical_qubits_active, time_baseline, time_active, reaction_limit = aux.resources_no_qram(D, List_size[sieve][int((D-dimensions[0])/len(dimensions))], 0, 'GaussSieve2')
		resources[sieve+'_no_qram']['physical_qubits_baseline'][-1] = max(resources[sieve+'_no_qram']['physical_qubits_baseline'][-1], physical_qubits_baseline)
		resources[sieve+'_no_qram']['physical_qubits_active'][-1] = max(resources[sieve+'_no_qram']['physical_qubits_active'][-1], physical_qubits_active)
		resources[sieve+'_no_qram']['time_baseline'][-1] += I_GaussSieve * time_baseline
		resources[sieve+'_no_qram']['time_active'][-1] += I_GaussSieve * time_active
		resources[sieve+'_no_qram']['reaction_limit'][-1] += I_GaussSieve * reaction_limit
		
		
		
		resources[sieve+'_nist']['physical_qubits_baseline'].append(0)
		resources[sieve+'_nist']['physical_qubits_active'].append(0)
		resources[sieve+'_nist']['time_baseline'].append(0)
		resources[sieve+'_nist']['time_active'].append(0)
		resources[sieve+'_nist']['reaction_limit'].append(0)
	
		physical_qubits_baseline, physical_qubits_active, time_baseline, time_active, reaction_limit = aux.resources_nist(D, List_size[sieve][int((D-dimensions[0])/len(dimensions))], 1, mp.power(2, 40), 'GaussSieve1')
		resources[sieve+'_nist']['physical_qubits_baseline'][-1] = max(resources[sieve+'_nist']['physical_qubits_baseline'][-1], physical_qubits_baseline)
		resources[sieve+'_nist']['physical_qubits_active'][-1] = max(resources[sieve+'_nist']['physical_qubits_active'][-1], physical_qubits_active)
		resources[sieve+'_nist']['time_baseline'][-1] += 9 * I_GaussSieve * time_baseline
		resources[sieve+'_nist']['time_active'][-1] += 9 * I_GaussSieve * time_active
		resources[sieve+'_nist']['reaction_limit'][-1] += 9 * I_GaussSieve * reaction_limit
		
		physical_qubits_baseline, physical_qubits_active, time_baseline, time_active, reaction_limit = aux.resources_nist(D, List_size[sieve][int((D-dimensions[0])/len(dimensions))], 0, mp.power(2, 40), 'GaussSieve1')
		resources[sieve+'_nist']['physical_qubits_baseline'][-1] = max(resources[sieve+'_nist']['physical_qubits_baseline'][-1], physical_qubits_baseline)
		resources[sieve+'_nist']['physical_qubits_active'][-1] = max(resources[sieve+'_nist']['physical_qubits_active'][-1], physical_qubits_active)
		resources[sieve+'_nist']['time_baseline'][-1] += I_GaussSieve * time_baseline
		resources[sieve+'_nist']['time_active'][-1] += I_GaussSieve * time_active
		resources[sieve+'_nist']['reaction_limit'][-1] += I_GaussSieve * reaction_limit
		
		physical_qubits_baseline, physical_qubits_active, time_baseline, time_active, reaction_limit = aux.resources_nist(D, List_size[sieve][int((D-dimensions[0])/len(dimensions))], 0, mp.power(2, 40), 'GaussSieve2')
		resources[sieve+'_nist']['physical_qubits_baseline'][-1] = max(resources[sieve+'_nist']['physical_qubits_baseline'][-1], physical_qubits_baseline)
		resources[sieve+'_nist']['physical_qubits_active'][-1] = max(resources[sieve+'_nist']['physical_qubits_active'][-1], physical_qubits_active)
		resources[sieve+'_nist']['time_baseline'][-1] += I_GaussSieve * time_baseline
		resources[sieve+'_nist']['time_active'][-1] += I_GaussSieve * time_active
		resources[sieve+'_nist']['reaction_limit'][-1] += I_GaussSieve * reaction_limit

		


for sieve in ['NVSieve_plain', 'NVSieve_angular', 'NVSieve_spherical', 'NVSieve_filter', 'GaussSieve_plain', 'GaussSieve_angular', 'GaussSieve_spherical', 'GaussSieve_filter']:
	for condition in ['', '_no_qram', '_nist']:
		with open('data/heuristics_' + sieve + condition + '.txt', "w") as data_file:
			for string in ['physical_qubits_baseline', 'physical_qubits_active', 'time_baseline', 'time_active', 'reaction_limit', 'hashing_time', 'classical_time']:
				data_file.write(string)
				data_file.write('\n')
				for i in range(len(dimensions)):
					data_file.write(str(dimensions[i]))
					data_file.write(' ')
					if string in ['hashing_time', 'classical_time']:
						data_file.write(str(resources[sieve][string][i]))
					else:
						data_file.write(str(resources[sieve+condition][string][i]))
					data_file.write('\n')
