import numpy as np
from mpmath import mp
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import auxiliary_functions as auxiliary
import scienceplots

plt.style.use('science')

dimensions = [20*i + 100 for i in range(46)]
dimensions2 = [20*i + 400 for i in range(6)]
limits = range(40,72)

heuristics = {
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
		'reaction_limit': [],
		'hashing_time': [],
		'classical_time': []
	},
	'NVSieve_angular_no_qram': {
		'physical_qubits_baseline': [],
		'physical_qubits_active': [],
		'time_baseline': [],
		'time_active': [],
		'reaction_limit': [],
		'hashing_time': [],
		'classical_time': []
	},
	'NVSieve_spherical_no_qram': {
		'physical_qubits_baseline': [],
		'physical_qubits_active': [],
		'time_baseline': [],
		'time_active': [],
		'reaction_limit': [],
		'hashing_time': [],
		'classical_time': []
	},
	'NVSieve_filter_no_qram': {
		'physical_qubits_baseline': [],
		'physical_qubits_active': [],
		'time_baseline': [],
		'time_active': [],
		'reaction_limit': [],
		'hashing_time': [],
		'classical_time': []
	},
	'GaussSieve_plain_no_qram': {
		'physical_qubits_baseline': [],
		'physical_qubits_active': [],
		'time_baseline': [],
		'time_active': [],
		'reaction_limit': [],
		'hashing_time': [],
		'classical_time': []
	},
	'GaussSieve_angular_no_qram': {
		'physical_qubits_baseline': [],
		'physical_qubits_active': [],
		'time_baseline': [],
		'time_active': [],
		'reaction_limit': [],
		'hashing_time': [],
		'classical_time': []
	},
	'GaussSieve_spherical_no_qram': {
		'physical_qubits_baseline': [],
		'physical_qubits_active': [],
		'time_baseline': [],
		'time_active': [],
		'reaction_limit': [],
		'hashing_time': [],
		'classical_time': []
	},
	'GaussSieve_filter_no_qram': {
		'physical_qubits_baseline': [],
		'physical_qubits_active': [],
		'time_baseline': [],
		'time_active': [],
		'reaction_limit': [],
		'hashing_time': [],
		'classical_time': []
	},
	'NVSieve_plain_nist': {
		'physical_qubits_baseline': [],
		'physical_qubits_active': [],
		'time_baseline': [],
		'time_active': [],
		'reaction_limit': [],
		'hashing_time': [],
		'classical_time': []
	},
	'NVSieve_angular_nist': {
		'physical_qubits_baseline': [],
		'physical_qubits_active': [],
		'time_baseline': [],
		'time_active': [],
		'reaction_limit': [],
		'hashing_time': [],
		'classical_time': []
	},
	'NVSieve_spherical_nist': {
		'physical_qubits_baseline': [],
		'physical_qubits_active': [],
		'time_baseline': [],
		'time_active': [],
		'reaction_limit': [],
		'hashing_time': [],
		'classical_time': []
	},
	'NVSieve_filter_nist': {
		'physical_qubits_baseline': [],
		'physical_qubits_active': [],
		'time_baseline': [],
		'time_active': [],
		'reaction_limit': [],
		'hashing_time': [],
		'classical_time': []
	},
	'GaussSieve_plain_nist': {
		'physical_qubits_baseline': [],
		'physical_qubits_active': [],
		'time_baseline': [],
		'time_active': [],
		'reaction_limit': [],
		'hashing_time': [],
		'classical_time': []
	},
	'GaussSieve_angular_nist': {
		'physical_qubits_baseline': [],
		'physical_qubits_active': [],
		'time_baseline': [],
		'time_active': [],
		'reaction_limit': [],
		'hashing_time': [],
		'classical_time': []
	},
	'GaussSieve_spherical_nist': {
		'physical_qubits_baseline': [],
		'physical_qubits_active': [],
		'time_baseline': [],
		'time_active': [],
		'reaction_limit': [],
		'hashing_time': [],
		'classical_time': []
	},
	'GaussSieve_filter_nist': {
		'physical_qubits_baseline': [],
		'physical_qubits_active': [],
		'time_baseline': [],
		'time_active': [],
		'reaction_limit': [],
		'hashing_time': [],
		'classical_time': []
	}
}

comparison = {
	'simulations': {
		'physical_qubits_baseline': [],
		'physical_qubits_active': [],
		'time_baseline': [],
		'time_active': [],
		'reaction_limit': []
	},
	'heuristics': {
		'physical_qubits_baseline': [],
		'physical_qubits_active': [],
		'time_baseline': [],
		'time_active': [],
		'reaction_limit': []
	}
}

for sieve in ['NVSieve_plain', 'NVSieve_angular', 'NVSieve_spherical', 'NVSieve_filter', 'GaussSieve_plain', 'GaussSieve_angular', 'GaussSieve_spherical', 'GaussSieve_filter']:
	for condition in ['', '_no_qram', '_nist']:

		with open('data/heuristics/heuristics_'+sieve+condition+'.txt', "r") as data_file:
			for line in data_file:
		    		if 'physical_qubits_baseline' in line:
		    			string = 'physical_qubits_baseline'
		    		elif 'physical_qubits_active' in line:
		    			string = 'physical_qubits_active'
		    		elif 'time_baseline' in line:
		    			string = 'time_baseline'
		    		elif 'time_active' in line:
		    			string = 'time_active'
		    		elif 'reaction_limit' in line:
		    			string = 'reaction_limit'
		    		elif 'hashing_time' in line:
		    			string = 'hashing_time'
		    		elif 'classical_time' in line:
		    			string = 'classical_time'
		    		else:
		    			heuristics[sieve+condition][string].append( float((line.strip('\n').split(" "))[1]) )
	 

#		plt.rcParams.update({'font.size': 12})
		fig = plt.figure()
		ax1 = fig.add_axes([0, 0, 1.5, 1.5])
		plt.grid(color='k', linestyle='--', linewidth=0.3)
		plt.plot(dimensions, heuristics[sieve+condition]['physical_qubits_baseline'], label='Baseline')
		plt.plot(dimensions, heuristics[sieve+condition]['physical_qubits_active'], label='Active-volume')
		plt.yscale("log")
		plt.xlim([100,1000])
		ax1.tick_params(axis='x', labelsize=13)
		ax1.tick_params(axis='y', labelsize=14)
		plt.xlabel('Dimension $D$', fontsize=13)
		plt.ylabel('Physical qubits', fontsize=13)
		plt.legend(frameon=True)
		if 'NVSieve' in sieve:
			plt.title('NVSieve', fontsize=13)
		elif 'GaussSieve' in sieve:
			plt.title('GaussSieve', fontsize=13)
			
		if condition == '_no_qram':
			plt.savefig('plots/heuristics_extra/without_qram/heuristics_'+sieve+'_physical_qubits.eps', format='eps', bbox_inches='tight')
		elif condition == '_nist':
			plt.savefig('plots/heuristics_extra/nist/heuristics_'+sieve+'_physical_qubits.eps', format='eps', bbox_inches='tight')
		else:
			plt.savefig('plots/heuristics_extra/with_qram/heuristics_'+sieve+'_physical_qubits.eps', format='eps', bbox_inches='tight')
			
		plt.close()

#		plt.rcParams.update({'font.size': 12})
		fig = plt.figure()
		ax1 = fig.add_axes([0, 0, 1.5, 1.5])
		plt.grid(color='k', linestyle='--', linewidth=0.3)
		plt.plot(dimensions, heuristics[sieve+condition]['time_baseline'], label='Baseline')
		plt.plot(dimensions, heuristics[sieve+condition]['time_active'], label='Active-volume')
		plt.plot(dimensions, heuristics[sieve+condition]['reaction_limit'], label='Reaction limit')
		plt.plot(dimensions, heuristics[sieve+condition]['hashing_time'], label='Classical hashing')
		plt.yscale("log")
		plt.xlim([100,1000])
		ax1.tick_params(axis='x', labelsize=13)
		ax1.tick_params(axis='y', labelsize=14)
		plt.xlabel('Dimension $D$', fontsize=13)
		plt.ylabel('Time (years)', fontsize=13)
		plt.legend(frameon=True)
		if 'NVSieve' in sieve:
			plt.title('NVSieve', fontsize=13)
		elif 'GaussSieve' in sieve:
			plt.title('GaussSieve', fontsize=13)
			
		if condition == '_no_qram':
			plt.savefig('plots/heuristics_extra/without_qram/heuristics_'+sieve+'_quantum_time.eps', format='eps', bbox_inches='tight')
		elif condition == '_nist':
			plt.savefig('plots/heuristics_extra/nist/heuristics_'+sieve+'_quantum_time.eps', format='eps', bbox_inches='tight')
		else:
			plt.savefig('plots/heuristics_extra/with_qram/heuristics_'+sieve+'_quantum_time.eps', format='eps', bbox_inches='tight')

		plt.close()
		
#		plt.rcParams.update({'font.size': 12})
		fig = plt.figure()
		ax1 = fig.add_axes([0, 0, 1.5, 1.5])
		plt.grid(color='k', linestyle='--', linewidth=0.3)
		plt.plot(dimensions, [x + y for x, y in zip(heuristics[sieve+condition]['time_baseline'], heuristics[sieve+condition]['hashing_time'])], label='Baseline')
		plt.plot(dimensions, [x + y for x, y in zip(heuristics[sieve+condition]['time_active'], heuristics[sieve+condition]['hashing_time'])], label='Active-volume')
		plt.plot(dimensions, [x + y for x, y in zip(heuristics[sieve+condition]['reaction_limit'], heuristics[sieve+condition]['hashing_time'])], label='Reaction limit')
		plt.yscale("log")
		plt.xlim([100,1000])
		ax1.tick_params(axis='x', labelsize=13)
		ax1.tick_params(axis='y', labelsize=14)
		plt.xlabel('Dimension $D$', fontsize=13)
		plt.ylabel('Time (years)', fontsize=13)
		plt.legend(frameon=True)
		if 'NVSieve' in sieve:
			plt.title('NVSieve', fontsize=13)
		elif 'GaussSieve' in sieve:
			plt.title('GaussSieve', fontsize=13)
			
		if condition == '_no_qram':
			plt.savefig('plots/heuristics_extra/without_qram/heuristics_'+sieve+'_total_time.eps', format='eps', bbox_inches='tight')
		elif condition == '_nist':
			plt.savefig('plots/heuristics_extra/nist/heuristics_'+sieve+'_total_time.eps', format='eps', bbox_inches='tight')
		else:
			plt.savefig('plots/heuristics_extra/with_qram/heuristics_'+sieve+'_total_time.eps', format='eps', bbox_inches='tight')

		plt.close()
	

#######################################################################################################

with open('data/simulations/simulations_GaussSieve.txt', "r") as data_file:
	for line in data_file:
    		if 'physical_qubits_baseline' in line:
    			string = 'physical_qubits_baseline'
    		elif 'physical_qubits_active' in line:
    			string = 'physical_qubits_active'
    		elif 'time_baseline' in line:
    			string = 'time_baseline'
    		elif 'time_active' in line:
    			string = 'time_active'
    		elif 'reaction_limit' in line:
    			string = 'reaction_limit'
    		else:
    			comparison['simulations'][string].append( float((line.strip('\n').split(" "))[1]) )


for D in limits:
	L_GaussSieve = mp.ceil(mp.power(2, 0.193*D + 2.325))
	iterations_GaussSieve = mp.ceil(mp.power(2, 0.283*D + 0.335))
	
	comparison['heuristics']['physical_qubits_baseline'].append(0)
	comparison['heuristics']['physical_qubits_active'].append(0)
	comparison['heuristics']['time_baseline'].append(0)
	comparison['heuristics']['time_active'].append(0)
	comparison['heuristics']['reaction_limit'].append(0)

	physical_qubits_baseline, physical_qubits_active, time_baseline, time_active, reaction_limit = auxiliary.resources(D, L_GaussSieve, 1, 'GaussSieve1')
	comparison['heuristics']['physical_qubits_baseline'][-1] = max(comparison['heuristics']['physical_qubits_baseline'][-1], physical_qubits_baseline)
	comparison['heuristics']['physical_qubits_active'][-1] = max(comparison['heuristics']['physical_qubits_active'][-1], physical_qubits_active)
	comparison['heuristics']['time_baseline'][-1] += 9 * iterations_GaussSieve * time_baseline
	comparison['heuristics']['time_active'][-1] += 9 * iterations_GaussSieve * time_active
	comparison['heuristics']['reaction_limit'][-1] += 9 * iterations_GaussSieve * reaction_limit
	
	physical_qubits_baseline, physical_qubits_active, time_baseline, time_active, reaction_limit = auxiliary.resources(D, L_GaussSieve, 0, 'GaussSieve1')
	comparison['heuristics']['physical_qubits_baseline'][-1] = max(comparison['heuristics']['physical_qubits_baseline'][-1], physical_qubits_baseline)
	comparison['heuristics']['physical_qubits_active'][-1] = max(comparison['heuristics']['physical_qubits_active'][-1], physical_qubits_active)
	comparison['heuristics']['time_baseline'][-1] += iterations_GaussSieve * time_baseline
	comparison['heuristics']['time_active'][-1] += iterations_GaussSieve * time_active
	comparison['heuristics']['reaction_limit'][-1] += iterations_GaussSieve * reaction_limit
	
	physical_qubits_baseline, physical_qubits_active, time_baseline, time_active, reaction_limit = auxiliary.resources(D, L_GaussSieve, 0, 'GaussSieve2')
	comparison['heuristics']['physical_qubits_baseline'][-1] = max(comparison['heuristics']['physical_qubits_baseline'][-1], physical_qubits_baseline)
	comparison['heuristics']['physical_qubits_active'][-1] = max(comparison['heuristics']['physical_qubits_active'][-1], physical_qubits_active)
	comparison['heuristics']['time_baseline'][-1] += iterations_GaussSieve * time_baseline
	comparison['heuristics']['time_active'][-1] += iterations_GaussSieve * time_active
	comparison['heuristics']['reaction_limit'][-1] += iterations_GaussSieve * reaction_limit

	
#######################################################################################################
#######################################################################################################

for sieve in ['NVSieve', 'GaussSieve']:
	for condition in ['', '_no_qram', '_nist']:
		for resource in ['physical_qubits_baseline', 'physical_qubits_active']:
		
			if condition == '_no_qram':
				style = 'dashed'
			else:
				style = 'solid'
	
			fig = plt.figure()
			ax1 = fig.add_axes([0, 0, 1.5, 1.5])
			plt.grid(color='k', linestyle='--', linewidth=0.3)
			plt.plot(dimensions, heuristics[sieve+'_plain'+condition][resource], label='No LSH/LSF', color='orange', linestyle=style)
			plt.plot(dimensions, heuristics[sieve+'_angular'+condition][resource], label='Angular LSH', color='blue', linestyle=style)
			plt.plot(dimensions, heuristics[sieve+'_spherical'+condition][resource], label='Spherical LSH', color='green', linestyle=style)
			plt.plot(dimensions, heuristics[sieve+'_filter'+condition][resource], label='Spherical LSF', color='red', linestyle=style)
			plt.yscale("log")
			plt.gca().yaxis.set_ticks_position('both')
			plt.xlim([100,1000])
			ax1.tick_params(axis='x', labelsize=13)
			ax1.tick_params(axis='y', labelsize=14)
			plt.xlabel('Dimension $D$', fontsize=13)
			plt.title(sieve, fontsize=13)
			
			if resource == 'physical_qubits_baseline':
				plt.ylabel('Physical qubits (baseline)', fontsize=13)
				
			elif resource == 'physical_qubits_active':
				plt.ylabel('Physical qubits (active-volume)', fontsize=13)
				
	
			plt.legend(frameon=True)
			if condition == '_no_qram':
				plt.savefig('plots/heuristics/without_qram/'+resource+'_'+sieve+condition+'.eps', format='eps', bbox_inches='tight')
			elif condition == '_nist':
				plt.savefig('plots/heuristics/nist/'+resource+'_'+sieve+condition+'.eps', format='eps', bbox_inches='tight')
			else:
				plt.savefig('plots/heuristics/with_qram/'+resource+'_'+sieve+condition+'.eps', format='eps', bbox_inches='tight')
				
			plt.close()
		
		
		for resource in ['time_baseline', 'time_active', 'reaction_limit']:
		
			if condition == '_no_qram':
				style = 'dashed'
			else:
				style = 'solid'
			
			fig = plt.figure()
			ax1 = fig.add_axes([0, 0, 1.5, 1.5])
			plt.grid(color='k', linestyle='--', linewidth=0.3)
			plt.plot(dimensions, [x + y for x, y in zip(heuristics[sieve+'_plain'+condition][resource], heuristics[sieve+'_plain']['hashing_time'])], label='No LSH/LSF (quantum)', color='orange', linestyle=style)
			plt.plot(dimensions, [x + y for x, y in zip(heuristics[sieve+'_angular'+condition][resource], heuristics[sieve+'_angular']['hashing_time'])], label='Angular LSH (quantum)', color='blue', linestyle=style)
			plt.plot(dimensions, [x + y for x, y in zip(heuristics[sieve+'_spherical'+condition][resource], heuristics[sieve+'_spherical']['hashing_time'])], label='Spherical LSH (quantum)', color='green', linestyle=style)
			plt.plot(dimensions, [x + y for x, y in zip(heuristics[sieve+'_filter'+condition][resource], heuristics[sieve+'_filter']['hashing_time'])], label='Spherical LSF (quantum)', color='red', linestyle=style)
			plt.yscale("log")
			plt.gca().yaxis.set_ticks_position('both')
			plt.xlim([100,1000])
			ax1.tick_params(axis='x', labelsize=13)
			ax1.tick_params(axis='y', labelsize=14)
			plt.xlabel('Dimension $D$', fontsize=13)
			plt.title(sieve, fontsize=13)
			
			if resource == 'time_baseline':
				plt.ylabel('Baseline time (years)', fontsize=13)
			
			elif resource == 'time_active':
				plt.ylabel('Active-volume time (years)', fontsize=13)			
				
			elif resource == 'reaction_limit':
				plt.ylabel('Reaction limit (years)', fontsize=13)
			
			plt.legend(frameon=True)				
			if condition == '_no_qram':
				plt.savefig('plots/heuristics/without_qram/'+resource+'_'+sieve+condition+'.eps', format='eps', bbox_inches='tight')
			elif condition == '_nist':
				plt.savefig('plots/heuristics/nist/'+resource+'_'+sieve+condition+'.eps', format='eps', bbox_inches='tight')
			else:
				plt.savefig('plots/heuristics/with_qram/'+resource+'_'+sieve+condition+'.eps', format='eps', bbox_inches='tight')
			
			plt.close()
		
	for resource in ['time_baseline', 'time_active', 'reaction_limit']:
	
	
		fig = plt.figure()
		ax1 = fig.add_axes([0, 0, 1.5, 1.5])
		
#		plt.rcParams.update({'font.size': 12})
		ax1.grid(color='k', linestyle='--', linewidth=0.3)
		ax1.plot(dimensions, [x + y for x, y in zip(heuristics[sieve+'_plain'][resource], heuristics[sieve+'_plain']['hashing_time'])], label='No LSH/LSF (quantum)', color='orange', linestyle='solid')
		ax1.plot(dimensions, heuristics[sieve+'_plain']['classical_time'], label='No LSH/LSF (classical)', color='orange', linestyle='dashed')
		ax1.plot(dimensions, [x + y for x, y in zip(heuristics[sieve+'_angular'][resource], heuristics[sieve+'_angular']['hashing_time'])], label='Angular LSH (quantum)', color='blue', linestyle='solid')
		ax1.plot(dimensions, heuristics[sieve+'_angular']['classical_time'], label='Angular LSH (classical)', color='blue', linestyle='dashed')
		ax1.plot(dimensions, [x + y for x, y in zip(heuristics[sieve+'_spherical'][resource], heuristics[sieve+'_spherical']['hashing_time'])], label='Spherical LSH (quantum)', color='green', linestyle='solid')
		ax1.plot(dimensions, heuristics[sieve+'_spherical']['classical_time'], label='Spherical LSH (classical)', color='green', linestyle='dashed')
		ax1.plot(dimensions, [x + y for x, y in zip(heuristics[sieve+'_filter'][resource], heuristics[sieve+'_filter']['hashing_time'])], label='Spherical LSF (quantum)', color='red', linestyle='solid')
		ax1.plot(dimensions, heuristics[sieve+'_filter']['classical_time'], label='Spherical LSF (classical)', color='red', linestyle='dashed')
		
		ax2 = fig.add_axes([0.03, 0.8, 0.68, 0.68])
		ax2.grid(color='k', linestyle='--', linewidth=0.5)
		
		ax2.plot(dimensions2, [heuristics[sieve+'_plain'][resource][i] + heuristics[sieve+'_plain']['hashing_time'][i] for i in range(14,20)], label='No LSH/LSF (quantum)', color='orange', linestyle='solid')
		ax2.plot(dimensions2, [heuristics[sieve+'_plain']['classical_time'][i] for i in range(14,20)], label='No LSH/LSF (classical)', color='orange', linestyle='dashed')
		ax2.plot(dimensions2, [heuristics[sieve+'_angular'][resource][i] + heuristics[sieve+'_plain']['hashing_time'][i] for i in range(14,20)], label='Angular LSH (quantum)', color='blue', linestyle='solid')
		ax2.plot(dimensions2, [heuristics[sieve+'_angular']['classical_time'][i] for i in range(14,20)], label='Angular LSH (classical)', color='blue', linestyle='dashed')
		ax2.plot(dimensions2, [heuristics[sieve+'_spherical'][resource][i] + heuristics[sieve+'_plain']['hashing_time'][i] for i in range(14,20)], label='Spherical LSH (quantum)', color='green', linestyle='solid')
		ax2.plot(dimensions2, [heuristics[sieve+'_spherical']['classical_time'][i] for i in range(14,20)], label='Spherical LSH (classical)', color='green', linestyle='dashed')
		ax2.plot(dimensions2, [heuristics[sieve+'_filter'][resource][i] + heuristics[sieve+'_plain']['hashing_time'][i] for i in range(14,20)], label='Spherical LSF (quantum)', color='red', linestyle='solid')
		ax2.plot(dimensions2, [heuristics[sieve+'_filter']['classical_time'][i] for i in range(14,20)], label='Spherical LSF (classical)', color='red', linestyle='dashed')
		
		ax1.set_xlabel('Dimension $D$', fontsize=13)
		ax1.set_yscale("log")

		ax1.tick_params(axis='x', labelsize=13)
		ax1.tick_params(axis='y', labelsize=14)
		ax2.tick_params(axis='x', labelsize=10)
		ax2.tick_params(axis='y', labelsize=11)
		ax2.set_yscale("log")
		ax2.yaxis.tick_right()
#		plt.gca().yaxis.set_ticks_position('both')
		ax1.set_xlim(100,1000)
		ax2.set_xlim(400,500)
		ax2.set_ylim(None, 2*heuristics[sieve+'_angular']['classical_time'][19])
		ax1.set_title(sieve, fontsize=13)
		
		if resource == 'time_baseline':
			ax1.set_ylabel('Baseline time (years)', fontsize=13)
		
		elif resource == 'time_active':
			ax1.set_ylabel('Active-volume time (years)', fontsize=13)			
			
		elif resource == 'reaction_limit':
			ax1.tick_params(labelright=True)
			ax1.set_ylabel('Reaction limit (years)', fontsize=13)
		
		ax1.legend(loc=4, fontsize=9.5, frameon=True)
		plt.savefig('plots/heuristics/with_qram/'+resource+'_'+sieve+'_complete.eps', format='eps', bbox_inches='tight')
		
		plt.close()

#################################################################3

#plt.rcParams.update({'font.size': 12})
fig = plt.figure()
ax1 = fig.add_axes([0, 0, 1.5, 1.5])
plt.grid(color='k', linestyle='--', linewidth=0.3)
plt.scatter(limits, comparison['heuristics']['physical_qubits_baseline'], marker='+', color='blue', label='Baseline (heuristic)')
plt.scatter(limits, comparison['simulations']['physical_qubits_baseline'], marker='+', color='orange', label='Baseline (numerical)')
plt.scatter(limits, comparison['heuristics']['physical_qubits_active'], marker='*', color='blue', label='Active-volume (heuristic)')
plt.scatter(limits, comparison['simulations']['physical_qubits_active'], marker='*', color='orange', label='Active-volume (numerical)')
plt.yscale("log")
plt.xlim([40,71])
plt.ylim([10**7,10**9])
ax1.tick_params(axis='x', labelsize=13)
ax1.tick_params(axis='y', labelsize=14)
plt.gca().yaxis.set_ticks_position('both')
plt.xlabel('Dimension $D$', fontsize=13)
plt.ylabel('Physical qubits', fontsize=13)
plt.legend(frameon=True, fontsize=11)
plt.savefig('plots/comparison/comparison_physical_qubits_active.eps', format='eps', bbox_inches='tight')
plt.close()

#plt.rcParams.update({'font.size': 12})
fig = plt.figure()
ax1 = fig.add_axes([0, 0, 1.5, 1.5])
plt.grid(color='k', linestyle='--', linewidth=0.3)
plt.scatter(limits, [365*comparison['heuristics']['time_baseline'][i] for i in range(32)], marker='+', color='green', label='Baseline (heuristic)')
plt.scatter(limits, [365*comparison['simulations']['time_baseline'][i] for i in range(32)], marker='^', color='green', label='Baseline (numerical)')
plt.scatter(limits, [365*comparison['heuristics']['time_active'][i] for i in range(32)], marker='+', color='blue', label='Active-volume (heuristic)')
plt.scatter(limits, [365*comparison['simulations']['time_active'][i] for i in range(32)], marker='^', color='blue', label='Active-volume (numerical)')
plt.scatter(limits, [365*comparison['heuristics']['reaction_limit'][i] for i in range(32)], marker='+', color='red', label='Reaction limit (heuristic)')
plt.scatter(limits, [365*comparison['simulations']['reaction_limit'][i] for i in range(32)], marker='^', color='red', label='Reaction limit (numerical)')
plt.yscale("log")
plt.xlim([40,71])
#plt.ylim([2*10**7,10**9])
ax1.tick_params(axis='x', labelsize=13)
ax1.tick_params(axis='y', labelsize=14)
plt.gca().yaxis.set_ticks_position('both')
plt.xlabel('Dimension $D$', fontsize=13)
plt.ylabel('Time (days)', fontsize=13)
plt.legend(frameon=True, fontsize=10.5)
plt.savefig('plots/comparison/comparison_time.eps', format='eps', bbox_inches='tight')
plt.close()

#######################################################################################################

#plt.rcParams.update({'font.size': 12})
fig = plt.figure()
ax1 = fig.add_axes([0, 0, 1.5, 1.5])
plt.grid(color='k', linestyle='--', linewidth=0.3)
#plt.plot(dimensions, heuristics['GaussSieve_plain']['physical_qubits_active'], label='No LSH/LSF (quantum)', color='orange', linestyle='solid')
plt.plot(dimensions, heuristics['GaussSieve_angular']['physical_qubits_active'], label='Angular LSH', color='blue', linestyle='solid')
plt.plot(dimensions, heuristics['GaussSieve_spherical']['physical_qubits_active'], label='Spherical LSH ', color='green', linestyle='solid')
plt.plot(dimensions, heuristics['GaussSieve_filter']['physical_qubits_active'], label='Spherical LSF', color='red', linestyle='solid')
plt.yscale("log")
plt.gca().yaxis.set_ticks_position('both')
plt.xlim(200,1000)
#plt.ylim(0.9*10**9, None)
#plt.tick_params(labelright=True)
plt.xlabel('Dimension $D$', fontsize=13)
plt.ylabel('Physical qubits (active-volume)', fontsize=13)
ax1.tick_params(axis='x', labelsize=13)
ax1.tick_params(axis='y', labelsize=14)
plt.legend(frameon=True, fontsize=13)
plt.title('GaussSieve', fontsize=13)
plt.savefig('plots/heuristics/with_qram/GaussSieve_physical_qubits_introduction.eps', format='eps', bbox_inches='tight')
plt.close()


fig = plt.figure()
ax1 = fig.add_axes([0, 0, 1.5, 1.5])
ax1.set_xlabel('Dimension $D$', fontsize=13)
ax1.set_ylabel('Execution time (years)', fontsize=13)
ax1.set_title('GaussSieve', fontsize=13)

ax1.grid(color='k', linestyle='--', linewidth=0.3)
ax1.plot(dimensions, [x + y for x, y in zip(heuristics['GaussSieve_angular']['reaction_limit'], heuristics['GaussSieve_angular']['hashing_time'])], label='Angular LSH (quantum)', color='blue', linestyle='solid')
ax1.plot(dimensions, heuristics['GaussSieve_angular']['classical_time'], label='Angular LSH (classical)', color='blue', linestyle='dashed')
ax1.plot(dimensions, [x + y for x, y in zip(heuristics['GaussSieve_spherical']['reaction_limit'], heuristics['GaussSieve_spherical']['hashing_time'])], label='Spherical LSH (quantum)', color='green', linestyle='solid')
ax1.plot(dimensions, heuristics['GaussSieve_spherical']['classical_time'], label='Spherical LSH (classical)', color='green', linestyle='dashed')
ax1.plot(dimensions, [x + y for x, y in zip(heuristics['GaussSieve_filter']['reaction_limit'], heuristics['GaussSieve_filter']['hashing_time'])], label='Spherical LSF (quantum)', color='red', linestyle='solid')
ax1.plot(dimensions, heuristics['GaussSieve_filter']['classical_time'], label='Spherical LSF (classical)', color='red', linestyle='dashed')

ax2 = fig.add_axes([0.04, 0.8, 0.65, 0.65])
ax2.grid(color='k', linestyle='--', linewidth=0.5)
ax2.plot(dimensions2, [heuristics['GaussSieve_angular']['reaction_limit'][i] + heuristics['GaussSieve_angular']['hashing_time'][i] for i in range(14,20)], label='Angular LSH (quantum)', color='blue', linestyle='solid')
ax2.plot(dimensions2, [heuristics['GaussSieve_angular']['classical_time'][i] for i in range(14,20)], label='Angular LSH (classical)', color='blue', linestyle='dashed')
ax2.plot(dimensions2, [heuristics['GaussSieve_spherical']['reaction_limit'][i] + heuristics['GaussSieve_plain']['hashing_time'][i] for i in range(14,20)], label='Spherical LSH (quantum)', color='green', linestyle='solid')
ax2.plot(dimensions2, [heuristics['GaussSieve_spherical']['classical_time'][i] for i in range(14,20)], label='Spherical LSH (classical)', color='green', linestyle='dashed')
ax2.plot(dimensions2, [heuristics['GaussSieve_filter']['reaction_limit'][i] + heuristics['GaussSieve_plain']['hashing_time'][i] for i in range(14,20)], label='Spherical LSF (quantum)', color='red', linestyle='solid')
ax2.plot(dimensions2, [heuristics['GaussSieve_filter']['classical_time'][i] for i in range(14,20)], label='Spherical LSF (classical)', color='red', linestyle='dashed')
ax1.set_yscale("log")
ax1.tick_params(axis='x', labelsize=13)
ax1.tick_params(axis='y', labelsize=14)
ax2.tick_params(axis='x', labelsize=10)
ax2.tick_params(axis='y', labelsize=11)
ax2.set_yscale("log")
ax2.yaxis.tick_right()
#plt.gca().yaxis.set_ticks_position('both')
ax1.set_xlim(200,1000)
ax2.set_xlim(400,500)
ax1.set_ylim(10**9, None)
ax2.set_ylim(None, 2*heuristics['GaussSieve_angular']['classical_time'][19])
ax1.legend(loc=4, fontsize=12.5, frameon=True)
plt.savefig('plots/heuristics/with_qram/GaussSieve_reaction_limit_introduction.eps', format='eps', bbox_inches='tight')
#plt.show()
plt.close()

