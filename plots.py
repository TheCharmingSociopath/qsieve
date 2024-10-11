import numpy as np
from mpmath import mp
import auxiliary_functions as auxiliary
import matplotlib.pyplot as plt

mp.dps = 30

dimensions = [50*i + 50 for i in range(10)]


NVSieve = {
	'physical_qubits_baseline': [],
	'physical_qubits_active': [],
	'time_baseline': [],
	'time_active': [],
	'reaction_limit': []
}
NVSieve_angular = {
	'physical_qubits_baseline': [],
	'physical_qubits_active': [],
	'time_baseline': [],
	'time_active': [],
	'reaction_limit': []
}
NVSieve_spherical = {
	'physical_qubits_baseline': [],
	'physical_qubits_active': [],
	'time_baseline': [],
	'time_active': [],
	'reaction_limit': []
}
NVSieve_filter = {
	'physical_qubits_baseline': [],
	'physical_qubits_active': [],
	'time_baseline': [],
	'time_active': [],
	'reaction_limit': []
}
GaussSieve = {
	'physical_qubits_baseline': [],
	'physical_qubits_active': [],
	'time_baseline': [],
	'time_active': [],
	'reaction_limit': []
}
GaussSieve_angular = {
	'physical_qubits_baseline': [],
	'physical_qubits_active': [],
	'time_baseline': [],
	'time_active': [],
	'reaction_limit': []
}
GaussSieve_spherical = {
	'physical_qubits_baseline': [],
	'physical_qubits_active': [],
	'time_baseline': [],
	'time_active': [],
	'reaction_limit': []
}
GaussSieve_filter = {
	'physical_qubits_baseline': [],
	'physical_qubits_active': [],
	'time_baseline': [],
	'time_active': [],
	'reaction_limit': []
}


for D in dimensions:

	L_NVSieve = 100*mp.exp(0.163*D + 0.102*mp.log(D) + 1.73)
	S_NVSieve = mp.exp(0.163*D + 0.102*mp.log(D) + 1.73)
	L_GaussSieve = mp.power(2, 0.199*D + 2.149)
	
	iterations_NVSieve = 100					# L_NVSieve / S_NVSieve
	iterations_GaussSieve = 0.9*mp.power(2, 0.283*D + 0.335)


	physical_qubits_baseline, physical_qubits_active, time_baseline, time_active, reaction_limit = auxiliary.resources(D, S_NVSieve, 1, 'NVSieve')
	
	NVSieve['physical_qubits_baseline'].append(physical_qubits_baseline)
	NVSieve['physical_qubits_active'].append(physical_qubits_active)
	NVSieve['time_baseline'].append(iterations_NVSieve * L_NVSieve * time_baseline)
	NVSieve['time_active'].append(iterations_NVSieve * L_NVSieve * time_active)
	NVSieve['reaction_limit'].append(iterations_NVSieve * L_NVSieve * reaction_limit)
	
	
	C_angular = auxiliary.angular_LSH(D, S_NVSieve, 'q')
	physical_qubits_baseline, physical_qubits_active, time_baseline, time_active, reaction_limit = auxiliary.resources(D, C_angular, 1, 'NVSieve')
	
	NVSieve_angular['physical_qubits_baseline'].append(physical_qubits_baseline)
	NVSieve_angular['physical_qubits_active'].append(physical_qubits_active)
	NVSieve_angular['time_baseline'].append(iterations_NVSieve * L_NVSieve * time_baseline)
	NVSieve_angular['time_active'].append(iterations_NVSieve * L_NVSieve * time_active)
	NVSieve_angular['reaction_limit'].append(iterations_NVSieve * L_NVSieve * reaction_limit)
	
	
	C_spherical = auxiliary.spherical_LSH(D, S_NVSieve, 'q')
	physical_qubits_baseline, physical_qubits_active, time_baseline, time_active, reaction_limit = auxiliary.resources(D, C_spherical, 1, 'NVSieve')
	
	NVSieve_spherical['physical_qubits_baseline'].append(physical_qubits_baseline)
	NVSieve_spherical['physical_qubits_active'].append(physical_qubits_active)
	NVSieve_spherical['time_baseline'].append(iterations_NVSieve * L_NVSieve * time_baseline)
	NVSieve_spherical['time_active'].append(iterations_NVSieve * L_NVSieve * time_active)
	NVSieve_spherical['reaction_limit'].append(iterations_NVSieve * L_NVSieve * reaction_limit)
	
	
	C_filter = auxiliary.spherical_LSF(D, S_NVSieve, 'q')
	physical_qubits_baseline, physical_qubits_active, time_baseline, time_active, reaction_limit = auxiliary.resources(D, C_filter, 1, 'NVSieve')
	
	NVSieve_filter['physical_qubits_baseline'].append(physical_qubits_baseline)
	NVSieve_filter['physical_qubits_active'].append(physical_qubits_active)
	NVSieve_filter['time_baseline'].append(iterations_NVSieve * L_NVSieve * time_baseline)
	NVSieve_filter['time_active'].append(iterations_NVSieve * L_NVSieve * time_active)
	NVSieve_filter['reaction_limit'].append(iterations_NVSieve * L_NVSieve * reaction_limit)
	
	
		
	GaussSieve['physical_qubits_baseline'].append(0)
	GaussSieve['physical_qubits_active'].append(0)
	GaussSieve['time_baseline'].append(0)
	GaussSieve['time_active'].append(0)
	GaussSieve['reaction_limit'].append(0)
	
	GaussSieve_angular['physical_qubits_baseline'].append(0)
	GaussSieve_angular['physical_qubits_active'].append(0)
	GaussSieve_angular['time_baseline'].append(0)
	GaussSieve_angular['time_active'].append(0)
	GaussSieve_angular['reaction_limit'].append(0)
	
	GaussSieve_spherical['physical_qubits_baseline'].append(0)
	GaussSieve_spherical['physical_qubits_active'].append(0)
	GaussSieve_spherical['time_baseline'].append(0)
	GaussSieve_spherical['time_active'].append(0)
	GaussSieve_spherical['reaction_limit'].append(0)
	
	GaussSieve_filter['physical_qubits_baseline'].append(0)
	GaussSieve_filter['physical_qubits_active'].append(0)
	GaussSieve_filter['time_baseline'].append(0)
	GaussSieve_filter['time_active'].append(0)
	GaussSieve_filter['reaction_limit'].append(0)
	
	C_angular = auxiliary.angular_LSH(D, L_GaussSieve, 'q')
	C_spherical = auxiliary.spherical_LSH(D, L_GaussSieve, 'q')
	C_filter = auxiliary.spherical_LSF(D, L_GaussSieve, 'q')
	
	for M in range(11):
		for i in [1,2]:
			physical_qubits_baseline, physical_qubits_active, time_baseline, time_active, reaction_limit = auxiliary.resources(D, L_GaussSieve, M, 'GaussSieve' + str(i))
				
			GaussSieve['physical_qubits_baseline'][-1] = max(GaussSieve['physical_qubits_baseline'][-1], physical_qubits_baseline)
			GaussSieve['physical_qubits_active'][-1] = max(GaussSieve['physical_qubits_active'][-1], physical_qubits_active)
			GaussSieve['time_baseline'][-1] += iterations_GaussSieve * time_baseline
			GaussSieve['time_active'][-1] += iterations_GaussSieve * time_active
			GaussSieve['reaction_limit'][-1] += iterations_GaussSieve * reaction_limit
			
			
			physical_qubits_baseline, physical_qubits_active, time_baseline, time_active, reaction_limit = auxiliary.resources(D, C_angular, M, 'GaussSieve' + str(i))
				
			GaussSieve_angular['physical_qubits_baseline'][-1] = max(GaussSieve_angular['physical_qubits_baseline'][-1], physical_qubits_baseline)
			GaussSieve_angular['physical_qubits_active'][-1] = max(GaussSieve_angular['physical_qubits_active'][-1], physical_qubits_active)
			GaussSieve_angular['time_baseline'][-1] += iterations_GaussSieve * time_baseline
			GaussSieve_angular['time_active'][-1] += iterations_GaussSieve * time_active
			GaussSieve_angular['reaction_limit'][-1] += iterations_GaussSieve * reaction_limit
			
			
			physical_qubits_baseline, physical_qubits_active, time_baseline, time_active, reaction_limit = auxiliary.resources(D, C_spherical, M, 'GaussSieve' + str(i))
				
			GaussSieve_spherical['physical_qubits_baseline'][-1] = max(GaussSieve_spherical['physical_qubits_baseline'][-1], physical_qubits_baseline)
			GaussSieve_spherical['physical_qubits_active'][-1] = max(GaussSieve_spherical['physical_qubits_active'][-1], physical_qubits_active)
			GaussSieve_spherical['time_baseline'][-1] += iterations_GaussSieve * time_baseline
			GaussSieve_spherical['time_active'][-1] += iterations_GaussSieve * time_active
			GaussSieve_spherical['reaction_limit'][-1] += iterations_GaussSieve * reaction_limit
			
			
			physical_qubits_baseline, physical_qubits_active, time_baseline, time_active, reaction_limit = auxiliary.resources(D, C_filter, M, 'GaussSieve' + str(i))
				
			GaussSieve_filter['physical_qubits_baseline'][-1] = max(GaussSieve_filter['physical_qubits_baseline'][-1], physical_qubits_baseline)
			GaussSieve_filter['physical_qubits_active'][-1] = max(GaussSieve_filter['physical_qubits_active'][-1], physical_qubits_active)
			GaussSieve_filter['time_baseline'][-1] += iterations_GaussSieve * time_baseline
			GaussSieve_filter['time_active'][-1] += iterations_GaussSieve * time_active
			GaussSieve_filter['reaction_limit'][-1] += iterations_GaussSieve * reaction_limit
		


		


plt.grid(color='k', linestyle='--', linewidth=0.5)
plt.plot(dimensions, NVSieve['physical_qubits_baseline'], label='Baseline')
plt.plot(dimensions, NVSieve['physical_qubits_active'], label='Active-volume')
plt.yscale("log")
plt.title("NVSieve without LSH")
plt.xlabel('Dimension $D$')
plt.ylabel('Physical qubits')
plt.legend()
plt.savefig('plots/NVSieve_plain_physical_qubits.png', bbox_inches='tight')
plt.close()


plt.grid(color='k', linestyle='--', linewidth=0.5)
plt.plot(dimensions, NVSieve['time_baseline'], label='Baseline')
plt.plot(dimensions, NVSieve['time_active'], label='Active-volume')
plt.plot(dimensions, NVSieve['reaction_limit'], label='Reaction limit')
plt.yscale("log")
plt.title("NVSieve without LSH")
plt.xlabel('Dimension $D$')
plt.ylabel('Time (years)')
plt.legend()
plt.savefig('plots/NVSieve_plain_time.png', bbox_inches='tight')
plt.close()


plt.grid(color='k', linestyle='--', linewidth=0.5)
plt.plot(dimensions, NVSieve_angular['physical_qubits_baseline'], label='Baseline')
plt.plot(dimensions, NVSieve_angular['physical_qubits_active'], label='Active-volume')
plt.yscale("log")
plt.title("NVSieve with angular LSH")
plt.xlabel('Dimension $D$')
plt.ylabel('Physical qubits')
plt.legend()
plt.savefig('plots/NVSieve_angular_physical_qubits.png', bbox_inches='tight')
plt.close()


plt.grid(color='k', linestyle='--', linewidth=0.5)
plt.plot(dimensions, NVSieve_angular['time_baseline'], label='Baseline')
plt.plot(dimensions, NVSieve_angular['time_active'], label='Active-volume')
plt.plot(dimensions, NVSieve_angular['reaction_limit'], label='Reaction limit')
plt.yscale("log")
plt.title("NVSieve with angular LSH")
plt.xlabel('Dimension $D$')
plt.ylabel('Time (years)')
plt.legend()
plt.savefig('plots/NVSieve_angular_time.png', bbox_inches='tight')
plt.close()


plt.grid(color='k', linestyle='--', linewidth=0.5)
plt.plot(dimensions, NVSieve_spherical['physical_qubits_baseline'], label='Baseline')
plt.plot(dimensions, NVSieve_spherical['physical_qubits_active'], label='Active-volume')
plt.yscale("log")
plt.title("NVSieve with spherical LSH")
plt.xlabel('Dimension $D$')
plt.ylabel('Physical qubits')
plt.legend()
plt.savefig('plots/NVSieve_spherical_physical_qubits.png', bbox_inches='tight')
plt.close()


plt.grid(color='k', linestyle='--', linewidth=0.5)
plt.plot(dimensions, NVSieve_spherical['time_baseline'], label='Baseline')
plt.plot(dimensions, NVSieve_spherical['time_active'], label='Active-volume')
plt.plot(dimensions, NVSieve_spherical['reaction_limit'], label='Reaction limit')
plt.yscale("log")
plt.title("NVSieve with spherical LSH")
plt.xlabel('Dimension $D$')
plt.ylabel('Time (years)')
plt.legend()
plt.savefig('plots/NVSieve_spherical_time.png', bbox_inches='tight')
plt.close()


plt.grid(color='k', linestyle='--', linewidth=0.5)
plt.plot(dimensions, NVSieve_filter['physical_qubits_baseline'], label='Baseline')
plt.plot(dimensions, NVSieve_filter['physical_qubits_active'], label='Active-volume')
plt.yscale("log")
plt.title("NVSieve with spherical LSF")
plt.xlabel('Dimension $D$')
plt.ylabel('Physical qubits')
plt.legend()
plt.savefig('plots/NVSieve_filter_physical_qubits.png', bbox_inches='tight')
plt.close()


plt.grid(color='k', linestyle='--', linewidth=0.5)
plt.plot(dimensions, NVSieve_filter['time_baseline'], label='Baseline')
plt.plot(dimensions, NVSieve_filter['time_active'], label='Active-volume')
plt.plot(dimensions, NVSieve_filter['reaction_limit'], label='Reaction limit')
plt.yscale("log")
plt.title("NVSieve with spherical LSF")
plt.xlabel('Dimension $D$')
plt.ylabel('Time (years)')
plt.legend()
plt.savefig('plots/NVSieve_filter_time.png', bbox_inches='tight')
plt.close()


plt.grid(color='k', linestyle='--', linewidth=0.5)
plt.plot(dimensions, GaussSieve['physical_qubits_baseline'], label='Baseline')
plt.plot(dimensions, GaussSieve['physical_qubits_active'], label='Active-volume')
plt.yscale("log")
plt.title("GaussSieve without LSH")
plt.xlabel('Dimension $D$')
plt.ylabel('Physical qubits')
plt.legend()
plt.savefig('plots/GaussSieve_plain_physical_qubits.png', bbox_inches='tight')
plt.close()


plt.grid(color='k', linestyle='--', linewidth=0.5)
plt.plot(dimensions, GaussSieve['time_baseline'], label='Baseline')
plt.plot(dimensions, GaussSieve['time_active'], label='Active-volume')
plt.plot(dimensions, GaussSieve['reaction_limit'], label='Reaction limit')
plt.yscale("log")
plt.title("GaussSieve without LSH")
plt.xlabel('Dimension $D$')
plt.ylabel('Time (years)')
plt.legend()
plt.savefig('plots/GaussSieve_plain_time.png', bbox_inches='tight')
plt.close()


plt.grid(color='k', linestyle='--', linewidth=0.5)
plt.plot(dimensions, GaussSieve_angular['physical_qubits_baseline'], label='Baseline')
plt.plot(dimensions, GaussSieve_angular['physical_qubits_active'], label='Active-volume')
plt.yscale("log")
plt.title("GaussSieve with angular LSH")
plt.xlabel('Dimension $D$')
plt.ylabel('Physical qubits')
plt.legend()
plt.savefig('plots/GaussSieve_angular_physical_qubits.png', bbox_inches='tight')
plt.close()


plt.grid(color='k', linestyle='--', linewidth=0.5)
plt.plot(dimensions, GaussSieve_angular['time_baseline'], label='Baseline')
plt.plot(dimensions, GaussSieve_angular['time_active'], label='Active-volume')
plt.plot(dimensions, GaussSieve_angular['reaction_limit'], label='Reaction limit')
plt.yscale("log")
plt.title("GaussSieve with angular LSH")
plt.xlabel('Dimension $D$')
plt.ylabel('Time (years)')
plt.legend()
plt.savefig('plots/GaussSieve_angular_time.png', bbox_inches='tight')
plt.close()


plt.grid(color='k', linestyle='--', linewidth=0.5)
plt.plot(dimensions, GaussSieve_spherical['physical_qubits_baseline'], label='Baseline')
plt.plot(dimensions, GaussSieve_spherical['physical_qubits_active'], label='Active-volume')
plt.yscale("log")
plt.title("GaussSieve with spherical LSH")
plt.xlabel('Dimension $D$')
plt.ylabel('Physical qubits')
plt.legend()
plt.savefig('plots/GaussSieve_spherical_physical_qubits.png', bbox_inches='tight')
plt.close()


plt.grid(color='k', linestyle='--', linewidth=0.5)
plt.plot(dimensions, GaussSieve_spherical['time_baseline'], label='Baseline')
plt.plot(dimensions, GaussSieve_spherical['time_active'], label='Active-volume')
plt.plot(dimensions, GaussSieve_spherical['reaction_limit'], label='Reaction limit')
plt.yscale("log")
plt.title("GaussSieve with spherical LSH")
plt.xlabel('Dimension $D$')
plt.ylabel('Time (years)')
plt.legend()
plt.savefig('plots/GaussSieve_spherical_time.png', bbox_inches='tight')
plt.close()


plt.grid(color='k', linestyle='--', linewidth=0.5)
plt.plot(dimensions, GaussSieve_filter['physical_qubits_baseline'], label='Baseline')
plt.plot(dimensions, GaussSieve_filter['physical_qubits_active'], label='Active-volume')
plt.yscale("log")
plt.title("GaussSieve with spherical LSF")
plt.xlabel('Dimension $D$')
plt.ylabel('Physical qubits')
plt.legend()
plt.savefig('plots/GaussSieve_filter_physical_qubits.png', bbox_inches='tight')
plt.close()


plt.grid(color='k', linestyle='--', linewidth=0.5)
plt.plot(dimensions, GaussSieve_filter['time_baseline'], label='Baseline')
plt.plot(dimensions, GaussSieve_filter['time_active'], label='Active-volume')
plt.plot(dimensions, GaussSieve_filter['reaction_limit'], label='Reaction limit')
plt.yscale("log")
plt.title("GaussSieve with spherical LSF")
plt.xlabel('Dimension $D$')
plt.ylabel('Time (years)')
plt.legend()
plt.savefig('plots/GaussSieve_filter_time.png', bbox_inches='tight')
plt.close()


