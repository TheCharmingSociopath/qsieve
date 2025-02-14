import numpy as np
from mpmath import mp
import os
import auxiliary_functions as aux

mp.dps = 30
dimensions = [20*i + 100 for i in range(46)]


#for sieve in ['NVSieve_spherical', 'GaussSieve_spherical']:
#for sieve in ['NVSieve_angular', 'NVSieve_spherical', 'GaussSieve_angular', 'GaussSieve_spherical']:
#for sieve in ['NVSieve_filter']:
#for sieve in ['GaussSieve_angular', 'GaussSieve_spherical']:
for sieve in ['GaussSieve_filter']:
	with open('parameters_hashing_quantum.txt', "a") as data_file:
		data_file.write(sieve)
		data_file.write('\n')
		
		for D in dimensions:
			data_file.write(str(D))
			data_file.write(' ')
			
			if sieve == 'NVSieve_angular':
				mp.dps = D/10 + 20
				L = mp.ceil(mp.exp(0.163*D + 0.102*mp.log(D) + 1.73))
				t_angular = mp.findroot(lambda t: L * aux.p2_angular(t, D) - mp.power(aux.k_angular(t) * t / D/D, 2), D*mp.power(2, 0.078430*D), verify=False)
				print('D:', D)
				print('Angular Error p2:', L * aux.p2_angular(t_angular, D) - mp.power(aux.k_angular(t_angular) * t_angular / D/D, 2))
				data_file.write(str( t_angular ))
				
			elif sieve == 'NVSieve_spherical':
				mp.dps = D/10 + 20
				L = mp.ceil(mp.exp(0.163*D + 0.102*mp.log(D) + 1.73))
				t_spherical = mp.findroot(lambda t: L * aux.p2_spherical(t, D) - mp.power(mp.power(2, mp.sqrt(D)) * aux.k_spherical(t, D) * t / D, 2), 2*mp.power(2, 0.054*D), verify=False)
				print('D:', D)
				print('Spherical Error p2:', L * aux.p2_spherical(t_spherical, D) - mp.power(mp.power(2, mp.sqrt(D)) * aux.k_spherical(t_spherical, D) * t_spherical / D, 2) )
				data_file.write(str( t_spherical ))
				
			elif sieve == 'NVSieve_filter':
				mp.dps = D/20 + 35
				L = mp.ceil(mp.exp(0.163*D + 0.102*mp.log(D) + 1.73))
				alpha = mp.findroot(lambda alpha: aux.time_q_diff_nvsieve(alpha, D, L), mp.pi/3, verify=False)
				print('D:', D)
				print('alpha:', alpha)
				print('Filter Time Error:', aux.time_q_diff_nvsieve(alpha, D, L))
				if alpha < mp.pi/3:
					alpha = mp.pi/3
				data_file.write(str( alpha ))
				
			elif sieve == 'GaussSieve_angular':
				mp.dps = D/5 + 25
				L = mp.ceil(mp.power(2, 0.193*D + 2.325))
				I = mp.ceil(mp.power(2, 0.283*D + 0.335))
				t_angular = mp.findroot(lambda t: I * I / L * aux.p2_angular(t, D) - mp.power(aux.k_angular(t) * t / D, 2), 100*D*D*mp.power(2, 0.1*D), verify=False)
				print('D:', D)
				print('Angular Error p2:', I * I / L * aux.p2_angular(t_angular, D) - mp.power(aux.k_angular(t_angular) * t_angular / D, 2))
				data_file.write(str( t_angular ))
				
			elif sieve == 'GaussSieve_spherical':
				mp.dps = D/5 + 15
				L = mp.ceil(mp.power(2, 0.193*D + 2.325))
				I = mp.ceil(mp.power(2, 0.283*D + 0.335))
				t_spherical = mp.findroot(lambda t: I * I / L * aux.p2_spherical(t, D) - mp.power(mp.power(2, mp.sqrt(D)) * aux.k_spherical(t, D) * t, 2), 10*D*D*mp.power(2, 0.06*D), verify=False)
				print('D:', D)
				print('Spherical Error p2:', I * I / L * aux.p2_spherical(t_spherical, D) - mp.power(mp.power(2, mp.sqrt(D)) * aux.k_spherical(t_spherical, D) * t_spherical, 2) )
				data_file.write(str( t_spherical ))
				
			elif sieve == 'GaussSieve_filter':
				mp.dps = D/15 + 40
				L = mp.ceil(mp.power(2, 0.193*D + 2.325))
				I = mp.ceil(mp.power(2, 0.283*D + 0.335))
				alpha = mp.findroot(lambda alpha: aux.time_q_diff_gausssieve(alpha, D, L, I), mp.pi/3, verify=False)
				print('D:', D)
				print('alpha:', alpha)
				print('Filter Time Error:', aux.time_q_diff_gausssieve(alpha, D, L, I))
				if alpha < mp.pi/3:
					alpha = mp.pi/3
				data_file.write(str( alpha ))
			
			data_file.write('\n')
