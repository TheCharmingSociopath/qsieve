import numpy as np
from mpmath import mp
import os
import auxiliary_functions as aux

mp.dps = 30
dimensions = [10*i + 100 for i in range(41)]


#for sieve in ['NVSieve_spherical', 'GaussSieve_spherical']:
#for sieve in ['NVSieve_angular', 'NVSieve_spherical']:
#for sieve in ['NVSieve_filter']:
for sieve in ['GaussSieve_angular', 'GaussSieve_spherical']:
#for sieve in ['GaussSieve_filter']:
	with open('parameters_hashing_classical.txt', "a") as data_file:
		data_file.write(sieve)
		data_file.write('\n')
		
		for D in dimensions:
			data_file.write(str(D))
			data_file.write(' ')
			
			if sieve == 'NVSieve_angular':
				mp.dps = D/10 + 35
				L = D * mp.ceil(mp.exp(0.163*D + 0.102*mp.log(D) + 1.73))
				t_angular = mp.findroot(lambda t: D * L * aux.p2_angular(t, D) - 3 * aux.k_angular(t) * t, D*mp.power(2, 0.129043*D), verify=False)
				print('D:', D)
				print('Angular Error p2:', D * L * aux.p2_angular(t_angular, D) - 3 * aux.k_angular(t_angular) * t_angular)
				data_file.write(str( t_angular ))
				
			elif sieve == 'NVSieve_spherical':
				mp.dps = D/10 + 30
				L = D * mp.ceil(mp.exp(0.163*D + 0.102*mp.log(D) + 1.73))
				t_spherical = mp.findroot(lambda t: 3 * L * aux.p2_spherical(t, D) - 5 * mp.power(2, mp.sqrt(D)) * aux.k_spherical(t, D) * t, D*mp.power(2, 0.089624*D), verify=False)
				print('D:', D)
				print('Spherical Error p2:', 3 * L * aux.p2_spherical(t_spherical, D) - 5 * mp.power(2, mp.sqrt(D)) * aux.k_spherical(t_spherical, D) * t_spherical)
				data_file.write(str( t_spherical ))
				
			elif sieve == 'NVSieve_filter':
				mp.dps = D/15 + 40
				L = D * mp.ceil(mp.exp(0.163*D + 0.102*mp.log(D) + 1.73))
				alpha = mp.findroot(lambda alpha: aux.time_c_diff_nvsieve(alpha, D, L), mp.pi/3, verify=False)
				print('D:', D)
				print('alpha:', alpha)
				print('Filter Time Error:', aux.time_c_diff_nvsieve(alpha, D, L))
				if alpha < mp.pi/3:
					alpha = mp.pi/3
				data_file.write(str( alpha ))
				
			elif sieve == 'GaussSieve_angular':
				mp.dps = D/6 + 20
				L = mp.ceil(mp.power(2, 0.193*D + 2.325))
				I = mp.ceil(mp.power(2, 0.283*D + 0.335))
				t_angular = mp.findroot(lambda t: (125*D - 19) * I * aux.p2_angular(t, D) - 9 * aux.k_angular(t) * t, D*mp.power(2, 0.129043*D), verify=False)
				print('D:', D)
				print('Angular Error p2:', (125*D - 19) * I * aux.p2_angular(t_angular, D) - 9 * aux.k_angular(t_angular) * t_angular)
				data_file.write(str( t_angular ))
				
			elif sieve == 'GaussSieve_spherical':
				mp.dps = D/6 + 20
				L = mp.ceil(mp.power(2, 0.193*D + 2.325))
				I = mp.ceil(mp.power(2, 0.283*D + 0.335))
				t_spherical = mp.findroot(lambda t: (125 - 19/D) * I * p2_spherical(t, D) - 5 * mp.power(2, mp.sqrt(D)) * aux.k_spherical(t) * t, D*mp.power(2, 0.089624*D), verify=False)
				print('D:', D)
				print('Spherical Error p2:', (125 - 19/D) * I * aux.p2_spherical(t_spherical, D) - 5 * mp.power(2, mp.sqrt(D)) * aux.k_spherical(t_spherical) * t_spherical)
				data_file.write(str( t_spherical ))
				
			elif sieve == 'GaussSieve_filter':
				mp.dps = D/15 + 40
				L = mp.ceil(mp.power(2, 0.193*D + 2.325))
				I = mp.ceil(mp.power(2, 0.283*D + 0.335))
				alpha = mp.findroot(lambda alpha: aux.time_c_diff_gausssieve(alpha, D, L, I), mp.pi/3, verify=False)
				print('D:', D)
				print('alpha:', alpha)
				print('Filter Time Error:', aux.time_c_diff_gausssieve(alpha, D, L, I))
				if alpha < mp.pi/3:
					alpha = mp.pi/3
				data_file.write(str( alpha )) 
			
			data_file.write('\n')
