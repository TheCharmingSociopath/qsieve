from mpmath import mp

mp.dps = 30

K = 32					# Fixed-point precision
C_CCZ = 65
p_phy = mp.power(10, -5)
error = 0.001
epsilon = 0.001
delta = 0.001
reaction_time = mp.power(10, -6)
code_cycle = mp.power(10, -7)


def angular_LSH(D, L, case = 'q'):

	mp.dps = D/10 + 20

	normalisation = mp.quad(lambda theta: mp.power(mp.sin(theta), D-2), [mp.pi/3, mp.pi/2])

	k_angular = lambda t: mp.log(t)/mp.log(3/2) - mp.log(mp.log(1/epsilon))/mp.log(3/2)
	p2_angular = lambda t: mp.quad(lambda theta: mp.power(mp.sin(theta), D-2)*(1 - mp.power(1 - mp.power(1 - theta/mp.pi, k_angular(t) ), t) ), [mp.pi/3, mp.pi/2])/normalisation

	if case == 'q':
		t_angular = mp.findroot(lambda t: L*p2_angular(t) - t*t, 10*D*mp.power(2, 0.078430*D), verify=False)

	else:
		t_angular = mp.findroot(lambda t: L*p2_angular(t) - t,   10*D*mp.power(2, 0.129043*D), verify=False)
		
	print('D:', D)
	print('Angular Error p2:', L * p2_angular(t_angular) - t_angular*t_angular)
	print('')
		
	return L * p2_angular(t_angular)

########################################################

def spherical_LSH(D, L, case = 'q'):

	mp.dps = D/10 + 20

	normalisation = mp.quad(lambda theta: mp.power(mp.sin(theta), D-2), [mp.pi/3, mp.pi/2])

	k_spherical = lambda t: 6*(mp.log(t) - mp.log(mp.log(1/epsilon)))
	p2_spherical = lambda t: mp.quad(lambda theta: mp.power(mp.sin(theta),D-2)*(1 - mp.power(1 - mp.exp(-k_spherical(t)/2*mp.tan(theta/2)*mp.tan(theta/2)), t)), [mp.pi/3, mp.pi/2])/normalisation 

	if case == 'q':
		t_spherical = mp.findroot(lambda t: L*p2_spherical(t) - t*t, 10*D*mp.power(2, 0.059581*D), verify=False)

	else:
		t_spherical = mp.findroot(lambda t: L*p2_spherical(t) - t,   10*D*mp.power(2, 0.089624*D), verify=False)
		
	print('D:', D)
	print('Spherical Error p2:', L * p2_spherical(t_spherical) - t_spherical * t_spherical)
	print('')
	return L * p2_spherical(t_spherical)

########################################################

def spherical_LSF(D, L, case = 'q'):

	mp.dps = D/20 + 30

	Cap = lambda alpha: mp.re( 1/mp.sqrt(mp.pi)*mp.gamma(D/2)/mp.gamma(D/2-1/2) * mp.quad(lambda phi: mp.power(mp.sin(phi), D-2), [0, alpha]) )
	Cap_diff = lambda alpha: mp.re( 1/mp.sqrt(mp.pi)*mp.gamma(D/2)/mp.gamma(D/2-1/2) * mp.power(mp.sin(alpha), D-2) )
	Wedge = lambda alpha: mp.re( 2/mp.pi*(D/2-1) * mp.quad(lambda phi1 : mp.power(mp.sin(phi1), D-2) * mp.quad(lambda phi2: mp.power(mp.sin(phi2), D-3), [0, mp.acos(1/mp.sqrt(3)/mp.tan(phi1))]), [mp.pi/6, alpha]) )
	Wedge_diff = lambda alpha: mp.re( 2/mp.pi*(D/2-1) * mp.power(mp.sin(alpha), D-2) * mp.quad(lambda phi: mp.power(mp.sin(phi), D-3), [0, mp.acos(1/mp.sqrt(3)/mp.tan(alpha))]) )
	t_filter = lambda alpha: mp.log(1/epsilon) / Wedge(alpha)


	def time_c_diff(alpha):
		Cap_value = Cap(alpha)
		Cap_diff_value = Cap_diff(alpha)
		Wedge_value = Wedge(alpha)
		Wedge_diff_value = Wedge_diff(alpha)
		return Cap_diff_value / Wedge_value * (2 + 2 * L * Cap_value) - Cap_value * Wedge_diff_value / mp.power(Wedge_value, 2) * (2 + L * Cap_value)

	def time_q_diff(alpha):
		Cap_value = Cap(alpha)
		Cap_diff_value = Cap_diff(alpha)
		Wedge_value = Wedge(alpha)
		Wedge_diff_value = Wedge_diff(alpha)
		return 2 * mp.sqrt(mp.log(1/epsilon)) * Cap_diff_value / Wedge_value - 2 * mp.sqrt(mp.log(1/epsilon)) * Cap_value * Wedge_diff_value / mp.power(Wedge_value, 2) + mp.sqrt(L) * Cap_diff_value / mp.sqrt(Wedge_value) - mp.sqrt(L)/2 * Cap_value * Wedge_diff_value / mp.power(Wedge_value, 3/2)


	if case == 'q':
		alpha = mp.findroot(time_q_diff, mp.pi/3, verify=False)

	else:
		alpha = mp.findroot(time_c_diff, mp.pi/3, verify=False)
		
	print('D:', D)
	print('Filter Time Error:', time_q_diff(alpha))
	print('')
	
	if alpha < mp.pi/3:
		alpha = mp.pi/3
	
	return L * t_filter(alpha) * mp.power(Cap(alpha), 2)
	
########################################################	

def resources(D, C, M, case = 'NVSieve'):

	if M == 0:
		grover_iterations = mp.ceil(9.2 * mp.log(1/delta, 3) * mp.sqrt(C))
	else:
		grover_iterations = mp.ceil(7.67 * mp.sqrt(C/M))


	if case == 'NVSieve':
		toffoli_count = grover_iterations*(C - 2 + 2*D*(K-1) + D*(K*K - K + 1) + mp.ceil(mp.log(C,2)) - 1)
		logical_qubits = 2*(2*C + D*K - 1 + D*(3*K + 2*K*K) + K)
		reaction_depth = 2*grover_iterations*(mp.ceil(mp.log(C,2)) - 1 + K*mp.log(K,2) - K - mp.log(K,2) + 2 + (mp.ceil(mp.log(D,2)) + 2)*(K-1) + mp.ceil(mp.log(mp.ceil(mp.log(C,2)), 2)) )
		active_volume = grover_iterations*((25 + 1.5*K + C_CCZ)*C + 2*D*((K-1)*(39 + C_CCZ) + 7) + D*(28*K*K - 44*K + 30 + (K*K - K + 1)*C_CCZ) + (mp.ceil(mp.log(C,2)) - 1)*(18 + C_CCZ))
	elif case == 'GaussSieve1':
		toffoli_count = grover_iterations*(C - 2 + (4*D-2)*(K-1) + 2*D*(K*K - K + 1) + mp.ceil(mp.log(C,2)) - 1)
		logical_qubits = 2*(2*C + D*K - 1 + D*(3*K + 2*K*K) + K)
		reaction_depth = grover_iterations*(2*mp.ceil(mp.log(C,2)) - 2 + 2*K*mp.log(K,2) - 2*K - 2*mp.log(K,2) + 4 + 2*(mp.ceil(mp.log(D,2)) + 1)*(K-1) + 2*mp.ceil(mp.log(mp.ceil(mp.log(C,2)), 2)) )
		active_volume = grover_iterations*((25 + 1.5*K + C_CCZ)*C + 4*(2*K*D + 4) + (4*D-2)*((K-1)*(39 + C_CCZ) + 7) + 2*D*(28*K*K - 44*K + 30 + (K*K - K + 1)*C_CCZ) + (mp.ceil(mp.log(C,2)) - 1)*(18 + C_CCZ))
	elif case == 'GaussSieve2':
		toffoli_count = grover_iterations*(C - 2 + (D+1)*(K-1) + D*(0.5*K*K - 1.5*K + 1) + mp.ceil(mp.log(C,2)) - 1)
		logical_qubits = 2*(2*C + D*K - 1 + D*(3*K + 2*K*K) + K)
		reaction_depth = grover_iterations*(2*mp.ceil(mp.log(C,2)) - 2 + 2*K*mp.log(K,2) - 2*K - 2*mp.log(K,2) + 3 + 2*(mp.ceil(mp.log(D,2)) + 1)*(K-1) + 2*mp.ceil(mp.log(mp.ceil(mp.log(C,2)), 2)) )
		active_volume = grover_iterations*((25 + 1.5*K + C_CCZ)*C + 4*(2*K*D + 4) + (D-1)*((K-1)*(39 + C_CCZ) + 7) + D*(20.25*K*K - 48.75*K + 32 + (0.5*K*K - 1.5*K + 1)*C_CCZ) + (mp.ceil(mp.log(C,2)) - 1)*(18 + C_CCZ))
		


	########################################################
	# Code distance


	baseline_spacetime_volume = 4 * (reaction_depth/2) * logical_qubits


	distance_baseline = 1
	error_probabiltiy_baseline = baseline_spacetime_volume * distance_baseline * 0.1 * mp.power(100*p_phy, (distance_baseline + 1)/2) 


	while error_probabiltiy_baseline > error:
		distance_baseline += 1
		error_probabiltiy_baseline = baseline_spacetime_volume * distance_baseline * 0.1 * mp.power(100*p_phy, (distance_baseline + 1)/2)
		
		
	distance_active = 1
	error_probability_active = 2 * active_volume * distance_active * 0.1 * mp.power(100*p_phy, (distance_active + 1)/2) 


	while error_probability_active > error:
		distance_active += 1
		error_probability_active = 2 * active_volume * distance_active * 0.1 * mp.power(100*p_phy, (distance_active + 1)/2)


	if distance_active % 2 == 1:
		distance_active += 1 	 
		error_probability_active = 2 * active_volume * distance_active * 0.1 * mp.power(100*p_phy, (distance_active + 1)/2)


	#######################################################
	# Magic state distillation factories


	dx =  distance_baseline/4
	dz =  distance_baseline/8
	dm =  distance_baseline/8
	dx2 = distance_baseline/2
	dz2 = distance_baseline/4
	dm2 = distance_baseline/4
	dx3 = distance_baseline
	dz3 = distance_baseline/2
	dm3 = distance_baseline/2
	nl1 = 4
	nl2 = 4
	distillation_factories = 4 * max(6 * dm2 / (nl2 / 2), dm3) / (4 * distance_baseline) * C / 2
	qubits_factories = distillation_factories * (2 * mp.ceil((3 * dx3 + dz3) * 3 * dx3 + nl2 * (dx2 + 4 * dz2) * dm3 / 2 + 20 * dm3 * dm3 + 2 * dx3 * dm3) + 2 * nl2 * mp.ceil((3 * dx2 + dz2) * 3 * dx2 + nl1 * ((dx + 4 * dz) * (3 * dx + dm2 / 2) + 2 * dm) + 20 * dm2 * dm2 + 2 * dx2 * dm2) )


	#factory = cost_of_three_level_8toccz(p_phy, dx, dz, dm, dx2, dz2, dm2, dx3, dz3, dm3, 4, 4)
	#distillation_factories[i] += [factory.distillation_time_in_cycles / (4 * distance_baseline[i][j]) * S / 2]
	#qubits_factories[i] += [factory.qubits * distillation_factories[i][j]]		
	#print('Toffoli error:', '{:.2e}'.format(error/toffoli_count[i][j]) )
	#print('Distillation error:', '{:.2e}'.format(factory.distilled_magic_state_error_rate) )


	#######################################################


	physical_qubits_baseline = 2 * distance_baseline * distance_baseline * logical_qubits + qubits_factories
	physical_qubits_active = distance_active * distance_active * logical_qubits
	time_baseline = 4 * distance_baseline * (reaction_depth/2) * code_cycle / 3600/24/365
	time_active = active_volume / logical_qubits / 2 * distance_active * code_cycle / 3600/24/365
	reaction_limit =  reaction_depth * reaction_time / 3600/24/365
	
	return physical_qubits_baseline, physical_qubits_active, time_baseline, time_active, reaction_limit

