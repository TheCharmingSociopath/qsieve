from mpmath import mp

K = 32					# Fixed-point precision
C_CCZ = 65
p_phy = mp.power(10, -5)
error = 0.001
epsilon = 0.001
delta = 0.001
reaction_time = mp.power(10, -6)
code_cycle = mp.power(10, -7)

# Angular LSH functions
k_angular = lambda t: mp.log(t,3/2) - mp.log(mp.log(1/epsilon),3/2)

p2_angular = lambda t, D: mp.quad(lambda theta: mp.power(mp.sin(theta), D-2)*(1 - mp.power(1 - mp.power(1 - theta/mp.pi, k_angular(t) ), t) ), [mp.pi/3, mp.pi/2]) / mp.quad(lambda theta: mp.power(mp.sin(theta), D-2), [mp.pi/3, mp.pi/2])

# Spherical LSH functions
k_spherical = lambda t, D: 6*(mp.log(t) - mp.log(mp.log(1/epsilon)))/mp.sqrt(D)

p2_spherical = lambda t, D: mp.quad(lambda theta: mp.power(mp.sin(theta),D-2)*(1 - mp.power(1 - mp.exp(-k_spherical(t, D)*mp.sqrt(D)/2*mp.tan(theta/2)*mp.tan(theta/2)), t)), [mp.pi/3, mp.pi/2]) / mp.quad(lambda theta: mp.power(mp.sin(theta), D-2), [mp.pi/3, mp.pi/2])

# Spherical LSF functions
Cap = lambda alpha, D: mp.re( 1/mp.sqrt(mp.pi)*mp.gamma(D/2)/mp.gamma(D/2-1/2) * mp.quad(lambda phi: mp.power(mp.sin(phi), D-2), [0, alpha]) )

Cap_diff = lambda alpha, D: mp.re( 1/mp.sqrt(mp.pi)*mp.gamma(D/2)/mp.gamma(D/2-1/2) * mp.power(mp.sin(alpha), D-2) )

Wedge = lambda alpha, D: mp.re( 2/mp.pi*(D/2-1) * mp.quad(lambda phi1 : mp.power(mp.sin(phi1), D-2) * mp.quad(lambda phi2: mp.power(mp.sin(phi2), D-3), [0, mp.acos(1/mp.sqrt(3)/mp.tan(phi1))]), [mp.pi/6, alpha]) )

Wedge_diff = lambda alpha, D: mp.re( 2/mp.pi*(D/2-1) * mp.power(mp.sin(alpha), D-2) * mp.quad(lambda phi: mp.power(mp.sin(phi), D-3), [0, mp.acos(1/mp.sqrt(3)/mp.tan(alpha))]) )

t_filter = lambda alpha, D: mp.log(1/epsilon) / Wedge(alpha, D)

def time_c_diff_nvsieve(alpha, D, L):
	Cap_value = Cap(alpha, D)
	Cap_diff_value = Cap_diff(alpha, D)
	Wedge_value = Wedge(alpha, D)
	Wedge_diff_value = Wedge_diff(alpha, D)
	return (2 * mp.log(D,2) + 6 * D * L * Cap_value) * Cap_diff_value / Wedge_value - (2 * mp.log(D,2) + 3 * D * L * Cap_value) * Cap_value * Wedge_diff_value / mp.power(Wedge_value, 2)
	
def time_c_diff_gausssieve(alpha, D, L, I):
	Cap_value = Cap(alpha, D)
	Cap_diff_value = Cap_diff(alpha, D)
	Wedge_value = Wedge(alpha, D)
	Wedge_diff_value = Wedge_diff(alpha, D)
	return (2 * mp.log(D,2) + 2*(125*D - 19) * I * Cap_value) * Cap_diff_value/Wedge_value - (2*mp.log(D,2) + (125*D - 19) * I * Cap_value) * Cap_value * Wedge_diff_value/mp.power(Wedge_value, 2)

def time_q_diff_nvsieve(alpha, D, L):
	Cap_value = Cap(alpha, D)
	Cap_diff_value = Cap_diff(alpha, D)
	Wedge_value = Wedge(alpha, D)
	Wedge_diff_value = Wedge_diff(alpha, D)
	return mp.log(D,2) * mp.sqrt(mp.log(1/epsilon)) * Cap_diff_value / Wedge_value - mp.log(D,2) * mp.sqrt(mp.log(1/epsilon)) * Cap_value * Wedge_diff_value / mp.power(Wedge_value, 2) + D * D * mp.sqrt(L) * Cap_diff_value / mp.sqrt(Wedge_value) - D * D * mp.sqrt(L)/2 * Cap_value * Wedge_diff_value / mp.power(Wedge_value, 3/2)

def time_q_diff_gausssieve(alpha, D, L, I):
	Cap_value = Cap(alpha, D)
	Cap_diff_value = Cap_diff(alpha, D)
	Wedge_value = Wedge(alpha, D)
	Wedge_diff_value = Wedge_diff(alpha, D)
	return mp.log(D,2) * mp.sqrt(mp.log(1/epsilon)) * Cap_diff_value / Wedge_value - mp.log(D,2) * mp.sqrt(mp.log(1/epsilon)) * Cap_value * Wedge_diff_value / mp.power(Wedge_value, 2) + D * I / L * mp.sqrt(L) * Cap_diff_value / mp.sqrt(Wedge_value) - D * I / L * mp.sqrt(L)/2 * Cap_value * Wedge_diff_value / mp.power(Wedge_value, 3/2)


# Classical hashing cost functions
hashing_angular = lambda D, t, L: 9 * D * k_angular(t) * t * L

hashing_spherical = lambda D, t, L: 5 * D * mp.ceil(mp.power(2, mp.sqrt(D))) * k_spherical(t) * t * L

hashing_filter = lambda D, alpha, L: 2 * mp.ceil(mp.log(D,2)) * Cap(alpha, D) * L * mp.ceil( mp.log(1/epsilon) / Wedge(alpha, D) )

# Classical searching cost functions
searching_classical_NVSieve_plain = lambda D, L: 3 * D * L * L

searching_classical_NVSieve_angular = lambda D, t, L: 3 * D * L * L * p2_angular(t, D)

searching_classical_NVSieve_spherical = lambda D, t, L: 3 * D * L * L * p2_spherical(t, D)

searching_classical_NVSieve_filter = lambda D, alpha, L: 3 * D * mp.power(L * Cap(alpha, D), 2) * mp.ceil( mp.log(1/epsilon) / Wedge(alpha, D) )

searching_classical_GaussSieve_plain = lambda D, L, I: (125 * D - 19) * L * I

searching_classical_GaussSieve_angular = lambda D, t, L, I: (125 * D - 19) * L * I * p2_angular(t, D)

searching_classical_GaussSieve_spherical = lambda D, t, L, I: (125 * D - 19) * L * I * p2_spherical(t, D)

searching_classical_GaussSieve_filter = lambda D, alpha, L, I: (125 * D - 19) * L * I * mp.power(Cap(alpha, D), 2) * mp.ceil( mp.log(1/epsilon) / Wedge(alpha, D) )

	
########################################################	

def resources(D, C, M, case = 'NVSieve'):
	
	# Outputs number of physical qubits and circuit time required by a baseline and an active-volume architectures given dimension D, list size C, and number of solutions M for Grover's search

	if M == 0:
		grover_iterations = mp.ceil(9.2 * mp.log(1/delta, 3) * mp.sqrt(C))
	else:
		grover_iterations = mp.ceil(3.1 * mp.sqrt(C/M))


	if case == 'NVSieve':
		toffoli_count = grover_iterations*(C - 2 + 2*D*(K-1) + D*(K*K - K + 1) + mp.ceil(mp.log(C,2)) - 1)
		logical_qubits = 2*(2*C + D*K - 1 + D*(3*K + 2*K*K) + K)
		reaction_depth = 2*grover_iterations*(mp.ceil(mp.log(C,2)) - 1 + K*mp.log(K,2) - K - mp.log(K,2) + 2 + (mp.ceil(mp.log(D,2)) + 2)*(K-1) + mp.ceil(mp.log(mp.ceil(mp.log(C,2)), 2)) )
		active_volume = grover_iterations*((25 + 1.5*K + C_CCZ)*C + 2*D*((K-1)*(39 + C_CCZ) + 7) + D*(28*K*K - 42*K + 28 + (K*K - K + 1)*C_CCZ) + (mp.ceil(mp.log(C,2)) - 1)*(18 + C_CCZ))
	elif case == 'GaussSieve1':
		toffoli_count = grover_iterations*(C - 2 + (4*D-2)*(K-1) + 2*D*(K*K - K + 1) + mp.ceil(mp.log(C,2)) - 1)
		logical_qubits = 2*(2*C + D*K - 1 + D*(3*K + 2*K*K) + K)
		reaction_depth = grover_iterations*(2*mp.ceil(mp.log(C,2)) - 2 + 2*K*mp.log(K,2) - 2*K - 2*mp.log(K,2) + 4 + 2*(mp.ceil(mp.log(D,2)) + 1)*(K-1) + 2*mp.ceil(mp.log(mp.ceil(mp.log(C,2)), 2)) )
		active_volume = grover_iterations*((25 + 1.5*K + C_CCZ)*C + 4*(2*K*D + 4) + (4*D-2)*((K-1)*(39 + C_CCZ) + 7) + 2*D*(28*K*K - 42*K + 28 + (K*K - K + 1)*C_CCZ) + (mp.ceil(mp.log(C,2)) - 1)*(18 + C_CCZ))
	elif case == 'GaussSieve2':
		toffoli_count = grover_iterations*(C - 2 + (D+1)*(K-1) + D*(0.5*K*K - 1.5*K + 1) + mp.ceil(mp.log(C,2)) - 1)
		logical_qubits = 2*(2*C + D*K - 1 + D*(3*K + 2*K*K) + K)
		reaction_depth = grover_iterations*(2*mp.ceil(mp.log(C,2)) + 2*K*mp.log(K,2) - 2*K - 2*mp.log(K,2) + 1 + 2*(mp.ceil(mp.log(D,2)) + 1)*(K-1) + 2*mp.ceil(mp.log(mp.ceil(mp.log(C,2)), 2)) )
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


	dx =  mp.ceil(distance_baseline/4)
	dz =  mp.ceil(distance_baseline/8)
	dm =  mp.ceil(distance_baseline/8)
	dx2 = mp.ceil(distance_baseline/2)
	dz2 = mp.ceil(distance_baseline/4)
	dm2 = mp.ceil(distance_baseline/4)
	dx3 = distance_baseline
	dz3 = mp.ceil(distance_baseline/2)
	dm3 = mp.ceil(distance_baseline/2)
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
	
	

def resources_no_qram(D, C, M, case = 'NVSieve'):

	# Outputs number of physical qubits and circuit time required by a baseline and an active-volume architectures given dimension D, list size C, and number of solutions M for Grover's search
	# The QRAM cost is ignored

	if M == 0:
		grover_iterations = mp.ceil(9.2 * mp.log(1/delta, 3) * mp.sqrt(C))
	else:
		grover_iterations = mp.ceil(3.1 * mp.sqrt(C/M))


	if case == 'NVSieve':
		toffoli_count = grover_iterations*(2*D*(K-1) + D*(K*K - K + 1) + mp.ceil(mp.log(C,2)) - 1)
		logical_qubits = 2*(D*(3*K + 2*K*K) + K)
		reaction_depth = 2*grover_iterations*(K*mp.log(K,2) - K - mp.log(K,2) + 2 + (mp.ceil(mp.log(D,2)) + 2)*(K-1) + mp.ceil(mp.log(mp.ceil(mp.log(C,2)), 2)) )
		active_volume = grover_iterations*(2*D*((K-1)*(39 + C_CCZ) + 7) + D*(28*K*K - 42*K + 28 + (K*K - K + 1)*C_CCZ) + (mp.ceil(mp.log(C,2)) - 1)*(18 + C_CCZ))
	elif case == 'GaussSieve1':
		toffoli_count = grover_iterations*((4*D-2)*(K-1) + 2*D*(K*K - K + 1) + mp.ceil(mp.log(C,2)) - 1)
		logical_qubits = 2*(D*(3*K + 2*K*K) + K)
		reaction_depth = grover_iterations*(2*K*mp.log(K,2) - 2*K - 2*mp.log(K,2) + 4 + 2*(mp.ceil(mp.log(D,2)) + 1)*(K-1) + 2*mp.ceil(mp.log(mp.ceil(mp.log(C,2)), 2)) )
		active_volume = grover_iterations*(4*(2*K*D + 4) + (4*D-2)*((K-1)*(39 + C_CCZ) + 7) + 2*D*(28*K*K - 42*K + 28 + (K*K - K + 1)*C_CCZ) + (mp.ceil(mp.log(C,2)) - 1)*(18 + C_CCZ))
	elif case == 'GaussSieve2':
		toffoli_count = grover_iterations*((D+1)*(K-1) + D*(0.5*K*K - 1.5*K + 1) + mp.ceil(mp.log(C,2)) - 1)
		logical_qubits = 2*(D*(3*K + 2*K*K) + K)
		reaction_depth = grover_iterations*(2*K*mp.log(K,2) - 2*K - 2*mp.log(K,2) + 3 + 2*(mp.ceil(mp.log(D,2)) + 1)*(K-1) + 2*mp.ceil(mp.log(mp.ceil(mp.log(C,2)), 2)) )
		active_volume = grover_iterations*(4*(2*K*D + 4) + (D-1)*((K-1)*(39 + C_CCZ) + 7) + D*(20.25*K*K - 48.75*K + 32 + (0.5*K*K - 1.5*K + 1)*C_CCZ) + (mp.ceil(mp.log(C,2)) - 1)*(18 + C_CCZ))
		


	########################################################
	# Code distance


	baseline_spacetime_volume = 4 * (reaction_depth/2) * logical_qubits


	distance_baseline = 1
	error_probability_baseline = baseline_spacetime_volume * distance_baseline * 0.1 * mp.power(100*p_phy, (distance_baseline + 1)/2) 


	while error_probability_baseline > error:
		distance_baseline += 1
		error_probability_baseline = baseline_spacetime_volume * distance_baseline * 0.1 * mp.power(100*p_phy, (distance_baseline + 1)/2)
		
		
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


	dx =  mp.ceil(distance_baseline/4)
	dz =  mp.ceil(distance_baseline/8)
	dm =  mp.ceil(distance_baseline/8)
	dx2 = mp.ceil(distance_baseline/2)
	dz2 = mp.ceil(distance_baseline/4)
	dm2 = mp.ceil(distance_baseline/4)
	dx3 = distance_baseline
	dz3 = mp.ceil(distance_baseline/2)
	dm3 = mp.ceil(distance_baseline/2)
	nl1 = 4
	nl2 = 4
	
	if case == 'NVSieve':
		distillation_factories = 4 * max(6 * dm2 / (nl2 / 2), dm3) / (4 * distance_baseline) * (D * (K*K + K)/2)
	elif case == 'GaussSieve1':
		distillation_factories = 4 * max(6 * dm2 / (nl2 / 2), dm3) / (4 * distance_baseline) * (2*D * (K*K + K)/2)
	elif case == 'GaussSieve2':
		distillation_factories = 4 * max(6 * dm2 / (nl2 / 2), dm3) / (4 * distance_baseline) * (D * K/2)
	
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





reaction_depth_NVSieve = lambda C, D: 2 * mp.ceil(9.2 * mp.log(1/delta, 3) * mp.sqrt(C))*(mp.ceil(mp.log(C,2)) + K*mp.log(K,2) - K - mp.log(K,2) + 1 + (mp.ceil(mp.log(D,2)) + 2)*(K-1) + mp.ceil(mp.log(mp.ceil(mp.log(C,2)), 2)) )

reaction_depth_GaussSieve_search1 = lambda C, D: mp.ceil(9.2 * mp.log(1/delta, 3) * mp.sqrt(C))*(2*mp.ceil(mp.log(C,2)) + 2*K*mp.log(K,2) - 2*K - 2*mp.log(K,2) + 2 + 2*(mp.ceil(mp.log(D,2)) + 1)*(K-1) + 2*mp.ceil(mp.log(mp.ceil(mp.log(C,2)), 2)) )

reaction_depth_GaussSieve_search2 = lambda C, D: mp.ceil(9.2 * mp.log(1/delta, 3) * mp.sqrt(C))*(2*mp.ceil(mp.log(C,2)) + 2*K*mp.log(K,2) - 2*K - 2*mp.log(K,2) + 1 + 2*(mp.ceil(mp.log(D,2)) + 1)*(K-1) + 2*mp.ceil(mp.log(mp.ceil(mp.log(C,2)), 2)) )
	
def resources_nist(D, C, M, I, case = 'NVSieve'):

	# Outputs number of physical qubits and circuit time required by a baseline and an active-volume architectures given dimension D, list size C, and number of solutions M for Grover's search
	# The reaction limit of each Grover's search is capped at input I

	if M == 0:
		grover_iterations = mp.ceil(9.2 * mp.log(1/delta, 3) * mp.sqrt(C))
	else:
		grover_iterations = mp.ceil(3.1 * mp.sqrt(C/M))

	if case == 'NVSieve':
	
		if reaction_depth_NVSieve(C, D) < I:
		
			F = 1
			toffoli_count = grover_iterations*(C - 2 + 2*D*(K-1) + D*(K*K - K + 1) + mp.ceil(mp.log(C,2)) - 1)
			logical_qubits = 2*(2*C + D*K - 1 + D*(3*K + 2*K*K) + K)
			reaction_depth = 2*grover_iterations*(mp.ceil(mp.log(C,2)) + K*mp.log(K,2) - K - mp.log(K,2) + 1 + (mp.ceil(mp.log(D,2)) + 2)*(K-1) + mp.ceil(mp.log(mp.ceil(mp.log(C,2)), 2)) )
			active_volume = grover_iterations*((25 + 1.5*K + C_CCZ)*C + 2*D*((K-1)*(39 + C_CCZ) + 7) + D*(28*K*K - 42*K + 28 + (K*K - K + 1)*C_CCZ) + (mp.ceil(mp.log(C,2)) - 1)*(18 + C_CCZ))
			
		else:
			F = mp.floor( mp.findroot(lambda F: reaction_depth_NVSieve(C/F, D) - I, (0.05*9.2*mp.log(1/delta, 3))**2 * C/I**2, verify=False) )
			if reaction_depth_NVSieve(C/F, D) > I:
				while reaction_depth_NVSieve(C/F, D) > I:
					F *= 1.00001
			else:
				while reaction_depth_NVSieve(C/F, D) < I:
					F /= 1.00001

			F = mp.ceil(F)
			
#			print('Error:', reaction_depth_nvsieve(C/F) - I )
			
			reaction_depth = reaction_depth_NVSieve(C/F, D)
			toffoli_count = mp.ceil(9.2 * mp.log(1/delta, 3) * mp.sqrt(C/F))*(C/F - 2 + 2*D*(K-1) + D*(K*K - K + 1) + mp.ceil(mp.log(C/F,2)) - 1)
			logical_qubits = 2*(2*C/F + D*K - 1 + D*(3*K + 2*K*K) + K)
			active_volume = mp.ceil(9.2 * mp.log(1/delta, 3) * mp.sqrt(C/F))*((25 + 1.5*K + C_CCZ)*C/F + 2*D*((K-1)*(39 + C_CCZ) + 7) + D*(28*K*K - 42*K + 28 + (K*K - K + 1)*C_CCZ) + (mp.ceil(mp.log(C/F,2)) - 1)*(18 + C_CCZ))
			

	elif case == 'GaussSieve1':
	
		if reaction_depth_GaussSieve_search1(C, D) < I:
		
			F = 1
			toffoli_count = grover_iterations*(C - 2 + (4*D-2)*(K-1) + 2*D*(K*K - K + 1) + mp.ceil(mp.log(C,2)) - 1)
			logical_qubits = 2*(2*C + D*K - 1 + D*(3*K + 2*K*K) + K)
			reaction_depth = grover_iterations*(2*mp.ceil(mp.log(C,2)) + 2*K*mp.log(K,2) - 2*K - 2*mp.log(K,2) + 2 + 2*(mp.ceil(mp.log(D,2)) + 1)*(K-1) + 2*mp.ceil(mp.log(mp.ceil(mp.log(C,2)), 2)) )
			active_volume = grover_iterations*((25 + 1.5*K + C_CCZ)*C + 4*(2*K*D + 4) + (4*D-2)*((K-1)*(39 + C_CCZ) + 7) + 2*D*(28*K*K - 42*K + 28 + (K*K - K + 1)*C_CCZ) + (mp.ceil(mp.log(C,2)) - 1)*(18 + C_CCZ))

		else:
			F = mp.floor( mp.findroot(lambda F: reaction_depth_GaussSieve_search1(C/F, D) - I, (0.05*9.2*mp.log(1/delta, 3))**2 * C/I**2, verify=False) )
			if reaction_depth_GaussSieve_search1(C/F, D) > I:
				while reaction_depth_GaussSieve_search1(C/F, D) > I:
					F *= 1.00001
			else:
				while reaction_depth_GaussSieve_search1(C/F, D) < I:
					F /= 1.00001

			F = mp.ceil(F)
				
#			print('Error:', reaction_depth_GaussSieve_search1(C/F) - I )
			
			reaction_depth = reaction_depth_GaussSieve_search1(C/F, D)
			toffoli_count = mp.ceil(9.2 * mp.log(1/delta, 3) * mp.sqrt(C/F))*(C/F - 2 + (4*D-2)*(K-1) + 2*D*(K*K - K + 1) + mp.ceil(mp.log(C/F,2)) - 1)
			logical_qubits = 2*(2*C/F + D*K - 1 + D*(3*K + 2*K*K) + K)
			active_volume = mp.ceil(9.2 * mp.log(1/delta, 3) * mp.sqrt(C/F))*((25 + 1.5*K + C_CCZ)*C/F + 4*(2*K*D + 4) + (4*D-2)*((K-1)*(39 + C_CCZ) + 7) + 2*D*(28*K*K - 42*K + 28 + (K*K - K + 1)*C_CCZ) + (mp.ceil(mp.log(C/F,2)) - 1)*(18 + C_CCZ))
		
	elif case == 'GaussSieve2':
	
		if reaction_depth_GaussSieve_search2(C, D) < I:
		
			F = 1
			toffoli_count = grover_iterations*(C - 2 + (D+1)*(K-1) + D*(0.5*K*K - 1.5*K + 1) + mp.ceil(mp.log(C,2)) - 1)
			logical_qubits = 2*(2*C + D*K - 1 + D*(3*K + 2*K*K) + K)
			reaction_depth = grover_iterations*(2*mp.ceil(mp.log(C,2)) - 2 + 2*K*mp.log(K,2) - 2*K - 2*mp.log(K,2) + 3 + 2*(mp.ceil(mp.log(D,2)) + 1)*(K-1) + 2*mp.ceil(mp.log(mp.ceil(mp.log(C,2)), 2)) )
			active_volume = grover_iterations*((25 + 1.5*K + C_CCZ)*C + 4*(2*K*D + 4) + (D-1)*((K-1)*(39 + C_CCZ) + 7) + D*(20.25*K*K - 48.75*K + 32 + (0.5*K*K - 1.5*K + 1)*C_CCZ) + (mp.ceil(mp.log(C,2)) - 1)*(18 + C_CCZ))
			
		else:
			F = mp.floor( mp.findroot(lambda F: reaction_depth_GaussSieve_search2(C/F, D) - I, (0.05*9.2*mp.log(1/delta, 3))**2 * C/I**2, verify=False) )
			if reaction_depth_GaussSieve_search2(C/F, D) > I:
				while reaction_depth_GaussSieve_search2(C/F, D) > I:
					F *= 1.00001
			else:
				while reaction_depth_GaussSieve_search2(C/F, D) < I:
					F /= 1.00001
			
			F = mp.ceil(F)
				
#			print('Error:', reaction_depth_GaussSieve_search2(C/F) - I )
		
			reaction_depth = reaction_depth_GaussSieve_search2(C/F, D)
			toffoli_count = mp.ceil(9.2 * mp.log(1/delta, 3) * mp.sqrt(C/F))*(C/F - 2 + (D+1)*(K-1) + D*(0.5*K*K - 1.5*K + 1) + mp.ceil(mp.log(C/F,2)) - 1)
			logical_qubits = 2*(2*C/F + D*K - 1 + D*(3*K + 2*K*K) + K)
			active_volume = mp.ceil(9.2 * mp.log(1/delta, 3) * mp.sqrt(C/F))*((25 + 1.5*K + C_CCZ)*C/F + 4*(2*K*D + 4) + (D-1)*((K-1)*(39 + C_CCZ) + 7) + D*(20.25*K*K - 48.75*K + 32 + (0.5*K*K - 1.5*K + 1)*C_CCZ) + (mp.ceil(mp.log(C/F,2)) - 1)*(18 + C_CCZ))
		


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


	dx =  mp.ceil(distance_baseline/4)
	dz =  mp.ceil(distance_baseline/8)
	dm =  mp.ceil(distance_baseline/8)
	dx2 = mp.ceil(distance_baseline/2)
	dz2 = mp.ceil(distance_baseline/4)
	dm2 = mp.ceil(distance_baseline/4)
	dx3 = distance_baseline
	dz3 = mp.ceil(distance_baseline/2)
	dm3 = mp.ceil(distance_baseline/2)
	nl1 = 4
	nl2 = 4
	distillation_factories = 4 * max(6 * dm2 / (nl2 / 2), dm3) / (4 * distance_baseline) * C/F/ 2
	qubits_factories = distillation_factories * (2 * mp.ceil((3 * dx3 + dz3) * 3 * dx3 + nl2 * (dx2 + 4 * dz2) * dm3 / 2 + 20 * dm3 * dm3 + 2 * dx3 * dm3) + 2 * nl2 * mp.ceil((3 * dx2 + dz2) * 3 * dx2 + nl1 * ((dx + 4 * dz) * (3 * dx + dm2 / 2) + 2 * dm) + 20 * dm2 * dm2 + 2 * dx2 * dm2) )

	#factory = cost_of_three_level_8toccz(p_phy, dx, dz, dm, dx2, dz2, dm2, dx3, dz3, dm3, 4, 4)
	#distillation_factories[i] += [factory.distillation_time_in_cycles / (4 * distance_baseline[i][j]) * S / 2]
	#qubits_factories[i] += [factory.qubits * distillation_factories[i][j]]		
	#print('Toffoli error:', '{:.2e}'.format(error/toffoli_count[i][j]) )
	#print('Distillation error:', '{:.2e}'.format(factory.distilled_magic_state_error_rate) )

	#######################################################

	physical_qubits_baseline = 2 * distance_baseline * distance_baseline * logical_qubits + qubits_factories
	physical_qubits_active = distance_active * distance_active * logical_qubits
	time_baseline = F * 4 * distance_baseline * (reaction_depth/2) * code_cycle / 3600/24/365
	time_active = F * active_volume / logical_qubits / 2 * distance_active * code_cycle / 3600/24/365
	reaction_limit = F * reaction_depth * reaction_time / 3600/24/365
	
	return physical_qubits_baseline, physical_qubits_active, time_baseline, time_active, reaction_limit

