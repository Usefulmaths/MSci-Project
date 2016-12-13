from differential_evolution import differential_evolution
from network_functions import *
from hamiltonian_functions import *
from fidelity_functions import *
import time
from qutip import *
import numpy as np
import time

# Optimised parameters for hamiltonian.
lambd8 = [-6.765761074623748e-05, 4.8925584986167264e-05, -5.0, -0.0006515643909790049, 8.03265659236998e-05, 0.41522759546199006, 0.001128760474954534, -0.001959500831485896, 0.3888728922665736, -0.00024110569483596955, 0.002339789280788476, -1.3602831171903347, 0.0280247974930564, 0.012128822243884983, 4.535040541723734e-05, 6.05918798004767e-06, -0.00010168566700763935, -5.452112800657558e-05, 0.003144002388730177, 0.00047905546995522376, -3.223377696727116, -0.05755801030587989, -0.0018161248949912757, -4.9432613968343935, 1.3412674048237916, -2.841347967983532, -3.5512444524655473, 2.797108088956102, -0.03665068450783443, 4.750717834497733]

## This function carries out differential evolution on the fidelity.
#print(differential_evolution(Fidelity, 10000, 16, len(interactions)*len(combination_of_interactions) + N * len(spin_matrices), 0.5, 0.5, [-1.2, 1.2], [-5, 5], -0.99))

# Here are some tests to compare the theoretical gate operation and the optimised gate operation on all binary states.
H = hamiltonian(lambd8)
U = (-j * H).expm()
G = quantum_half_adder()

for i in range(0, 2):
	for j in range(0, 2):
		for k in range(0, 2):

			psi0 = tensor(basis(2, i), basis(2, j), basis(2, k))
			psi = tensor(psi0, sin(0.847) * basis(2,1) + cos(0.847) * basis(2,0))

			rho = psi * psi.dag()
			trans_state = (U * rho * U.dag()).ptrace([0, 1, 2])
			trans_G_state = G * rho.ptrace([0, 1, 2]) * G.dag()
			print(fidelity(trans_G_state, trans_state))


