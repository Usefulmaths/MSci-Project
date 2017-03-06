from math import log, sin, cos
from qutip import *
from time import time
j = complex(0, 1)

class Fidelity:

	def __init__(self, network, ideal_gate):
		self.ideal_gate = ideal_gate
		self.network = network

	def get_gate_dimension(self):
		return self.ideal_gate.shape[0]

	def number_region_qubits(self):
		return int(log(self.get_gate_dimension(), 2))


	def set_generate_dont_care_states(self):
		h = self.network.get_number_qubits() - self.number_region_qubits()
		states_dont_care = [sin(0.847) * basis(2,1) + cos(0.847) * basis(2,0)] * (h)
		self.states_dont_care = tensor(states_dont_care)

	# Decomposes the ideal gate into the positions of the non-zero elements
	def get_gate(self):
		s = []

		rows = self.ideal_gate.shape[0]
		columns = self.ideal_gate.shape[1]

		for i in range(rows):
			for j in range(columns):
				if self.ideal_gate[i][0][j] != 0:
					s.append([i, j])
		return s

	def set_decomposed_gate(self):
		self.decomposed_gate = self.get_gate()

	def get_bin(self, a):
		s = bin(a)[2:]
		l = len(s)
		if l < 4:
			r = '0' * (4 - l)
			s = r + s

		return s


	def get_basis(self, a):
		if a == 0:
			B = [basis(2, 0)] * self.number_region_qubits()
			return tensor(B)

		c = self.get_bin(a)
		return tensor(basis(2, int(c[-3])), basis(2, int(c[-2])), basis(2, int(c[-1])))

	def set_basis_array(self):
		basis = []
		for i in range(8):
			basis = [self.get_basis(i) for i in range(self.get_gate_dimension())]
		self.basis_array = basis

	def basis_8(self, i):
		a = basis(8, i)
		a.dims = [[2, 2, 2], [1, 1, 1]]
		return a


	def instantiate_xi_and_gv(self):
		xi = 0 

		for i in range(8):
			xi += 1./8**0.5 * tensor(self.basis_8(i), self.basis_8(i), self.states_dont_care)

		gv = 0

		for i in range(8):
			for p in range(8):
				gv += self.ideal_gate.data[i, p] * tensor(self.basis_8(p), self.basis_8(i))
				
		self.xi = xi
		self.gv = gv

	def average_gate_fidelity(self, params):
		H = self.network.hamiltonian(params)
		U = (-j * H).expm()
		D = self.ideal_gate.shape[0]

		e = qeye(D)
		e.dims = [[2, 2, 2], [2, 2, 2]]

		xi2 = tensor(e, U) * self.xi

		rho = xi2.ptrace([0, 1, 2, 3, 4, 5])

		overlap = self.gv.dag() * rho * self.gv

		return -(1./(D + 1) + 1./(D + 1) * abs(overlap[0][0][0])), params

	def deterministic_gate_fidelity(self, params):
		H = self.network.hamiltonian(params)
		U = (-j * H).expm()
		Udag = U.dag()
		G = self.ideal_gate

		gate_element_positions = self.decomposed_gate
		gate_dimension = float(self.get_gate_dimension())
		Fid = 1./(gate_dimension + 1)

		for x in range(len(gate_element_positions)):
		
			bra_i = self.basis_array[gate_element_positions[x][0]].dag()
			ket_k = self.basis_array[gate_element_positions[x][1]]
			for y in range(len(gate_element_positions)):
				ket_j = self.basis_array[gate_element_positions[y][0]]
				bra_l = self.basis_array[gate_element_positions[y][1]].dag()
				epsilon = U * tensor(ket_k * bra_l, self.states_dont_care * self.states_dont_care.dag()) * Udag
				
				eps_ijkl = bra_i * epsilon.ptrace([0, 1, 2]) * ket_j

				g_star_ik = G[gate_element_positions[x][0], gate_element_positions[x][1]].conjugate()
				g_jl = G[gate_element_positions[y][0], gate_element_positions[y][1]]

				fidStep = (1. / (gate_dimension * (gate_dimension + 1))) * g_star_ik * eps_ijkl * g_jl
				Fid += fidStep[0][0][0]
		return -abs(Fid)

	def likelihood(self, params, rho_0):
		rho = tensor(rho_0, self.generate_dont_care_states())
		H = self.network.hamiltonian(params)
		U = (-j * H).expm()

		A = (self.ideal_gate * rho_0) * (self.ideal_gate * rho_0).dag()
		B = ((U * rho) * (U * rho).dag()).ptrace([0, 1, 2])

		return abs((A * B).tr())

	def average_fidelity(self, params, number_iters):
		H = self.network.hamiltonian(params)
		U = (-j * H).expm()
		G = self.ideal_gate

		total_fid = 0

		for iters in range(number_iters):
			psi0 = tensor(rand_ket(N = 2), rand_ket(N = 2), rand_ket(N = 2))
			psi_state = tensor(psi0, self.states_dont_care)

			epsilon = (U * ket2dm(psi_state) * U.dag()).ptrace([0, 1, 2])

			fid_contribution = psi0.dag() * G.dag() * epsilon * G * psi0
			total_fid += fid_contribution[0][0][0]

		return abs(total_fid) / number_iters


