from qutip import *
import numpy as np

spin_matrices = [sigmax(), sigmay(), sigmaz()]

class Network: 

	def __init__(self, number_qubits, interactions):
		self.number_qubits = number_qubits
		self.interactions = interactions

	def get_interactions(self):
		return self.interactions

	def get_number_qubits(self):
		return self.number_qubits

	def number_of_magnetic_params(self):
		return len(spin_matrices) * self.number_qubits
	
	def number_of_interaction_params(self):
		return len(spin_matrices) * len(self.interactions)

	def number_of_params(self):
		return self.number_of_magnetic_params() + self.number_of_interaction_params()

	def hamiltonian(self, params):
		J = params[:self.number_of_interaction_params()]
		h = params[self.number_of_interaction_params():self.number_of_params()]
		H = 0
		
		j_iter = 0
		for interaction in self.get_interactions():
			H_contribution = 0
			for interaction_type in spin_matrices:
				OpChain = [qeye(2)] * self.get_number_qubits()
				OpChain[interaction[0]] = interaction_type
				OpChain[interaction[1]] = interaction_type

				H_contribution += J[j_iter] * tensor(OpChain)
		
				j_iter += 1 

			H += 1./4 * H_contribution

		h_iter = 0
		for qubit_number in range(self.get_number_qubits()):
			H_contribution = 0
			for interaction_type in spin_matrices:
				OpChain = [qeye(2)] * self.get_number_qubits()
				OpChain[qubit_number] = interaction_type

				H_contribution += h[h_iter] * tensor(OpChain)

				h_iter += 1

			H += 1./2 * H_contribution

		return H
