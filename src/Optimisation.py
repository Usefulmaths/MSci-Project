from differential_evolution import differential_evolution
from hamiltonian_functions import *
from math import log, sin, cos
from Fidelity import Fidelity

j = complex(0, 1)

class Optimisation:
	def __init__(self, ideal_gate, number_qubits_network, filename):
		self.ideal_gate = ideal_gate
		self.number_qubits_network = number_qubits_network
		self.filename = filename

		self.params = []
		self.minimum = 0
		self.fidelity_obj = Fidelity(self.ideal_gate, self.number_qubits_network)

	def set_params(self, params):
		self.params = params

	def get_number_params(self):
		return len(self.params)

	def get_closed_solution_fidelity(self):
		return Fidelity(self.ideal_gate, self.number_qubits_network).deterministic_gate_fidelity(self.params)

	def get_open_solution_fidelity(self):
		return Fidelity(self.ideal_gate, self.number_qubits_network).average_fidelity(self.params, 200)
	
	def optimise(self, func, iterations=100, number_of_agents=6, CR=0.5, F=0.5, bounds=[-10,10], hard_bounds=[-10,10]):
		return differential_evolution(func, iterations, number_of_agents, 30, CR, F, bounds, hard_bounds, self.filename, threshold_value=-0.99, thread=True)
