from differential_evolution import adaptive_differential_evolution
from Fidelity import Fidelity
from qutip import *
import numpy as np
from random import randrange

'''
This class provides methods to optimise a quantum network to simulate a 
given ideal quantum gate.
'''

j = complex(0, 1)

class Optimisation:

	def __init__(self, network, ideal_gate, filename):
		self.ideal_gate = ideal_gate
		self.network = network
		self.filename = filename

		''' Default random parameters to avoid errors. '''
		self.params = [float(randrange(-2000, 2000))/100 for param in range(network.number_of_params())]

		''' Fidelity object between the network and ideal gate '''
		self.fidelity_obj = Fidelity(network, ideal_gate)

	'''Set initial parameters for the network'''
	def set_params(self, params):
		self.params = params

	'''Optimise the network via differential evolution'''
	def optimise_de(self, func, bounds, number_of_processors, iterations=10000000, CR=0.5, F=0.9):
		hard_bounds = [bound for bound in bounds]
		return adaptive_differential_evolution(func, iterations, number_of_processors, self.network.number_of_params(), CR, F, bounds, hard_bounds, self.filename, threshold_value=-0.99)

	def write_to_file(self, file_name, iters, value, params):
		file_open = open(file_name + ".txt", "a")
		file_open.write("\n" + str(iters) + "," + str(value))
		file_open.close()	
