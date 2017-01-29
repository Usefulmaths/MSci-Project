from differential_evolution import differential_evolution
from math import log, sin, cos, sqrt
from Fidelity import Fidelity
from qutip import *
from random import randrange
import numpy as np
from multiprocessing import Process, Manager
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
		self.params = [float(randrange(-1000, 1000))/1000 for param in range(network.number_of_params())]

		''' Fidelity object between the network and ideal gate '''
		self.fidelity_obj = Fidelity(network, ideal_gate)

	'''Set initial parameters for the network'''
	def set_params(self, params):
		self.params = params

	'''Optimise the network via differential evolution'''
	def optimise_de(self, func, iterations=10000, number_of_agents=8, CR=0.5, F=0.5):
		bounds = [[-5, 5] for params in range(self.network.number_of_params())]
		hard_bounds = [bound for bound in bounds]
		return differential_evolution(func, iterations, number_of_agents, self.network.number_of_params(), CR, F, bounds, hard_bounds, self.filename, threshold_value=-0.99)

	def write_to_file(self, file_name, iters, value, params):
		file_open = open(file_name + ".txt", "a")
		file_open.write("\n" + str(iters) + "," + str(value))
		file_open.close()	


	def gradient(self, func, params_increment, rho_0, delta, par, return_dict):
		return_dict[par] = (func(params_increment, rho_0) - func(params, rho_0))/delta


	'''STILL IN WORKING PROGRESS'''
	def optimise_sgd(self, func, step, delta, learning_rate_convergence, iterations, iteration_length):
		print("Beginning Optimisation")
		time = step

		for i in range(iterations):
			walk = step * learning_rate_convergence * (i + 1)
			while True:
				rho_0 = tensor(rand_ket(N = 2), rand_ket(N = 2), rand_ket(N = 2))

				index = randrange(len(self.params))
				params_increment = [x for x in self.params]
				params_increment[index] += delta
				grad = (func(params_increment, rho_0) - func(self.params, rho_0))/delta

				learning_rate = 1/sqrt(walk)

				self.params[index] = self.params[index] + learning_rate * grad

				function_evaluation = func(self.params, rho_0)				
				self.write_to_file(self.filename, int(time/step), function_evaluation, self.params)

			
				if(int(time/step) % 10 < 1):
					print(function_evaluation, learning_rate, time, int(time/step))
					
				time += step
				walk += step

				if(int(time/step) > iteration_length * (i + 1)):
					print('HIT')
					print(self.params)
					break

		print(self.params)

'''
		for i in range(100):
			walk = step * 500 * i

			while True:
				walk += step
				time += step
				rho_0 = tensor(rand_ket(N = 2), rand_ket(N = 2), rand_ket(N = 2))
				
				index = randrange(len(self.params))
				params_increment = [x for x in self.params]
				params_increment[index] += delta
				grad = (func(params_increment) - func(params))/delta
				self.params[index] = self.params[index] + grad/sqrt(walk)
#				grad_vec = []
#				jobs = []
#
#				manager = Manager()
#				return_dict = manager.dict()
#				for par in range(len(self.params)):
#					params_increment = [u for u in self.params]
#					params_increment[par] += delta
#					p = Process(target=self.gradient, args=(func, params_increment, rho_0, delta, par, return_dict,))
#					jobs.append(p)
#					p.start()
#
#				for job in jobs:
#					job.join()
#
#				grad_vec = np.array(return_dict.values(), dtype=float)
#
#				grad_vec = [m/sqrt(walk) for m in grad_vec]
#				self.params = np.add(self.params, grad_vec)
				
				s = func(params)				
				self.write_to_file(self.filename, iters, s)
				iters += 1
				print(walk)		
				if(iters > iterations * (i + 1)):
					print(self.params)
					print(walk)
					break
		print(self.params)
		'''
