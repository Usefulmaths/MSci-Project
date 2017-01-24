from differential_evolution import differential_evolution
from math import log, sin, cos, sqrt
from Fidelity import Fidelity
from qutip import *
from random import randrange
import numpy as np
from multiprocessing import Process, Manager
from random import randrange

j = complex(0, 1)

class Optimisation:
	def __init__(self, network, ideal_gate, filename):
		self.ideal_gate = ideal_gate
		self.network = network
		self.filename = filename

		self.params = [float(randrange(-1000, 1000))/100 for param in range(network.number_of_params())]
		self.minimum = 0
		self.fidelity_obj = Fidelity(self.ideal_gate, self.network.get_number_qubits)

	def set_params(self, params):
		self.params = params

	def optimise_de(self, func, iterations=10000, number_of_agents=4, CR=0.5, F=0.5, bounds=[-10,10], hard_bounds=[-10,10]):
		return differential_evolution(func, iterations, number_of_agents, self.network.number_of_params(), CR, F, bounds, hard_bounds, self.filename, threshold_value=-0.99)

	def write_to_file(self, file_name, iters, value):
		file_open = open(file_name + ".txt", "a")
		file_open.write("\n" + str(iters) + "," + str(value))
		file_open.close()	

	def gradient(self, func, params_increment, rho_0, delta, par, return_dict):
		return_dict[par] = (func(params_increment, rho_0) - func(self.params, rho_0))/delta


	def optimise_sgd(self, func, iterations, delta):
		step = 1./iterations
		iters = 1

		for i in range(100):
			walk = step * 4000 * i

			while True:
				walk += step
				rho_0 = tensor(rand_ket(N = 2), rand_ket(N = 2), rand_ket(N = 2))

				grad_vec = []
				jobs = []

				manager = Manager()
				return_dict = manager.dict()
				for par in range(len(self.params)):
					params_increment = [u for u in self.params]
					params_increment[par] += delta
					p = Process(target=self.gradient, args=(func, params_increment, rho_0, delta, par, return_dict,))
					jobs.append(p)
					p.start()

				for job in jobs:
					job.join()

				grad_vec = np.array(return_dict.values(), dtype=float)

				grad_vec = [m/sqrt(walk) for m in grad_vec]
				self.params = np.add(self.params, grad_vec)

				
				self.write_to_file("sgd_toffoli2", iters, func(self.params, rho_0))
				iters += 1

				print(func(self.params, rho_0))

				if(iters > iterations * (i + 1)):
					print(self.params)
					break

		print(self.params)