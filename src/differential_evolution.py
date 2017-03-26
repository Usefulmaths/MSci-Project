import numpy as np
from random import uniform
from random import randrange
from multiprocessing import Process, Manager
from mpi4py import MPI
from time import time

# Generates an initial population of agents, given bounds.
def generate_agents(N, number_of_parameters, bounds):
	agents = []
	for i in range(N):
		agent = []
		for j in range(number_of_parameters):
			agent.append(uniform(bounds[j][0], bounds[j][1]))

		
		agents.append(agent)
	return agents

# Finds three unique random indexes in an array.
def three_agents(x_index, population):
	index1 = x_index
	index2 = x_index
	index3 = x_index
	
	while(x_index == index1):
		index1 = randrange(0, len(population))

	while(x_index == index2 or index1 == index2):
		index2 = randrange(0, len(population))

	while(x_index == index3 or index1 == index3 or index2 == index3):
		index3 = randrange(0, len(population))

	return index1, index2, index3

def write_to_file(file_name, iters, minimum, agent, time):
	file_open = open(file_name + ".txt", "a")
	file_open.write(str(iters) + ", " + str(minimum) + ", " + str(time) + ", " + str(agent) + "\n")
	file_open.close()	

def agent_procedure(func, agents, number_of_parameters, hard_bounds, F, CR, j):
	x = agents[j]
	a_index, b_index, c_index = three_agents(j, agents)
	
	a = agents[a_index]
	b = agents[b_index]
	c = agents[c_index]

	R = randrange(0, number_of_parameters)

	# Adaptive Parameters
	r1 = uniform(0, 1)
	r2 = uniform(0, 1)
	r3 = uniform(0, 1)
	r4 = uniform(0, 1)
	mu_l = 0.1
	mu_u = 0.9
	kappa1 = 0.1
	kappa2 = 0.1
	if (r2 < kappa1):
		F = mu_l + r1 * mu_u
	if(r4 < kappa2):
		CR = r3

	# Mutation and Recombination
	y = []
	for k in range(number_of_parameters):
		r_k = uniform(0, 1)
		if(r_k < CR or k == R):
			y.append(a[k] + F * (b[k] - c[k]))
		else:
			y.append(x[k])

	# Bound parameters within a region.
	for y_i in range(len(y)):
		if(y[y_i] > hard_bounds[y_i][1]):
			y[y_i] = uniform(hard_bounds[y_i][0], hard_bounds[y_i][1])
		elif(y[y_i] < hard_bounds[y_i][0]):
			y[y_i] = uniform(hard_bounds[y_i][0], hard_bounds[y_i][1])

	eval_y, _ = func(y)
	eval_x, _ = func(x)
	if(eval_y < eval_x):
		return y
	else:
		return x

def adaptive_differential_evolution(func, iterations, number_of_agents, number_of_parameters, CR, F, bounds, hard_bounds, file_name, threshold_value=10e100, thread=False):
	comm = MPI.COMM_WORLD
	rank = comm.Get_rank()
	size = comm.Get_size()
	minimum = 10E10
	optimised_agent = []
	iters = 0

	if(rank == 0):
		agents = generate_agents(number_of_agents, number_of_parameters, bounds)
		agents_array = [agents for i in range(size)]
	else:
		agents_array = []

	agents = comm.scatter(agents_array, root=0)
	

	while(True):
		start = time()
		data = agent_procedure(func, agents, number_of_parameters, hard_bounds, F, CR, rank)
		comm.Barrier()
		end = time()

		data = comm.gather(data, root=0)

		if(rank == 0):
			agents = data
			iters += 1
			print(iters)
			agents_array = [agents for i in range(size)]

		else:
			agents_array = []

		comm.Barrier()
		
		# Update and print the minimum every 10 iterations to keep track.

		iters = comm.bcast(iters, root=0)
		if(iters % 50 == 0):
			if(rank == 0):
				data = np.array_split(data, size)
			else:
				data = []
			agents_to_compare = comm.scatter(data, root=0)
			f, params = func(agents_to_compare[0])
			evaluated_agents = comm.gather((f, params), root=0)
			
			if(rank == 0):
				minimum, optimised_agent = min(evaluated_agents, key = lambda t: t[0])
				print(minimum, optimised_agent)
				write_to_file(file_name, iters, minimum, optimised_agent, end - start)  

		agents = comm.scatter(agents_array, root=0)
		comm.Barrier()

