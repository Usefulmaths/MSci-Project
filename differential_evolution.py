from random import uniform
from random import randrange
import time
import multiprocessing
import operator

# Generates an initial population of agents, given bounds.
def generate_agents(N, number_of_parameters, bounds):
	agents = []
	for i in range(N):
		agent = []
		for j in range(number_of_parameters):
			agent.append(uniform(bounds[0], bounds[1]))
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

def write_to_file(file_name, iters, minimum, agent):
	file_open = open(file_name + ".txt", "a")
	file_open.write(str(iters) + ", " + str(minimum) + ", " + str(agent) + "\n")
	file_open.close()	


def thread_minimum(func, agents, limits, return_dict):
	minimum = 10E10
	optimised_index = -1

	for n in range(limits[0], limits[1] + 1):
		f = func(agents[n])
		if(f < minimum):
			minimum = f
			optimised_index = n
			return_dict[n] = minimum


def find_maximum_index(return_dict):  
     v=list(return_dict.values())
     k=list(return_dict.keys())
     return k[v.index(max(v))]


def differential_evolution(func, iterations, number_of_agents, number_of_parameters, CR, F, bounds, hard_bounds, file_name, threshold_value=10e100):
	# Initialisation
	agents = generate_agents(number_of_agents, number_of_parameters, bounds)

	optimised_agent = []
	minimum = 10E10
	iters = 0

	while(True):
		for j in range(len(agents)):
			x = agents[j]
			a_index, b_index, c_index = three_agents(j, agents)
			
			a = agents[a_index]
			b = agents[b_index]
			c = agents[c_index]

			R = randrange(0, number_of_parameters)

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
					if(y[y_i] > hard_bounds[1]):
						y[y_i] = hard_bounds[1]

					elif(y[y_i] < hard_bounds[0]):
						y[y_i] = hard_bounds[0]

					else:
						continue

			# Selection
			if(func(y) < func(x)):
				agents[j] = y
			else:
				continue

		# Update and print the minimum every 10 iterations to keep track.
		if(iters % 10 == 0):
			manager = multiprocessing.Manager()
			return_dict = manager.dict()

			jobs = []
			for i in range(0, 8, 2):
				p = multiprocessing.Process(target=thread_minimum, args=(func, agents, [i, i + 1], return_dict))
				jobs.append(p)
				p.start()

			for p in jobs:
				p.join()

			best_index = find_maximum_index(return_dict)
			minimum = return_dict[best_index]
			optimised_agent = agents[best_index]

			print("Minimum after " + str((iters + 1)) + " iterations: " + str(minimum))  
			write_to_file(file_name, iters, minimum, optimised_agent)   

		iters += 1

		# If minimum is less than threshold_value, print details and continue running the code.
		if(minimum < threshold_value):
			print(minimum, optimised_agent, iters)

		# If iteration criteria is met, return values and end code.
		if(iters > iterations):
			return minimum, optimised_agent, iters

		
	return minimum, optimised_agent, iters