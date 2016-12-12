from random import uniform
from random import randrange

def generate_agents(N, number_of_parameters, bounds):
	agents = []
	for i in range(N):
		agent = []
		for j in range(number_of_parameters):
			agent.append(uniform(bounds[0], bounds[1]))
		agents.append(agent)
	return agents


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

def differential_evolution(func, iterations, number_of_agents, number_of_parameters, CR, F, bounds, threshold_value=10e-100):
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

			y = []
			for k in range(number_of_parameters):
				r_k = uniform(0, 1)
				if(r_k < CR or k == R):
					y.append(a[k] + F * (b[k] - c[k]))
				else:
					y.append(x[k])

			if(func(y) < func(x)):
				agents[j] = y
			else:
				continue

		if(iters % 10 == 0):
			for agent in agents:
				f = func(agent)
				if(f < minimum):
					minimum = f
					optimised_agent = agent
				else: 
					continue

			print("Minimum after " + str((iters + 1)) + " iterations: " + str(minimum))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
		iters += 1
		if(minimum < threshold_value):
			print(iters, optimised_agent)
		if(iters > iterations):
			break
			

		

		
	return minimum, optimised_agent, iters