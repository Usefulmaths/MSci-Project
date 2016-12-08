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
	
	while(not ((x_index != index1) and (x_index != index2) and (x_index != index3) and (index1 != index2) and (index1 != index3) and (index2 != index3))):
		index1 = randrange(0, len(population))
		index2 = randrange(0, len(population))
		index3 = randrange(0, len(population))

	return index1, index2, index3

def differential_evolution(func, number_of_agents, number_of_parameters, CR, F, bounds):
	agents = generate_agents(number_of_agents, number_of_parameters, bounds)

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

	minimum = 1E10
	optimised_agent = []
	for agent in agents:
		f = func(agent)
		if(f < minimum):
			minimum = f
			optimised_agent = agent
		else: 
			continue

	return minimum, optimised_agent

def func(x):
	return (x[0] - 2)**2 + 2
