from Network import Network
from Fidelity import Fidelity
from Optimisation import Optimisation
from qutip import *

if __name__ == '__main__':

	''' Define properties of the quantum network '''
	number_of_qubits = 5
	qubit_interations = [[0, 1], [0, 2], [0, 3], [0, 4], [1, 2], [1, 3], [1, 4], [2, 3], [2, 4], [3, 4]]
	interaction_types = [[sigmax(), sigmax()], [sigmay(), sigmay()], [sigmaz(), sigmaz()]]

	''' Define the gate you want to simulate '''
	ideal_gate = tensor(cnot(), qeye(2)) * toffoli()

	''' Create the quantum network '''
	quantum_network = Network(number_of_qubits, qubit_interations, interaction_types)

	''' Create the optimisation object between the network and ideal gate '''
	optimise_network = Optimisation(quantum_network, ideal_gate, "QHA_output")

	''' Run the optimisation on an objective function of your choice '''
	likelihood_function = optimise_network.fidelity_obj.likelihood

	''' SGD parameters '''
	step_size = 0.0008
	delta = 0.0001
	learning_rate_convergence = 300
	iterations = 1000
	iteration_length = 1600

	''' Start the optimisation '''
	optimise_network.optimise_sgd(likelihood_function, step_size, delta, learning_rate_convergence, iterations, iteration_length)

