from Network import Network
from Optimisation import Optimisation
from qutip import *

if __name__ == '__main__':

    # Define properties of the quantum network
    number_of_qubits = 4
    qubit_interations = [[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]]
    interaction_types = [[sigmax(), sigmax()], [sigmay(), sigmay()], [
        sigmaz(), sigmaz()]]

    # Define the gate you want to simulate
    ideal_gate = tensor(cnot(), qeye(2)) * toffoli()

    # Create the quantum network
    quantum_network = Network(
        number_of_qubits, qubit_interations, interaction_types)

    # Create the optimisation object between the network and ideal gate
    optimise_network = Optimisation(
        quantum_network, ideal_gate,
        "qha4_xxyyzz_de_output_vectorised2_approx")

    # Run the optimisation on an objective function of your choice
    average_gate_fidelity_function = optimise_network.fidelity_obj.average_gate_fidelity

    # Set up network dependencies
    optimise_network.fidelity_obj.set_basis_array()
    optimise_network.fidelity_obj.set_generate_dont_care_states()
    optimise_network.fidelity_obj.set_decomposed_gate()
    optimise_network.fidelity_obj.instantiate_xi_and_gv()

    # Set bounds for optimisation
    bounds = [[-20, 20]
              for params in range(quantum_network.number_of_params())]

    # Number of processors to run over
    number_of_processors = 4

    # Optimise Network
    optimise_network.optimise_de(
        average_gate_fidelity_function, bounds, number_of_processors)
