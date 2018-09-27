from qutip import *

# Define the spin matrices as a global variable.
spin_matrices = [sigmax(), sigmay(), sigmaz()]


class Network:
    '''
    A class that represents a quantum network.

    Specify the number of numbers, which qubits will be interactions,
    and the type of interactions they will experience.
    '''

    def __init__(self, number_qubits, interactions, interaction_types):
        self.number_qubits = number_qubits
        self.interactions = interactions
        self.interaction_types = interaction_types

    def get_interactions(self):
        '''
        Returns:
                iteractions: the array of interactions in the network.
        '''
        iteractions = self.interactions

        return iteractions

    def get_number_qubits(self):
        '''
        Returns:
                number_of_qubits: the number of qubits in the quantum
                                  network.
        '''
        number_of_qubits = self.number_qubits

        return number_of_qubits

    def number_of_magnetic_params(self):
        '''
        Returns:
                number_of_h_params: the number of magnetic spin parameters, h.
        '''
        number_of_h_params = len(spin_matrices) * self.number_qubits
        return number_of_h_params

    def number_of_interaction_params(self):
        '''
        Returns:
                number_of_j_params: the number of pairwise spin
                                    interactions, J.
        '''
        number_of_j_params = len(
            self.interaction_types) * len(self.interactions)

        return number_of_j_params

    def number_of_params(self):
        '''
        Returns:
                number_of_parameters: the total number of parameters in the
                                  quantum network
        '''
        number_of_parameters = self.number_of_magnetic_params(
        ) + self.number_of_interaction_params()

        return number_of_parameters

    def hamiltonian(self, params):
        '''
        Calculates the total hamiltonian of the network.

        Arguments:
                params: the number of parameters in the network.

        Returns:
                H: the total hamiltonian of the quantum system.
        '''
        J = params[:self.number_of_interaction_params()]
        h = params[self.number_of_interaction_params():self.number_of_params()]
        H = 0

        j_iter = 0
        for interaction in self.get_interactions():
            H_contribution = 0
            for interaction_type in self.interaction_types:
                OpChain = [qeye(2)] * self.get_number_qubits()
                OpChain[interaction[0]] = interaction_type[0]
                OpChain[interaction[1]] = interaction_type[1]

                H_contribution += J[j_iter] * tensor(OpChain)

                j_iter += 1

            H += 1. / 4 * H_contribution

        h_iter = 0
        for qubit_number in range(self.get_number_qubits()):
            H_contribution = 0
            for interaction_type in spin_matrices:
                OpChain = [qeye(2)] * self.get_number_qubits()
                OpChain[qubit_number] = interaction_type

                H_contribution += h[h_iter] * tensor(OpChain)

                h_iter += 1

            H += 1. / 2 * H_contribution

        return H
