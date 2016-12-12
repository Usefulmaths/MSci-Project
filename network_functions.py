from qutip import *

## Number of qubits in the network.
def number_of_qubits():
	return 4	

## Array of spin matrices.
def get_spin_matrices():
	return [sigmax(), sigmay(), sigmaz()]

## Array of the different combinations of spin interactions.
def interaction_combination():
	return [[sigmax(), sigmax()], [sigmay(), sigmay()], [sigmaz(), sigmaz()]]

## Calculates all possible interactions between qubits.
def interaction_between_qubits_array(number_of_qubits):
    interactions = []

    for i in range(0, number_of_qubits - 1):
        for j in range(i + 1, number_of_qubits):
            interactions.append([i, j])
    return interactions

## Definition of the quantum half-adder.
def quantum_half_adder():
	return tensor(cnot(), qeye(2)) * toffoli()