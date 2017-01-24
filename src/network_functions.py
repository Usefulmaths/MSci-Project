from qutip import *

## Definition of the quantum half-adder.
def quantum_half_adder():
	return tensor(cnot(), qeye(2)) * toffoli()