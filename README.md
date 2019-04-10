# MSci Project
## Quantum gate learning in engineered qubit networks: quantum half-adder

The source code required to optimise a qubit network, using self-adaptive differential evolution, to perform a given quantum gate within its unmodulated dynamics.

The code was used to optimised the pairwise Ising interactions in a quantum network of four qubits to implement a quantum half-adder within the network's time evolution. This procedure was successful, resulting in a quantum half-adder with a fidelity of ~0.98. 

Additionally, this method was able to encode a Toffoli gate within the unmodulated dynamics of a four qubit network with a fidelity of 0.998.

The code can be used to construct an arbitary number of qubits in a network to implement a quantum gate of your choice. Note that this is an incredibly expensive task and will require parallelisation over a large cluster. By specifying the number of cores you would like to use, this code automatically implements a Message Passing Interface parallelisation of differential evolution over those cores.

Please see the following published article for further discussion about the methods used and the results obtained:
https://journals.aps.org/pra/abstract/10.1103/PhysRevA.97.062321
