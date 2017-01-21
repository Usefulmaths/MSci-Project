from network_functions import *

N = number_of_qubits()
spin_matrices = get_spin_matrices()
interactions = interaction_between_qubits_array(N)
combination_of_interactions = interaction_combination()

# Calculates the hamiltonian of qubits without interactions
def hamiltonian_without_interaction_qubit(h):
    hamiltonian = 0
    h_iter = 0

    for i in range(N):
        operator = [qeye(2)] * N
        
        for j in range(len(spin_matrices)):
            operator[i] += h[h_iter] * spin_matrices[j]
            h_iter += 1

        hamiltonian += tensor(operator)/2
        
    return hamiltonian

# Calculates hamiltonian of interactions
def hamiltonian_interaction(J):
    hamiltonian = 0
    j_iter = 0
    for p in interactions: 
        
        temp = 0
        
        for S in combination_of_interactions: 
        
            OpChain = [qeye(2)]*N
            OpChain[p[0]] = S[0]
            OpChain[p[1]] = S[1]
            
            temp += J[j_iter]*tensor(OpChain)
            j_iter += 1

        hamiltonian += temp 

    return hamiltonian

def hamiltonian(params):
    J_ = params[0:len(interactions)*len(combination_of_interactions)]
    h_ = params[len(interactions)*len(combination_of_interactions):len(interactions)*len(combination_of_interactions) + len(N * spin_matrices)]
    return hamiltonian_interaction(J_) + hamiltonian_without_interaction_qubit(h_)
