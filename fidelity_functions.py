from qutip import *
from hamiltonian_functions import *
from network_functions import *
from math import log
from math import sin
from math import cos

j = complex(0, 1)
G = quantum_half_adder()
gate_dimension = G.shape[0]

# Dimension of states we care about (Region qubits)
dimension_state_care = int(log(gate_dimension, 2))

# Generates a tensor state and identity for the states we don't care about (ancillae qubits)
def generate_dont_care_states(N, dimension_state_care):
    h = N - dimension_state_care
    
    states_dont_care = [sin(0.847) * basis(2,1) + cos(0.847) * basis(2,0)] * (h)
    identity_dont_care = [qeye(2)] * (h)
    return tensor(states_dont_care), tensor(identity_dont_care)

# Extracts non zero element positions of the gate
def getGate(G): 
    
    non_zero_elements = []

    rows = G.shape[0]
    colums = G.shape[1]
    
    for i in range(rows):
        for j in range(colums):
            if G[i][0][j] != 0: 
                non_zero_elements.append([i,j])
    return non_zero_elements

# Converts a decimal number into its binary representation (computation basis).    
def get_bin(decimal_number): 
    binary_rep = bin(decimal_number)[2:]
    length_binary_rep = len(binary_rep)
    if length_binary_rep < 4: 
        append_zeros = '0' * (4 - length_binary_rep)
        binary_rep = append_zeros + binary_rep
        
    return binary_rep

# Converts decimal numbers into quantum basis states.
def getBasis(decimal_number):  
    if decimal_number == 0:
        B = [basis(2, 0)]*(dimension_state_care)
        return tensor(B)
    
    binary_rep = get_bin(decimal_number)
    if gate_dimension == 8:
        return tensor(basis(2, int(binary_rep[0])), basis(2, int(binary_rep[1])), basis(2, int(binary_rep[2])))

states_dont_care, identity_dont_care = generate_dont_care_states(N, dimension_state_care)

# Calculates the average fidelity.
def Fidelity (J):  
    gate_element_positions = getGate(G)
    H = hamiltonian(J)
    Fid = 1/(gate_dimension + 1)
    U = (-j*H).expm()
    Udag = U.dag()

    
    for x in gate_element_positions: 
        bra_i = getBasis(x[0]).dag()
        ket_k = getBasis(x[1])
        for y in gate_element_positions:
            ket_j = getBasis(y[0])
            bra_l = getBasis(y[1]).dag()
            
            Epsilon = U*tensor(ket_k*bra_l, states_dont_care * states_dont_care.dag()) * Udag
            Eps_ijkl = bra_i*(Epsilon.ptrace([0,1,2]))*ket_j
            
            Gstar_ik = G[x[0],x[1]].conjugate()
            G_jl = G[y[0],y[1]]
            
            fidStep = (1/(gate_dimension*(gate_dimension+1)))*Gstar_ik*Eps_ijkl*G_jl     
            Fid += abs(fidStep[0][0][0])
            
    return -Fid