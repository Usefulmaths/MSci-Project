from qutip import *
from hamiltonian_functions import *
from network_functions import *
from math import log
from math import sin
from math import cos
import time
import threading
import multiprocessing
import queue

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

gate_element_positions = getGate(G)


def gate_operators(x, y):
    Gstar_ik = G[x[0], x[1]].conjugate()
    G_jl = G[y[0], y[1]]
    return Gstar_ik, G_jl

def fidelity_contribution(gate_dimension, Gstar, G, Eps):
    return (1/(gate_dimension*(gate_dimension+1)))*Gstar*Eps*G     

def evolve_system(ket_k, bra_l, states_dont_care, U, Udag):
    return U * tensor(ket_k * bra_l, states_dont_care *states_dont_care.dag()) * Udag

def get_eps_ijkl(bra_i, ket_j, Epsilon):
    return bra_i * Epsilon.ptrace([0, 1, 2]) * ket_j

def increment_of_fidelity(x, y, bra_i, ket_k, gate_dimension, states_dont_care, U, Udag):
    ket_j = getBasis(y[0])
    bra_l = getBasis(y[1]).dag()
    
    Epsilon = evolve_system(ket_k, bra_l, states_dont_care, U, Udag)
    Eps_ijkl = get_eps_ijkl(bra_i, ket_j, Epsilon)
            
    Gstar_ik, G_jl = gate_operators(x, y)
            
    fidStep = fidelity_contribution(gate_dimension, Gstar_ik, G_jl, Eps_ijkl)
    return abs(fidStep[0][0][0])

def increment_fidelity(limits, gate_element_positions, x, bra_i, ket_k, gate_dimension, states_dont_care, U, Udag, queue):
    Fid = 0
    for y in range(limits[0], limits[1]):  
        increment = increment_of_fidelity(x, gate_element_positions[y], bra_i, ket_k, gate_dimension, states_dont_care, U, Udag)
        Fid += increment
    queue.put(Fid)
    return Fid


# Calculates the average fidelity.
def Fidelity (J):  
    gate_element_positions = getGate(G)
    H = hamiltonian(J)
    Fid = 1/(gate_dimension + 1)
    U = (-j*H).expm()
    Udag = U.dag()
    end = time.time()

    start = time.time()

    for x in gate_element_positions: 
        bra_i = getBasis(x[0]).dag()
        ket_k = getBasis(x[1])

        q = multiprocessing.Queue()
        t1 = multiprocessing.Process(target=increment_fidelity, args=([0, int(len(gate_element_positions)/2)], gate_element_positions, x, bra_i, ket_k, gate_dimension, states_dont_care, U, Udag, q))
        t2 = multiprocessing.Process(target=increment_fidelity, args=([int(len(gate_element_positions)/2), len(gate_element_positions)], gate_element_positions, x, bra_i, ket_k, gate_dimension, states_dont_care, U, Udag, q))

        pool = [t1, t2]

        for thread in pool:    
            thread.start()
            response = q.get()
            Fid += response
            
        for thread in pool:
            thread.join()

    print(Fid)
    end = time.time()
    print(end - start)
    return -Fid