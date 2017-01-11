from qutip import *
from hamiltonian_functions import *
from network_functions import *
from math import log
from math import sin
from math import cos

j = complex(0, 1)
G = quantum_half_adder()
gate_dimension = float(G.shape[0])

# Dimension of states we care about (Region qubits)
dimension_state_care = int(log(gate_dimension, 2))

# Generates a tensor state and identity for the states we don't care about (ancillae qubits)
def generate_dont_care_states(N, dimension_state_care):
    h = N - dimension_state_care
    
    states_dont_care = [sin(0.847) * basis(2,1) + cos(0.847) * basis(2,0)] * (h)
    identity_dont_care = [qeye(2)] * (h)
    return tensor(states_dont_care), tensor(identity_dont_care)

# Extracts non zero element positions of the gate
def getGate (G): #extract non zero elements of the gate and save them in s
    
    s = []

    rows = G.shape[0]
    colums = G.shape[1]
    
    for i in range(rows):
        for j in range(colums):
            if G[i][0][j] != 0 :
                s.append([i,j])
    return s

def get_bin(a):  #binary representation of a number a: useful to write the computational basis 
    
    s = bin(a)[2:]
    l = len(s)
    if l < 4 :
        r = '0'*(4-l)
        s = r + s
        
    return s

def getBasis (a) : #get the basis states according to the binary of a: 10 -> |10>
    
    if a == 0:
          B = [basis(2,0)]*(dimension_state_care)
          return tensor(B)
    
    c = get_bin(a)
    if gate_dimension == 16:
        return tensor(basis(2,int(c[0])) , basis(2,int(c[1])), basis(2,int(c[2])), basis(2,int(c[3])))
    if gate_dimension == 8:
        return tensor(basis(2,int(c[0])) , basis(2,int(c[1])), basis(2,int(c[2])))
    if gate_dimension == 4:
        return tensor(basis(2,int(c[1])) , basis(2,int(c[2])))

states_dont_care, identity_dont_care = generate_dont_care_states(N, dimension_state_care)

# Calculates the average fidelity.
def Fidelity (J):  
    H = hamiltonian(J)
    U = (-j*H).expm()
    Udag = U.dag()

    gate_element_positions = getGate(G)
    Fid = 1.0/(gate_dimension + 1)

    
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
            
            fidStep = (1.0/(gate_dimension*(gate_dimension+1)))*Gstar_ik*Eps_ijkl*G_jl     
            Fid += fidStep[0][0][0]
            
    return -abs(Fid)

def average_fidelity(J):
    H = hamiltonian(J)
    U = (-j * H).expm()
    Udag = U.dag()

    total_fid = 0

    for iters in range(10):
        psi0 = tensor(rand_ket(N = 2), rand_ket(N = 2), rand_ket(N = 2))
        psi_state = tensor(psi0, rand_ket(N = 2))
    
        epsilon = (U * ket2dm(psi_state) * U.dag()).ptrace([0, 1, 2])

        fid_contribution = psi0.dag() * G.dag() * epsilon * G * psi0
        total_fid += fid_contribution
    
    return total_fid / 10
        


