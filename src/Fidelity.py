from math import log, sin, cos
from qutip import *

# Define a global variable representing the imaginary number 'i'.
j = complex(0, 1)


class Fidelity:
    '''
    A class used in order to calculate different types of
    fidelity between two quantum gates
    '''

    def __init__(self, network, ideal_gate):
        self.ideal_gate = ideal_gate
        self.network = network

    def get_gate_dimension(self):
        '''
        Returns:
            ideal_gate_shape: the number of dimensions the ideal
            quantum gate has.
        '''
        ideal_gate_shape = self.ideal_gate.shape[0]
        return ideal_gate_shape

    def number_region_qubits(self):
        return int(log(self.get_gate_dimension(), 2))

    def set_generate_dont_care_states(self):
        h = self.network.get_number_qubits() - self.number_region_qubits()
        states_dont_care = [sin(0.847) * basis(2, 1) +
                            cos(0.847) * basis(2, 0)] * (h)
        self.states_dont_care = tensor(states_dont_care)

    def get_gate(self):
        '''
        Returns:
                s: the positions of the non-zero elements in the
                ideal quantum gate.
        '''
        s = []

        rows = self.ideal_gate.shape[0]
        columns = self.ideal_gate.shape[1]

        for i in range(rows):
            for j in range(columns):
                if self.ideal_gate[i][0][j] != 0:
                    s.append([i, j])
        return s

    def set_decomposed_gate(self):
        '''
        Sets the decomposed gate property to the
        positions of the non-zero elements in the
        ideal quantum gate.
        '''
        self.decomposed_gate = self.get_gate()

    def get_bin(self, a):
        '''
        Arguments:
                a: an integer number

        Returns:
                s: the binary representation of s
        '''
        s = bin(a)[2:]
        l = len(s)
        if l < 4:
            r = '0' * (4 - l)
            s = r + s

        return s

    def get_basis(self, a):
        '''
        Arguments:
                a: an integer number

        Returns:
                b: the corresponding quantum basis to a.

        '''
        if a == 0:
            B = [basis(2, 0)] * self.number_region_qubits()
            return tensor(B)

        c = self.get_bin(a)

        b = tensor(basis(2, int(c[-3])),
                   basis(2, int(c[-2])), basis(2, int(c[-1])))

        return b

    def set_basis_array(self):
        '''
        Sets the basis array
        '''
        basis = []
        for i in range(8):
            basis = [self.get_basis(i)
                     for i in range(self.get_gate_dimension())]

        self.basis_array = basis

    def basis_8(self, i):
        a = basis(8, i)
        a.dims = [[2, 2, 2], [1, 1, 1]]
        return a

    def instantiate_xi_and_gv(self):
        '''
        Instantiates xi and gv
        '''
        xi = 0

        for i in range(8):
            xi += 1. / 8**0.5 * \
                tensor(self.basis_8(i), self.basis_8(i), self.states_dont_care)

        gv = 0

        for i in range(8):
            for p in range(8):
                gv += self.ideal_gate.data[i, p] * \
                    tensor(self.basis_8(p), self.basis_8(i))

        self.xi = xi
        self.gv = gv

    def average_gate_fidelity(self, params):
        '''
        Calculates the average gate fidelity, given some
        interaction parameters of the network.

        Arguments:
                params: the interactions parameters of the
                quantum networks.

        Returns:
                fidelity: the average gate fidelity between
                          the quantum gate network and the
                          ideal quantum gate.

                params: the parameters of the network.
        '''
        H = self.network.hamiltonian(params)
        U = (-j * H).expm()
        D = self.ideal_gate.shape[0]

        e = qeye(D)
        e.dims = [[2, 2, 2], [2, 2, 2]]

        xi2 = tensor(e, U) * self.xi

        rho = xi2.ptrace([0, 1, 2, 3, 4, 5])

        overlap = self.gv.dag() * rho * self.gv

        fidelity = -(1. / (D + 1) + 1. / (D + 1) * abs(overlap[0][0][0]))

        return fidelity, params

    def deterministic_gate_fidelity(self, params):
        '''
        Calculates the deterministc gate fidelity, given
        the parameters of the quantum network between the
        ideal quantum gate.

        Arguments:
                params: the parameters of the quantum network

        Returns:
                fidelity: the deterministic gate fidelity
        '''
        H = self.network.hamiltonian(params)
        U = (-j * H).expm()
        Udag = U.dag()
        G = self.ideal_gate

        gate_element_positions = self.decomposed_gate
        gate_dimension = float(self.get_gate_dimension())
        Fid = 1. / (gate_dimension + 1)

        for x in range(len(gate_element_positions)):

            bra_i = self.basis_array[gate_element_positions[x][0]].dag()
            ket_k = self.basis_array[gate_element_positions[x][1]]
            for y in range(len(gate_element_positions)):
                ket_j = self.basis_array[gate_element_positions[y][0]]
                bra_l = self.basis_array[gate_element_positions[y][1]].dag()
                epsilon = U * \
                    tensor(ket_k * bra_l, self.states_dont_care *
                           self.states_dont_care.dag()) * Udag

                eps_ijkl = bra_i * epsilon.ptrace([0, 1, 2]) * ket_j

                g_star_ik = G[gate_element_positions[x][0],
                              gate_element_positions[x][1]].conjugate()
                g_jl = G[gate_element_positions[y][
                    0], gate_element_positions[y][1]]

                fidStep = (1. / (gate_dimension * (gate_dimension + 1))
                           ) * g_star_ik * eps_ijkl * g_jl
                Fid += fidStep[0][0][0]

        fidelity = -abs(Fid)

        return fidelity

    def likelihood(self, params, rho_0):
        '''
        Calculates the likelihood between the quantum network and
        a quantum density, given the network parameters.

        Arguments:
                params: the parameters of the quantum networks
                rho_0: the quantum state density.

        Returns:
                llh: the likelihood between the quantum network
                     and the quantum state density.
        '''
        rho = tensor(rho_0, self.generate_dont_care_states())
        H = self.network.hamiltonian(params)
        U = (-j * H).expm()

        A = (self.ideal_gate * rho_0) * (self.ideal_gate * rho_0).dag()
        B = ((U * rho) * (U * rho).dag()).ptrace([0, 1, 2])

        llh = abs((A * B).tr())

        return llh

    def average_fidelity(self, params, number_iters):
        '''
        Calculates the approximate average fidelity between
        the quantum network and the ideal quantum gate by iterating
        over random quantum states.

        Arguments:
                params: the parameters of the quantum network
                number: the number of quantum states to approximate with

        Returns:
                fidelity: the average fidelity of the quantum network.
        '''
        H = self.network.hamiltonian(params)
        U = (-j * H).expm()
        G = self.ideal_gate

        total_fid = 0

        for iters in range(number_iters):
            psi0 = tensor(rand_ket(N=2), rand_ket(N=2), rand_ket(N=2))
            psi_state = tensor(psi0, self.states_dont_care)

            epsilon = (U * ket2dm(psi_state) * U.dag()).ptrace([0, 1, 2])

            fid_contribution = psi0.dag() * G.dag() * epsilon * G * psi0
            total_fid += fid_contribution[0][0][0]

        fidelity = abs(total_fid) / number_iters

        return fidelity
