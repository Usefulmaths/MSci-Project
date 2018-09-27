from differential_evolution import adaptive_differential_evolution
from Fidelity import Fidelity
from qutip import *
from random import randrange

# A global variable that represents the imaginary number, i.
j = complex(0, 1)


class Optimisation:
    '''
    This class provides methods to optimise a quantum network
    to simulate a given ideal quantum gate.
    '''

    def __init__(self, network, ideal_gate, filename):
        self.ideal_gate = ideal_gate
        self.network = network
        self.filename = filename

        # Default random parameters to avoid errors.
        self.params = [float(randrange(-2000, 2000)) /
                       100 for param in range(network.number_of_params())]

        # Fidelity object between the network and ideal gate '''
        self.fidelity_obj = Fidelity(network, ideal_gate)

    def set_params(self, params):
        '''
        Set initial parameters for the network

        Arguments:
                params: the initial parameters of the network.
        '''
        self.params = params

    def optimise_de(self, func, bounds, number_of_processors,
                    iterations=10000000, CR=0.5, F=0.9):
        '''
        Optimise the network via adaptive differential evolution.

        Arguments:
                func: the functions to optimise
                bounds: the parameter bounds to search
                number_of_processors: the number of processes
                                      to run optimisation over

                iterations: the number of iterations to run DE for.
                CR: cross-over rate
                F: mutation factor

        Returns:
                the output of adaptive differential evolution.
        '''
        hard_bounds = [bound for bound in bounds]

        ade = adaptive_differential_evolution(func, iterations,
                                              number_of_processors,
                                              self.network.number_of_params(),
                                              CR, F, bounds,
                                              hard_bounds,
                                              self.filename,
                                              threshold_value=-0.99)
        return ade

    def write_to_file(self, file_name, iters, value, params):
        '''
        Logs the progress of the optimisation over
        each iteration.

        Arguments:
                file_name: the name of the file to log to.
                iters: the current iteration.
                value: the metric value.
                params: the params of the system.
        '''
        file_open = open(file_name + ".txt", "a")
        file_open.write("\n" + str(iters) + "," + str(value))
        file_open.close()
