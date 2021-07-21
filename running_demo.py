## IMPORTS

from system import System
from estimator import *
from dataprocessing import DataCollector

## MAIN

if __name__ == "__main__":

    # Simulation parameters
    D = 11.5e3 ** 2 / 1e-6
    T = 50e-6 # Timestep of each experiment
    alpha = 0.25
    beta = 0.67
    B0 = 100e6

    # System definition
    sys = System(D, T, alpha, beta, B0)

    # Estimator parameters
    tau = 1e-9
    M = 20 
    SD_B = 10e6

    # Estimator definition
    estimator = BayesianEstimator(sys, M, tau, SD_B, B0)
    #estimator = MeanFilter(sys, 20, B0)

    # Data collector definition
    dc = DataCollector(1)

    for _ in range(2000):
        # Update the system and estimator
        sys.update()
        estimator.update()

        # Get estimation result
        est = estimator.estimate()
        
        dc.append(0, sys.time, sys.B, est)

    print("MSE({}): {:.03f}MHz".format(type(estimator).__name__, dc.mse(0) / 1e6))
