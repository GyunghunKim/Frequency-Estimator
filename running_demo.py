## IMPORTS

from system import System
from estimator import *
from dataprocessing import DataCollector

from tqdm import tqdm
import numpy as np

## MAIN

if __name__ == "__main__":

    # Simulation parameters
    D = 11.5e3 ** 2 / 1e-6
    #D = 0
    T = 50e-6 # Timestep of each experiment
    alpha = 0.25
    beta = 0.67
    B0 = 100e6

    # System definition
    sys = System(D, T, alpha, beta, B0)

    # Estimator parameters
    tau = 1e-9
    M = 500
    SD_B = 100e6
    sig = 1.00693166885

    SD_B_between_steps = np.sqrt(2*D*T)

    # Estimator definition
    estimators = [#BayesianEstimator(sys, M, tau, SD_B, B0),
                #BayesianSmoothEstimator(sys, M, tau, SD_B, B0),
                #BayesianExpSamplingEstimator(sys, M, tau, SD_B, B0, sig),
                BaumWelchFilter(sys, tau, SD_B_between_steps, B0)
                #MeanFilter(sys, M, B0)
                ]

    # Data collector definition
    dc = DataCollector()

    for _ in tqdm(range(1000)):
        # Update the system and estimator
        sys.update()

        for estimator in estimators:
            estimator.update()
        
        dc.append(0, sys.time, sys.B,
                *[estimator.estimate() for estimator in estimators])

    # Print the results out
    for i, estimator in zip(range(len(estimators)), estimators):
        print("MSE({}): {:.03f}MHz".format(type(estimator).__name__, dc.mse(i) / 1e6))

    dc.visualize_exp(0, 'ref', *range(len(estimators)))
