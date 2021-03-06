# IMPORTS

from system import System
from estimator import *
from dataprocessing import DataCollector

from tqdm import tqdm
import numpy as np
import multiprocessing


# For multiprocessing

def estimator_update(est: Estimator):
    est.update()


def estimator_estimate(est: Estimator):
    return est.estimate()


# MAIN

if __name__ == "__main__":

    # Simulation parameters
    D = 11.5e3 ** 2 / 1e-6
    # D = 0
    T = 50e-6  # Timestep of each experiment
    alpha = 0.25
    beta = 0.67
    B0 = 100e6

    # System definition
    sys = System(D, T, alpha, beta, B0)

    # Estimator parameters
    tau = 1e-9
    M = 100
    SD_B = 50e6
    sig = 1.00693166885
    SD_B_between_steps = np.sqrt(2 * D * T)

    # Estimator definition

    estimators = [ BayesianEstimator(sys, M, tau, SD_B, B0),
        # BayesianSmoothEstimator(sys, M, tau, SD_B, B0),
        # BayesianExpSamplingEstimator(sys, M, tau, SD_B, B0, sig),
        # ContinuousBaumWelchEstimator(sys, 50, tau, B0, D, T, 1e6, time='linear'),
        # ContinuousBaumWelchEstimator(sys, 50, tau, B0, D, T, 1e6, time='exponential'),
        # OptBayesianEstimator(sys, D, T, B0),
        # CombinedBayesianEstimator(sys, D, T, B0, M, tau, SD_B, 10e6),
        # MeanFilter(sys, 50, B0)
    ]

    # Data collector definition
    name_list = [type(estimator).__name__ for estimator in estimators]
    dc = DataCollector(name_list)

    for _ in tqdm(range(300000)):
        # Update the system and estimator
        sys.update()

        for estimator in estimators:
            estimator.update()

        # estimation = pool.map(estimator_estimate, estimators)
        estimation = [estimator.estimate() for estimator in estimators]

        dc.append(0, sys.time, sys.B, estimation)

    # Print the results out
    for i, estimator in enumerate(estimators):
        print("MSE({}): {:.03f}MHz".format(
            type(estimator).__name__, dc.mse(i) / 1e6))

    dc.visualize_exp(0, ['ref', *range(len(estimators))])
