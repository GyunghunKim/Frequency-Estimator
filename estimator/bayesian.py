## IMPORTS

from .estimator_base import Estimator
import numpy as np
import scipy.stats

## CLASSES

class BayesianEstimator(Estimator):
    def __init__(self, sys, M, tau, SD_B, B0):
        """
        Bayesian Estiamtor using MAP
        :param M: Block size used in estimation.
        :param tau: timestep of t_meas = k * tau
        :param SD_B: prior pdf SD value
        """
        Estimator.__init__(self, sys)

        self.M = M
        self.meas_time = 0
        self.tau = tau
        self.previous_estimation = B0
        self.current_estimation = B0
        self.data_list = []
        self.SD_B = SD_B

    def update(self):
        def posterior(B):
            res = 1

            for data in self.data_list:
                alpha = self.sys.alpha
                beta = self.sys.beta
                time = data['meas_time']
                d = data['d']
                if d == 1:
                    res *= 0.5 * (1 - (alpha + beta * np.cos(2 * np.pi * B *
                                    time)))
                elif d == 0:
                    res *= 0.5 * (1 + (alpha + beta * np.cos(2 * np.pi * B *
                                    time)))
            res *= scipy.stats.norm(self.previous_estimation, self.SD_B).pdf(B)

            return res

        self.data_list.append({'meas_time': self.meas_time,
                'd':self.sys.measure(self.meas_time)})
        self.meas_time += self.tau

        # A block of data is completely collected.
        if self.meas_time > self.M * self.tau:
            self.previous_estimation = self.current_estimation

            B_list_precision = 200
            B_list = np.linspace(self.previous_estimation * 0.9,
                    self.previous_estimation * 1.1, B_list_precision)
            self.current_estimation = B_list[np.argmax([posterior(B) for B in
                    B_list])]

            self.meas_time = 0
            self.data_list = []

    def estimate(self):
        return self.current_estimation
