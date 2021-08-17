## IMPORTS

from .estimator_base import Estimator
import numpy as np
import scipy.stats
from scipy.optimize import minimize

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

    @staticmethod
    def findmax(func, initial_guess, algorithm='scipy'):
        ans = 0

        if algorithm == 'brute_force':
            B_list_upper_limit = 1.1
            B_list_lower_limit = 0.9
            B_list_precision = 100

            B_list = np.linspace(initial_guess * B_list_upper_limit,
                    initial_guess * B_list_lower_limit, B_list_precision)
            ans = B_list[np.argmax([func(B) for B in B_list])]

        if algorithm == 'scipy':
            ans = minimize(lambda x: -func(x), initial_guess,
                    method='Powell').x[0]

        return ans
            

    def update(self):
        def logposterior(B):
            res = 0

            for data in self.data_list:
                alpha = self.sys.alpha
                beta = self.sys.beta
                time = data['meas_time']
                d = data['d']
                if d == 1:
                    res += np.log(0.5 * (1 - (alpha + beta * np.cos(2 * np.pi *
                                        B * time))))
                elif d == 0:
                    res += np.log(0.5 * (1 + (alpha + beta * np.cos(2 * np.pi *
                                        B * time))))
            res += scipy.stats.norm(self.previous_estimation,
                    self.SD_B).logpdf(B)

            return res

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

            self.current_estimation = BayesianEstimator.findmax(posterior,
                    self.previous_estimation)

            self.meas_time = 0
            self.data_list = []

    def estimate(self):
        return self.current_estimation

class BayesianSmoothEstimator(BayesianEstimator):
    def __init__(self, sys, M, tau, SD_B, B0):
        BayesianEstimator.__init__(self, sys, M ,tau, SD_B, B0)
        self.B0 = B0
        self.tau = tau

    def update(self):
        def logposterior(B):
            res = 0

            for data in self.data_list:
                alpha = self.sys.alpha
                beta = self.sys.beta
                time = data['meas_time']
                d = data['d']
                if d == 1:
                    res += np.log(0.5 * (1 - (alpha + beta * np.cos(2 * np.pi *
                                        B * time))))
                elif d == 0:
                    res += np.log(0.5 * (1 + (alpha + beta * np.cos(2 * np.pi *
                                        B * time))))
            res += scipy.stats.norm(self.previous_estimation,
                    self.SD_B).logpdf(B)

            return res

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

        if len(self.data_list) > self.M:
            self.data_list.pop(0)

        self.previous_estimation = self.current_estimation
        self.current_estimation = BayesianEstimator.findmax(posterior,
                self.previous_estimation)

        if len(self.data_list) < self.M:
            self.current_estimation = self.B0
    
        if self.meas_time > self.tau * self.M:
            self.meas_time = 0

    def estimate(self):
        return self.current_estimation

class BayesianExpSamplingEstimator(BayesianEstimator):
    def __init__(self, sys, M, tau, SD_B, B0, sig):
        BayesianEstimator.__init__(self, sys, M, tau, SD_B, B0)

        self.meas_time = tau
        self.sig = sig

    def update(self):
        def logposterior(B):
            res = 0

            for data in self.data_list:
                alpha = self.sys.alpha
                beta = self.sys.beta
                time = data['meas_time']
                d = data['d']
                if d == 1:
                    res += np.log(0.5 * (1 - (alpha + beta * np.cos(2 * np.pi *
                                        B * time))))
                elif d == 0:
                    res += np.log(0.5 * (1 + (alpha + beta * np.cos(2 * np.pi *
                                        B * time))))
            res += scipy.stats.norm(self.previous_estimation,
                    self.SD_B).logpdf(B)

            return res

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
        self.meas_time *= self.sig

        # A block of data is completely collected.
        if self.meas_time > (self.sig ** self.M) * self.tau:
            self.meas_time = self.tau
            self.data_list.pop(0)

        self.previous_estimation = self.current_estimation
        self.current_estimation = BayesianEstimator.findmax(posterior,
                self.previous_estimation)
        
    def estimate(self):
        return self.current_estimation
