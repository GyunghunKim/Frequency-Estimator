# IMPORTS

from .estimator_base import Estimator
import numpy as np
import scipy.stats
from scipy import integrate
from scipy.optimize import minimize
from matplotlib import pyplot as plt

# CLASSES


class ContinuousBaumWelchFilter(Estimator):
    def __init__(self, sys, tau, SD_B, B0):
        Estimator.__init__(self, sys)
        self.B0 = B0
        self.previousB = B0
        self.time = 0
        self.tau = tau

        self.T = -1  # Number of samples used for estimation
        self.seq_len = 100
        self.SD_B = SD_B
        self.data = []

    def update(self):
        # append new measurement data and update the time
        self.time += self.tau
        if self.time > self.tau * self.seq_len:
            self.time = self.tau

        self.data.append({'time': self.time,
                          'meas': self.sys.measure(self.time)})
        self.T += 1

    def integrate(self, func):
        return integrate.fixed_quad(func, -self.B0, 3 * self.B0)[0]

    def A(self, B1, B2):
        return scipy.stats.norm(B1, self.SD_B).pdf(B2)

    def B(self, B, data):
        alpha = self.sys.alpha
        beta = self.sys.beta
        if data['meas'] == 1:
            prob = 0.5 * (1 - (alpha + beta * np.cos(2 * np.pi * B *
                                                     data['time'])))
        elif data['meas'] == 0:
            prob = 0.5 * (1 + (alpha + beta * np.cos(2 * np.pi * B *
                                                     data['time'])))
        else:
            print("Error: BaumWelchFilter.B")
            return 0

        return prob

    def pi(self, B):
        return self.A(self.B0, B)

    def alpha(self, B, t):
        if t == 0:
            return self.pi(B) * self.B(B, self.data[0])
        else:
            return self.B(B, self.data[t]) * self.integrate(lambda x:
                                                            self.alpha(x, t-1) * self.A(B, x))

    def beta(self, B, t):
        return 1

    def gamma(self, B):
        num = self.alpha(B, self.T) * self.beta(B, self.T)
        # den = self.integrate(lambda x: self.alpha(x, self.T) *
        #                                self.beta(x, self.T))
        den = 1

        return num / den

    def estimate(self):
        # Estimation

        # x = np.linspace(self.B0*0.95, self.B0*1.05, 100)
        # plt.plot(x, [self.gamma(_x) for _x in x])
        # plt.show()
        res = minimize(lambda x: -self.gamma(x), self.previousB,
                       method='Powell').x[0]
        self.previousB = res

        # print(f"Estimation: {res}, T: {self.T}")
        return res

# TODO: Implement a discrete version of Baum-Welch filter


class BaumWelchFilter(Estimator):
    def __init__(self, sys, tau, SD_B, B0):
        Estimator.__init__(self, sys)
        self.B0 = B0
        self.previousB = B0
        self.time = 0
        self.tau = tau

        self.T = -1  # Number of samples used for estimation
        self.seq_len = 100
        self.SD_B = SD_B
        self.data = []

    def update(self):
        # append new measurement data and update the time
        self.time += self.tau
        if self.time > self.tau * self.seq_len:
            self.time = self.tau

        self.data.append({'time': self.time,
                          'meas': self.sys.measure(self.time)})
        self.T += 1

    def integrate(self, func):
        return integrate.fixed_quad(func, -self.B0, 3 * self.B0)[0]

    def A(self, B1, B2):
        return scipy.stats.norm(B1, self.SD_B).pdf(B2)

    def B(self, B, data):
        alpha = self.sys.alpha
        beta = self.sys.beta
        if data['meas'] == 1:
            prob = 0.5 * (1 - (alpha + beta * np.cos(2 * np.pi * B *
                                                     data['time'])))
        elif data['meas'] == 0:
            prob = 0.5 * (1 + (alpha + beta * np.cos(2 * np.pi * B *
                                                     data['time'])))
        else:
            print("Error: BaumWelchFilter.B")
            return 0

        return prob

    def pi(self, B):
        return self.A(self.B0, B)

    def alpha(self, B, t):
        if t == 0:
            return self.pi(B) * self.B(B, self.data[0])
        else:
            return self.B(B, self.data[t]) * self.integrate(lambda x:
                                                            self.alpha(x, t-1) * self.A(B, x))

    def beta(self, B, t):
        return 1

    def gamma(self, B):
        num = self.alpha(B, self.T) * self.beta(B, self.T)
        # den = self.integrate(lambda x: self.alpha(x, self.T) *
        #                                self.beta(x, self.T))
        den = 1

        return num / den

    def estimate(self):
        # Estimation

        # x = np.linspace(self.B0*0.95, self.B0*1.05, 100)
        # plt.plot(x, [self.gamma(_x) for _x in x])
        # plt.show()
        res = minimize(lambda x: -self.gamma(x), self.previousB,
                       method='Powell').x[0]
        self.previousB = res

        # print(f"Estimation: {res}, T: {self.T}")
        return res
