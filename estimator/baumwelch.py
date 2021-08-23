# IMPORTS

from .estimator_base import Estimator
import numpy as np
from numpy import sqrt, exp, log, pi, cos
from scipy import integrate
from scipy.optimize import minimize
from matplotlib import pyplot as plt


# CLASSES

class ContinuousBaumWelchFilter(Estimator):
    def __init__(self, sys, tau, B0, D, T):
        Estimator.__init__(self, sys)
        self.alpha = self.sys.alpha
        self.beta = self.sys.beta
        self.B0 = B0
        self.previousB = B0
        self.DT = D * T

        self.tau = tau
        self.meas_time_index = 0
        self.meas_times = list(np.linspace(tau, 101 * tau, 100))

        self.data = []

    def update(self):
        self.meas_time_index += 1
        self.meas_time_index %= len(self.meas_times)
        meas_time = self.meas_times[self.meas_time_index]

        self.data.append({'meas_time': meas_time,
                          'result': self.sys.measure(meas_time)})
        if len(self.data) > 100:
            self.data.pop(0)

    def coin_prob(self, t, f, r):
        if r == 1:
            prob = 0.5 * (1 - (self.alpha + self.beta * cos(2 * pi * f * t)))
        elif r == 0:
            prob = 0.5 * (1 + (self.alpha + self.beta * cos(2 * pi * f * t)))

        return prob

    def posterior(self, datum, n, i):
        if n == 0:
            n = 0.01

        def integrand(j):
            res = self.coin_prob(datum['meas_time'], j, datum['result']) * 1 / sqrt(4 * pi * self.DT * n) * exp(
                -(i - j) ** 2 / (4 * self.DT * n))
            return res

        integral = integrate.quad(integrand, i - 5 * sqrt(2 * self.DT), i + 5 * sqrt(2 * self.DT))[0]
        return integral

    def log_likelihood(self, x):
        log_prob = 0

        for s, datum in enumerate(self.data):
            log_prob += log(self.posterior(datum, len(self.data) - s - 1, x))

        log_prob += -(x - self.B0) ** 2 / (4 * len(self.data) * self.DT)

        return log_prob

    def estimate(self):
        res = minimize(lambda x: -self.log_likelihood(x), self.previousB,
                       method='Powell').x[0]
        self.previousB = res
        return res