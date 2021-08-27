# IMPORTS

from .estimator_base import Estimator
import numpy as np
from numpy import sqrt, exp, log, pi, cos, sin
from scipy.optimize import minimize_scalar


# CLASSES

class ContinuousBaumWelchEstimator(Estimator):
    def __init__(self, sys, M, tau, B0, D, T, SD, time='linear'):
        Estimator.__init__(self, sys)
        self.alpha = self.sys.alpha
        self.beta = self.sys.beta
        self.B0 = B0
        self.previousB = B0
        self.DT = D * T
        self.SD = SD

        self.tau = tau
        self.meas_time_index = 0
        if time == 'linear':
            self.meas_times = list(np.linspace(tau, (M+1) * tau, M))
        elif time == 'exponential':
            self.meas_times = np.multiply([1.113 ** i for i in range(1, M+1)], 1e-9)
        self.M = M
        self.count = 0
        self.precision = 10000

        self.data = []

    def update(self):
        self.meas_time_index += 1
        self.meas_time_index %= len(self.meas_times)
        meas_time = self.meas_times[self.meas_time_index]

        self.data.append(
            {'meas_time': meas_time, 'result': self.sys.measure(meas_time)})
        if len(self.data) > self.M:
            self.data.pop(0)

        self.count += 1

    def coin_prob(self, t, f, r):
        if r == 1:
            prob = 0.5 * (1 - (self.alpha + self.beta * cos(2 * pi * f * t)))
        elif r == 0:
            prob = 0.5 * (1 + (self.alpha + self.beta * cos(2 * pi * f * t)))

        return prob

    def posterior(self, datum, n, i):
        t = datum['meas_time']
        # Conventions are different..
        r = 1 - 2 * datum['result']
        sd = sqrt(2 * n * self.DT)

        return .5 * (1 + r * self.alpha + r * self.beta * cos(
            2 * pi * i * t) * exp(-2 * (pi * sd * t) ** 2))

    def log_likelihood(self, x):
        log_prob = 0

        for s, datum in enumerate(self.data):
            log_prob += log(self.posterior(datum, len(self.data) - s - 1, x))

        log_prob += -(x - self.B0) ** 2 / (4 * self.count * self.DT)
        # log_prob += -(x - self.B0) ** 2 / (4 * len(self.data) * self.DT)
        log_prob += -(x - self.previousB) ** 2 / (4 * 10 * self.DT)

        return log_prob

    def dp(self, x):  # Derivative of the above function
        res = 0

        for s, datum in enumerate(self.data):
            n = len(self.data) - s - 1
            sd = sqrt(2 * n * self.DT)
            p = self.posterior(datum, n, x)
            t = datum['meas_time']
            r = 1 - 2 * datum['result']
            res += 1 / p * .5 * r * self.beta * 2 * pi * t * sin(
                2 * pi * x * t) * exp(-2 * (pi * sd * t) ** 2)

        res -= (x - self.B0) / (2 * self.count * self.DT)
        # res -= (x - self.previousB) / (2 * self.DT)

        return res

    def estimate(self):
        # res = minimize_scalar(lambda x: -self.log_likelihood(x), bounds=[.5*self.B0, 1.5*self.B0]).x
        f = np.linspace(self.previousB-20*self.SD, self.previousB+20*self.SD, self.precision)
        res = f[np.argmax(self.log_likelihood(f))]

        self.previousB = res

        return res
