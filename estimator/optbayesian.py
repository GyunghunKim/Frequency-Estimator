# Imports

from .estimator_base import Estimator
import numpy as np
from numpy import pi, cos, exp, sqrt, log


# Class

class OptBayesianEstimator(Estimator):
    def __init__(self, sys, M, SD_B, B0):
        Estimator.__init__(self, sys)
        self.M = M

        self.sd0 = SD_B
        self.sd = SD_B
        self.mu = B0
        self.tl = []
        self.mul = []
        self.rl = []
        self.sdl = []

        self.range = 2  # find maximum B in the range of +/-self.range*self.sd
        self.precision = 10000
        self.state = 0  # runs from 0 to M-1

    def probcoin(self, t, r, f):
        alpha = self.sys.alpha
        beta = self.sys.beta
        return .5 * (1 + r * (alpha + beta * cos(2 * pi * f * t)))

    def prob(self, t, r, mu, sd):
        alpha = self.sys.alpha
        beta = self.sys.beta
        return .5 * (1 + r * alpha + r * beta * cos(2 * pi * mu * t) * exp(
            -2 * (pi * t * sd) ** 2))

    def var(self, t, r, mu, sd):
        alpha = self.sys.alpha
        beta = self.sys.beta
        v = sd ** 2
        num = .5 * v * (1 + r * alpha) - \
              .5 * r * beta * v * cos(2 * pi * t * mu) * (
                          (2 * pi * t * sd) ** 2 - 1) * exp(
            -2 * (pi * t * sd) ** 2)
        return num / self.prob(t, r, mu, sd)

    def update(self):
        if self.state == 15:
            self.sd = self.sd0
            self.state = 0
        if len(self.tl) == 15:
            self.tl.pop(0)
            self.mul.pop(0)
            self.sdl.pop(0)
            self.rl.pop(0)

        k = round(self.mu / self.sd * 0.86038 - .5)
        t = (k + .5) / pi / self.mu
        r = 2 * self.sys.measure(t)

        self.tl.append(t)
        self.mul.append(self.mu)
        self.sdl.append(self.sd)
        self.rl.append(r)

        alpha = self.sys.alpha
        beta = self.sys.beta

        def mle(x):
            res = 0
            for r, t in zip(self.rl, self.tl):
                res += log(1 + r * (alpha + beta * cos(2 * pi * x * t)))
            return res

        f = np.linspace(self.mu - self.range * self.sd,
                        self.mu + self.range * self.sd, self.precision)
        lk = mle(f)

        self.sd, self.mu = sqrt(
            sum([self.probcoin(t, r, self.mu) * self.var(t, r, self.mu, self.sd)
                 for r in [-1, 1]])), f[np.argmax(lk)]

        self.state += 1

    def estimate(self):
        return self.mu