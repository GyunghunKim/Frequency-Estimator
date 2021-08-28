# Imports

from .estimator_base import Estimator
import numpy as np
from numpy import pi, cos, exp, sqrt, sin
from scipy.optimize import minimize
import scipy.stats


# Class

class OptBayesianEstimator(Estimator):
    def __init__(self, sys, D, T, B0):
        Estimator.__init__(self, sys)

        step_sdf = sqrt(2*D*T)
        self.step_sdw = step_sdf * 2 * pi
        self.mu = B0 * 2 * pi
        self.var = (1e6 * 2 * pi) ** 2
        self.mul = 10
        self.t = 0
        self.r = 0

    def update(self):
        # Optimal estimation part
        k = round(self.mu/pi/sqrt(self.var)+.5)
        self.t = (k + .5)*pi/self.mu
        t = self.t
        self.r = 1 - 2 * self.sys.measure(t)
        r = self.r
        alpha = self.sys.alpha
        beta = self.sys.beta

        e = exp(-.5 * self.var * t ** 2)

        new_mu = self.mu - (r * beta * self.var * t * sin(self.mu * t) * e) / (
                1 + r * alpha + r * beta * cos(self.mu * t) * e)
        new_var = self.var - r * beta * self.var ** 2 * t ** 2 * e * (
                (1 + r * alpha) * cos(self.mu * t) + r * beta * e) / (
                          1 + r * alpha + r * beta * cos(self.mu * t) * e) ** 2

        self.mu = new_mu
        self.var = new_var + (self.step_sdw ** 2) * self.mul

    def estimate(self):
        return self.mu/2/pi

class CombinedBayesianEstimator(OptBayesianEstimator):
    def __init__(self, sys, D, T, B0, M, tau, SD_B, threshold):
        OptBayesianEstimator.__init__(self, sys, D, T, B0)

        self.M = M
        self.tau = tau
        self.SD_w = SD_B*2*pi
        self.thresholdf = threshold
        self.bayesian_mu = B0*2*pi
        self.count = 0
        self.rl = []
        self.tl = []

    def update(self):
        OptBayesianEstimator.update(self)

        alpha = self.sys.alpha
        beta = self.sys.beta

        self.tl.append(self.t)
        self.rl.append(self.r)
        self.count += 1 
        if self.count == self.M:
            def posterior(w):
                res = 1
                for t, r in zip(self.tl, self.rl):
                    res *= (1+r*(alpha+beta*cos(w*t)))
                res *= scipy.stats.norm(self.bayesian_mu, self.SD_w).pdf(w)

                return res

            # 평균으로 추정하는 것이 자연스러울 수 있다. Powell을 얼마나 믿을 수 있는가?
            self.bayesian_mu = minimize(lambda x: -posterior(x),
                                        100e6*pi*2, method='Powell').x[0]
            self.count = 0
            self.tl = []
            self.rl = []
        
        #if np.abs(self.mu - self.bayesian_mu) > self.thresholdf*2*pi:
        #    print("Ouch")
        #    self.mu = self.bayesian_mu

    def estimate(self):
        return self.mu/2/pi
        # return self.bayesian_mu/2/pi
