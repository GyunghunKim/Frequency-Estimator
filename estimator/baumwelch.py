## IMPORTS

from .estimator_base import Estimator
import numpy as np
import scipy.stats 

## CLASSES

class BaumWelchFilter(Estimator):
    def __init__(self, sys, tau, SD_B, B0):
        Estimator.__init__(self, sys)
        self.B0 = B0
        self.time = 0
        self.tau = tau

        self.seq_len = 0
        self.seq_len_max = 30
        self.B_count = 100
        self.num_iter = 5
        self.B_list = np.linspace(.5*self.B0, 1.5*self.B0, self.B_count)
        self.SD_B = SD_B
        self.gamma = np.zeros(self.B_count)
        self.data = []

    def update(self):
        # append new measurement data and update the time
        self.seq_len = min(self.seq_len+1, self.seq_len_max)
        self.time += self.tau
        if self.time > self.tau * self.seq_len_max:
            self.time = self.tau

        self.data.append({'time': self.time, 'data':
                self.sys.measure(self.time)})
        if len(self.data) > self.seq_len:
            self.data.pop(0)

        # transition and qubit measurement probabilities
        def transition_prob(i, j):
            return scipy.stats.norm(self.B_list[i],
                    self.SD_B).pdf(self.B_list[j]) * self.B0 / self.B_count

        def qubit_prob(i, j):
            # i is for the measurement results, and j is for B.
            alpha = self.sys.alpha
            beta = self.sys.beta
            if i == 1:
                prob = 0.5 * (1 - (alpha + beta * np.cos(2 * np.pi *
                                self.B_list[j] * self.time)))
            elif i == 0:
                prob = 0.5 * (1 + (alpha + beta * np.cos(2 * np.pi *
                                self.B_list[j] * self.time)))
            
            return prob

        # Initialization
        A = [[transition_prob(i, j) for j in range(self.B_count)] for i in
        range(self.B_count)]
        pi = [transition_prob(i, self.B_count//2) for i in range(self.B_count)]
        B = [[qubit_prob(i, j) for j in range(self.B_count)] for i in range(2)] 

        # E-M iterations
        for _ in range(self.num_iter):
            # E step
            # Forward
            alpha = []
            for i in range(self.B_count):
                alpha.append([pi[i] * B[self.data[0]['data']][i]])
            for t in range(1, self.seq_len):
                for i in range(self.B_count):
                    alpha[i].append(B[self.data[t]['data']][i] *
                            sum([alpha[i][t-1]*A[i][j] 
                                for j in range(self.B_count)]))

            # Backward
            beta = [[] for __ in range(self.B_count)]
            for i in range(self.B_count):
                beta[i].insert(0, 1)
            for t in range(self.seq_len-1, 0, -1):
                for i in range(self.B_count):
                    beta[i].insert(0, sum([beta[j][0] * A[i][j] *
                                B[self.data[t]['data']][j]
                                for j in range(self.B_count)]))

            # M step
            gamma = [[] for __ in range(self.B_count)]
            for i in range(self.B_count):
                gamma[i] = [float(alpha[i][t]*beta[i][t]) / sum([alpha[j][t] *
                        beta[j][t] for j in range(self.B_count)])
                        for t in range(self.seq_len)]
            xi = [[[alpha[i][t] * A[i][j] * beta[j][t+1] *
                B[self.data[t+1]['data']][i] for t in range(self.seq_len-1)]
                for j in range(self.B_count)] for i in range(self.B_count)]

            # Update step
            pi = [gamma[i][0] for i in range(self.B_count)]
            # if self.seq_len != 1:
            #     A = [[float(sum([xi[i][j][t] for t in range(self.seq_len-1)]))/sum([gamma[i][t] for t in
            #             range(self.seq_len-1)]) for j in range(self.B_count)] for i in range(self.B_count)]
            #     B = [[float(sum([gamma[i][t] for t in range(self.seq_len)
            #                 if v == self.data[t]['data']]))/sum([gamma[i][t] for t in
            #                 range(self.seq_len)]) for i in range(self.B_count)] for
            #                 v in range(2)]

        self.gamma = [gamma[i][self.seq_len-1] for i in range(self.B_count)]

    def estimate(self):
        # Estimation
        res = self.B_list[np.argmax(self.gamma)]
        print(f"Estimation: {res}")
        return res
