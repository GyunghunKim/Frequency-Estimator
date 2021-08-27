# IMPORTS
import numpy as np


# CLASSES


class System:
    """
    Class defining the system with qubit and B
    """

    def __init__(self, D, T, alpha, beta, B0):
        self.D = D
        self.time = 0
        self.T = T
        self.alpha = alpha
        self.beta = beta
        self.B = B0

    def update(self):
        self.B += np.random.normal(scale=np.sqrt(2 * self.D * self.T))
        self.time += self.T

    def measure(self, tau):
        prob = 0.5 * (1 - (self.alpha + self.beta * np.cos(2 * np.pi * self.B *
                                                           tau)))
        result = np.random.choice([1, 0], p=[prob, 1 - prob])

        return result
