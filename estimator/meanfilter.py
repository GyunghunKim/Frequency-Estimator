## IMPORTS

from .estimator_base import Estimator

## CLASSES

class MeanFilter(Estimator):
    def __init__(self, sys, M, B0):
        """
        Mean filter using M samples
        It directly accesses to the system's B
        :param M: number of samples
        """
        Estimator.__init__(self, sys)

        self.M = M
        self.data = [B0] * M 

    def update(self):
        self.data.pop(0)
        self.data.append(self.sys.B)

    def estimate(self):
        return sum(self.data)/self.M
    
