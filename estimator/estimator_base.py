## IMPORTS

from abc import abstractmethod

## CLASSES

class Estimator():
    def __init__(self, sys):
        self.sys = sys

    @abstractmethod
    def update(self):
        pass

    @abstractmethod
    def estimate(self):
        pass
