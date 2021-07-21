## IMPORTS

import numpy as np

## CLASSES

class DataCollector():
    def __init__(self, N):
        """
        :params N: number of estimators
        Notes
         * Data structure: {exp_num:int, exp_data:{time:float, ref:float,
             est:list(float)}}
        """
        self.N = N
        self.data = []

    def append(self, num_exp, time, ref, *est):
        """
        :params num_exp: Experiement number
        :params time: Current time
        :params ref: Reference data (exact one)
        :params *est: Estimations
        """
        self.data.append({'num_exp':num_exp, 'time':time, 'ref':ref, 'est_list':est})
        
    
    def mse(self, num_est):
        """
        :params num_est: MSE of which estimator?
        """
        res = 0
        for item in self.data:
            res += (item['ref']-item['est_list'][num_est])**2

        return np.sqrt(res/len(self.data))
