## IMPORTS

import numpy as np
import matplotlib.pyplot as plt

## CLASSES

class DataCollector():
    def __init__(self):
        """
        :params N: number of estimators
        Notes
         * Data structure: {exp_num:int, exp_data:{time:float, ref:float,
             est:list(float)}}
        """
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

    def visualize_exp(self, num_exp, *est):
        time = [item['time'] for item in self.data]
        ref = []
        B_list = []
        for elem in est:
            if elem == 'ref':
                ref = [item['ref'] for item in self.data]
            else:
                B_list.append([item['est_list'][elem] for item in self.data])

        plt.figure(dpi=100)
        plt.title('$\Delta B_{z}$')
        plt.xlabel('time [s]')
        plt.ylabel('$\Delta B_{z}$ [Hz]')
        plt.ticklabel_format(axis='y', style='sci', scilimits=(6, 6), useOffset=False, useMathText=True)
        if ref:
            plt.plot(time, ref, label='True')
        for i in range(len(B_list)):
            plt.plot(time, B_list[i], label='$\Delta B_{}$'.format(i))
        plt.legend()
        plt.show()
