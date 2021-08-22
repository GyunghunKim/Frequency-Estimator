# TODO
# 1. Envelope 찾기. 각 n에 대한 함수에서.
# 2. Envelope를 가정하고, 전략 찾기
# 3. 이것도 엄청난 수준의 근사이므로 그럴듯한 논거가 될 것이다. (fringe 때문에
#    가정 2의 성립이 사실상 불가능함)

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import *

fig = plt.figure()

# Parameters
beta = 0.67
f = 100e6
alpha = 0.25
D = 11.5e3 ** 2 / 1e-6
T = 50e-6


# Plot 1 (Information vs time)
def information0(t):
    return (beta*2*np.pi*t*np.sin(2*np.pi*f*t))**2/(1-(alpha+beta*np.cos(2*np.pi*f*t))**2)

# r belongs to {-1, 1}
def prob(t, n, r, f):
    def integrand(x):
        return 0.5*(1+r*(alpha+beta*np.cos(2*np.pi*x*t)))*1/np.sqrt(4*np.pi*D*T*n)*np.exp(-(x-f)**2/(4*D*T*n))

    n_sig = 10
    return quad(integrand, f-n_sig*np.sqrt(2*D*T*n), f+n_sig*np.sqrt(2*D*T*n))[0]

def dprob(t, n, r, f):
    def integrand(x):
        return 0.5*(1+r*(alpha+beta*np.cos(2*np.pi*x*t)))*1/np.sqrt(4*np.pi*D*T*n)*(x-f)/(2*D*T*n)*np.exp(-(x-f)**2/(4*D*T*n))

    n_sig = 10
    return quad(integrand, f-n_sig*np.sqrt(2*D*T*n), f+n_sig*np.sqrt(2*D*T*n))[0]

def information(t, n=0):

    if n == 0:
        n = 0.001
    # return prob(t, n, 1, f)
    return sum([1/prob(t, n, r, f)*(dprob(t, n, r, f))**2 for r in [-1, 1]])

t_list = np.linspace(1e-9, 100e-9, 100)

ax = fig.add_subplot(2, 2, 1)
ax.set_title("Fisher information\n($f=100\mathrm{MHz}$)")
for n in (0, 50, 100, 500, 1000):
    if n == 0:
        b_list = [information0(t) for t in t_list]
    else:
        b_list = [information(t, n) for t in t_list]
    ax.plot(t_list*1e9, b_list, label=f"$n={n}$")
ax.legend()
ax.set_xlabel("$t$ [ns]")
ax.set_ylabel("$\mathcal{I}(t)$")


# Plot 2
t_list = np.linspace(1e-9, 400e-9, 100)

ax = fig.add_subplot(2, 2, 2)
ax.set_title("Fisher information\n($f=100\mathrm{MHz}$)")
for n in (20, 50, 100, 500, 1000):
    if n == 0:
        b_list = [information0(t) for t in t_list]
    else:
        b_list = [information(t, n) for t in t_list]
    ax.plot(t_list*1e9, b_list, label=f"$n={n}$")
ax.legend()
ax.set_xlabel("$t$ [ns]")
ax.set_ylabel("$\mathcal{I}(t)$")


# Plot 3
t_list = np.linspace(1e-9, 40e-9, 100)

ax = fig.add_subplot(2, 2, 3)
ax.set_title("Fisher information\n($n=0$)")
for f in np.linspace(95e6, 105e6, 3):
    b_list = [information0(t) for t in t_list]
    ax.plot(t_list*1e9, b_list, label=f"$f={f/100e6:.2f}$MHz")
f = 100e6 # Don't forget this line
ax.legend()
ax.set_xlabel("$t$ [ns]")
ax.set_ylabel("$\mathcal{I}(t)$")

# Finding the maximum information
# t_list = np.linspace(1e-9, 400e-9, 300)
t_list = np.linspace(1e-9, 3e-6, 300)

ax = fig.add_subplot(2, 2, 4)
ax.set_title("Fisher information\n($f=100\mathrm{MHz}$)")
for n in (1,2,3,4,5,6,7,8,9,10):
    if n == 0:
        b_list = [information0(t) for t in t_list]
    else:
        b_list = [information(t, n) for t in t_list]
    ax.plot(t_list*1e9, b_list, label=f"$n={n}$")
ax.legend()
ax.set_xlabel("$t$ [ns]")
ax.set_ylabel("$\mathcal{I}(t)$")


# index_range = range(1, 100)
# 
# mse_by_average_list = []
# for i in index_range:
#     mse_by_average_list.append(np.sqrt(2 * T * D * i / 2))
# 
# ax2 = fig.add_subplot(2, 2, 2)
# ax2.set_title("MSE by fluctuation")
# ax2.plot(index_range, mse_by_average_list)
# ax2.set_xlabel("Number of samples")
# ax2.set_ylabel("MSE [Hz]")
# ax2.ticklabel_format(axis="y", style="sci", scilimits=(6, 6), useOffset=False)

# mse_list = []
# mse_using_max_info_list = []
# info_dict = dict(zip(t_list, b_list))
# max_info_dict = dict(zip(t_list, b_max_list))
# for i in index_range:
#     total_info = 0.0
#     total_info_using_max = 0.0
#     for t in t_list[:i+1]:
#         total_info += info_dict[t]
#         total_info_using_max += max_info_dict[t]
#     mse_list.append(1/np.sqrt(total_info))
#     mse_using_max_info_list.append(1/np.sqrt(total_info_using_max))
# 
# ax3 = fig.add_subplot(2, 2, 3)
# ax3.plot(index_range, mse_list)
# ax3.plot(index_range, mse_using_max_info_list)
# ax3.set_xlabel("Number of samples")
# ax3.set_ylabel("MSE [Hz]")
# ax3.ticklabel_format(axis="y", style="sci", scilimits=(6, 6), useOffset=False)

plt.tight_layout()
plt.show()

