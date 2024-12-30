import numpy as np
from numpy import sqrt, pi, sin, cos, log, log10, exp, tanh, sinh, cosh
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import matplotlib.font_manager as font_manager
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import scipy
import scipy.special as sc
import struct
import re
import os
from mpl_toolkits.axes_grid1 import make_axes_locatable

from multiprocessing import Pool, Queue
# import gc


# Function for loading parameters
def load_params(project_dir):
    with open(project_dir + "param.dat", mode="rb") as f:
        param_raw = f.read()
    with open(project_dir + "paramTypes.txt") as f:
        param_types = f.read().split('\n')[:-1]
    with open(project_dir + "paramNames.txt") as f:
        param_names = f.read().split('\n')[:-1]

    param_type_map = dict()
    param_type_map["Integer64"] = "q"
    param_type_map["Real64"] = "d"

    param_format_str = "".join(list(map(lambda t: param_type_map[t], param_types)))
    param_unpacked = struct.unpack_from(param_format_str, param_raw)
    param = {param_names[i] : param_unpacked[i] for i in range(len(param_unpacked))}

    return param
    
# Load data
project_dir = "/home/hypermania/Research/BHQuasinormalModes/output/outgoing_wavepacket_flat_2/"
param = load_params(project_dir)

h = (param['r_max'] - param['r_min']) / (param['N'] - 1)
r_list = param['r_min'] + h * (np.arange(0, param['N']+1, 1) - 0.5)
t_list = np.fromfile(project_dir + "t_list.dat", dtype=np.float64)



# Font Settings
font_path = font_manager.findfont("Latin Modern Roman")
font = matplotlib.font_manager.FontProperties(fname=font_path)
plt.rcParams.update({
    "text.usetex": True
})


# Plotting
x_bounds = [-50, r_list[-1]]
y_bounds = [-2.5, 2.5]

t_text_pos = [50, 0.7 * y_bounds[1]]

grid_size = len(r_list)


def save_plot(i):
    filename = "state_{}.dat".format(i)
    state = np.fromfile(project_dir + filename, dtype=np.float64)
    psi11 = state[0:grid_size]
    psi22 = state[grid_size:2*grid_size]
    t = t_list[i]

    project_dir_orig = "/home/hypermania/Research/BHQuasinormalModes/output/outgoing_wavepacket_2/"
    state = np.fromfile(project_dir_orig + filename, dtype=np.float64)
    psi22_orig = state[grid_size:2*grid_size]
    
    ax = plt.axes()
    ax.set_xlim(*x_bounds)
    ax.set_ylim(*y_bounds)
    ax.set_xlabel(r'$r_*$',fontsize=15)
    ax.set_ylabel(r'$\Psi_{lm}(t)$',fontsize=15)
    ax.plot(r_list, psi11, linewidth=1, color='tab:blue', label=r'$\Psi_{11}(t)$')
    ax.plot(r_list, 100*psi22, linewidth=1, color='tab:orange', label=r'$10^2\Psi_{22}(t)$')
    ax.plot(r_list, 100*psi22_orig, linewidth=1, color='tab:red', label=r'$10^2\Psi_{22}(t)$')
    ax.text(*t_text_pos,r'$t={:.3f}$'.format(t),fontsize=20,color='0')

    ax.legend()

    plt.grid()
    plt.savefig(project_dir + 'plot_{}.png'.format(i), bbox_inches='tight', dpi=300)
    plt.close()

    plt.clf()

save_plot(100)

with Pool(4) as p:
    p.map(save_plot, range(0, len(t_list)))

# for i in range(0, len(t_list)):
#     save_plot(i)
