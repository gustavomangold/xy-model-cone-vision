import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import io
import glob

def get_correlation(data, d):
    '''
    Calculate the average spin-spin correlation

    Args:
    data (np.array): Configuration Snapshot
    d (int): Max distance

    Returns:
     (float): Average correlation
    '''
    correlation = 0

    for index in range(len(data)):
        angle_of_spin = float(data[index])
        angle_of_next = float(data[(index + d) % len(data)])

        correlation += np.dot(np.array(np.cos(angle_of_spin), np.sin(angle_of_spin)),
                              np.array(np.cos(angle_of_next), np.sin(angle_of_next)))

    return correlation / len(data)

def get_lists_to_plot(address_to_open):
    '''
    Calculate the average correlation for many configurations

    Args:
    address_to_open (str): Path of configurations

    Returns:
    x_to_plot (list): List of distance
    y_to_plot (list): List of Average correlation as a function of distance
    '''
    list_of_angles = []

    with io.open(address_to_open, mode="r", encoding="utf-8") as f:
        for line in f:
            list_of_angles.append(line.split())

    data = (np.float64(np.array(list_of_angles)[0]))

    x_to_plot = []
    y_to_plot = []

    for d in range(12):
        x_to_plot.append(d)
        y_to_plot.append(get_correlation(data, d))

    return x_to_plot, y_to_plot

path = "./data_new_version/corr/"
colormap = ['red', 'blue', 'cyan', 'black', 'green', 'grey', 'hotpink']
markers  = ['1', 'p', 'o', 'v', '^', '>', '*']

angle = '315'

for filename in glob.glob(path + '/finalconfig_T*L32A' + angle + '*.dat'):
    index = 0
    for temp in np.arange(0.5, 1.15, 0.2):
        if 'T' + '{:.3f}'.format(temp) in filename:
            x, y = get_lists_to_plot(filename)
            plt.scatter(x, y, label = r'$T=$' + '{:.1f}'.format(temp), marker = markers[index % len(markers)],
                        color = colormap[index % len(colormap)])
            plt.plot(x, y, alpha = 0.2, color = colormap[index % len(colormap)])
        index += 1

#plt.legend(loc='upper right', fontsize = 'small')
plt.xlabel('d')
#plt.ylabel(r'$\langle S_i S_{i+d} \rangle$')
plt.yscale('log')
plt.xscale('log')
plt.savefig('correlation-' + angle + '.png', dpi=400)
