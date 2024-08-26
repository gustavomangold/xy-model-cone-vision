import pandas as pd
import numpy as np
import glob
import re
import scipy.interpolate
import matplotlib.pyplot as plt
from   matplotlib import cm

def plot_quantity_versus_temperature(dict_key_is_theta, quantity_id: str):
    if quantity_id == 'binder-cummulant':
        plt.ylim(0.63, 0.67)
        plt.xlim(0.3, 1.0)
        plt.title('Binder Cummulant for L=32')
        plt.xlabel('T')
        plt.ylabel('1-<m^4>/(3<m^2>^2)')
    else:
        plt.xlim(0.25, 1.2)
        plt.ylim(0.05, 1.05)
        plt.title('Magnetization for L=32')
        plt.xlabel('T')
        plt.ylabel('<m>')

    colormap = ['red', 'blue', 'cyan', 'black']
    markers  = ['*', 'p', 'o', 'v']

    index = 0
    for theta in dict_key_is_theta.keys():
        if not theta % 45 and len(dict_key_is_theta[theta]) > 50:
            plt.scatter(*zip(*dict_key_is_theta[theta]), marker = markers[index % len(markers)],
                color = colormap[index % len(colormap)], label = 'Cone angle: ' + str(theta) + 'º')
            #plt.plot(*zip(*dict_key_is_theta[theta]), marker = 'o', color = color, alpha = 0.4)
            plt.legend(fontsize = 'small')
            index += 1

    plt.savefig(quantity_id + '.png', dpi = 400)
    plt.clf()

def plot_heatmap(temperature_list, theta_list, magnetization_list):

    x = np.array(temperature_list)
    y = np.array(theta_list)
    z = np.array(magnetization_list)

    plt.xlim(0.1, 1.15)
    plt.ylim(90, 360)

    plt.yticks(np.arange(90, 365, 45))
    plt.gca().invert_yaxis()

    plt.tricontourf(x, y, z, levels = 50)
    cbar = plt.colorbar()
    cbar.set_ticks(np.arange(0, 1.1, 0.1))

    plt.savefig('heatmap_mag_versus_T_and_theta.png', dpi=400)
    plt.clf()

def get_data_for_heatmap(filename):
    total_values_for_mean = 5000
    mean_magnetization    = 0
    mean_binder_cummulant = 0
    # toda essa parte de dar match nas strings eh horroroso
    numbers = re.findall(r'\d+', filename)

    if 'temp' in filename:
        temperature  = (float("{}.{}".format(numbers[0], numbers[1])))
        theta        = int(numbers[2])
        lattice_size = int(numbers[3])
        seed         = int(numbers[-1])

        if lattice_size == 32:
            dataframe = pd.DataFrame(pd.read_csv(filename, dtype=str, skiprows = 1))

            last_values_for_mean_magnetization = dataframe['M'].iloc[:]
            last_values_for_binder_cummulant   = dataframe['U'].iloc[:]

            mean_magnetization    = np.mean(np.float64(last_values_for_mean_magnetization))
            mean_binder_cummulant = np.mean(np.float64(last_values_for_binder_cummulant))

        return theta, temperature, mean_magnetization, mean_binder_cummulant
    return 0, 0, 0, 0

magnetization_versus_temp_and_theta_dict = {}

theta_list            = []
temperature_list      = []
magnetization_list    = []

dict_to_plot_magnetization = {}
# nome suspeito de variavel
dict_to_plot_binder_cumm   = {}

list_of_filenames = []

# essa parte procura por todos os arquivos na pasta
# a variavel seedless_filename serve pra tirar simulações com mesma seed
# pq eh pra ser uma amostra e tem simulaçao repetida
for filename in glob.glob("data/temp_T0.002Theta\=90L32S231.dat"):
    seedless_filename = filename.split('S', 1)[0]

    if seedless_filename not in list_of_filenames:
        list_of_filenames.append(seedless_filename)

        theta, temperature, mean_magnetization, binder_cummulant = get_data_for_heatmap(filename)
        if (theta != 0):
            if not ((theta in theta_list) and (temperature in temperature_list)
                and (mean_magnetization in magnetization_list)):
                theta_list.append(theta)
                temperature_list.append(temperature)
                magnetization_list.append(mean_magnetization)

            if theta in dict_to_plot_magnetization.keys():
                dict_to_plot_magnetization[theta].append([temperature, mean_magnetization])
                dict_to_plot_binder_cumm[theta].append([temperature, binder_cummulant])
            else:
                dict_to_plot_magnetization[theta] = [[temperature, mean_magnetization]]
                dict_to_plot_binder_cumm[theta]   = [[temperature, binder_cummulant]]

#plot_heatmap(temperature_list, theta_list, magnetization_list)
plot_quantity_versus_temperature(dict_to_plot_magnetization, 'magnetization')
#plot_quantity_versus_temperature(dict_to_plot_binder_cumm, 'binder-cummulant')
