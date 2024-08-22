import pandas as pd
import numpy as np
import glob
import re
import matplotlib.pyplot as plt
from   matplotlib import cm

magnetization_versus_temp_and_theta_dict = {}

total_values_for_mean = 1000
theta_list         = []
temperature_list   = []
magnetization_list = []

for filename in glob.glob("heatmap_data/*.dat"):
    # toda essa parte de dar match nas strings eh horroroso, horrivel
    numbers = re.findall(r'\d+', filename)

    temperature  = (float("{}.{}".format(numbers[0], numbers[1])))
    theta        = int(numbers[2])
    lattice_size = int(numbers[3])
    seed         = int(numbers[-1])

    if lattice_size == 32:
        last_values_for_mean = np.array(pd.DataFrame(pd.read_csv(filename, dtype=str))
            ["#seed = {}".format(seed)].iloc[-total_values_for_mean:])

        mean_magnetization = 0

        # isso aqui tb ta mt estupido
        for value in last_values_for_mean:
            mean_magnetization += float(value)

        mean_magnetization = mean_magnetization / total_values_for_mean

        theta_list.append(theta)
        temperature_list.append(temperature)
        magnetization_list.append(mean_magnetization)

#for i in range(0, 100):
#    print(theta_list[i], temperature_list[i], magnetization_list[i])

#nao entendi pq deu o formato errado, mas isso resolve pra alguns casos..

x = np.unique(temperature_list)
y = np.unique(theta_list)
X, Y = np.meshgrid(x,y)

Z = np.array(magnetization_list).reshape(len(y),len(x))

fig1, ax2 = plt.subplots(layout='constrained')
CS = ax2.contourf(X, Y, Z, 10, cmap=plt.cm.bone)

ax2.set_xlabel('temperature')
ax2.set_ylabel('theta')

cbar = fig1.colorbar(CS)
cbar.ax.set_ylabel('magnetization')

plt.savefig('heatmap_mag_versus_T_and_theta.png', dpi=400)
plt.clf()

'''dict_mag_versus_theta_and_temperature = {}

for index in range(len(theta_list)):
    if theta_list[index] not in dict_mag_versus_theta_and_temperature[theta_list[index]].keys():
        dict_mag_versus_theta_and_temperature[theta_list[index]] = [temperature_list[index], magnetization_list[index]]
    else:
        dict_mag_versus_theta_and_temperature[theta_list[index]] = [temperature_list[index], magnetization_list[index]]


for index in range(len(theta_list)):
    plt.scatter(temperature_list[index], magnetization_list[index], label = str(theta_list[index]))
    plt.legend(str(theta_list[index]))

plt.legend(loc='upper right')

plt.savefig('simple_plot.png', dpi = 400)'''
