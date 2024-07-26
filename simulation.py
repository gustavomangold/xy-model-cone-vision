import numpy as np
import math
import matplotlib.pyplot as plt
import os
from   matplotlib.colors import ListedColormap as lcm

def simulation_for_temperature(temperature):
    def plot_snapshot(matrix, temperature=temperature):
        figure = plt.figure()
        axes = figure.add_subplot(111)

        # using the matshow() function
        cax = axes.matshow(matrix_of_spins, cmap='hot', interpolation='nearest', rasterized=True)
        figure.colorbar(cax)

        path_folder = '/Figures'
        if not os.path.exists(path_folder):
            os.makedirs(path_folder)

        plt.savefig(path_folder + '/initial_matrix-temperature=' + str(temperature) + '.png', dpi=400)

        plt.clf()

    def calculate_magnetization(matrix, length):
        #vector_form = np.vectorize(matrix)

        sum_cos = np.sum((np.cos(matrix)))
        sum_sin = np.sum((np.sin(matrix)))

        return math.sqrt(sum_cos**2 + sum_sin**2)/(length*length)

    def coupling_constant(spin_angle, angle_between_spin_neighbour, theta):
        see_spin_variable = min(2*math.pi - abs(spin_angle - angle_between_spin_neighbour),
            abs(spin_angle - angle_between_spin_neighbour))

        if see_spin_variable <= theta / 2:
            return 1
        else:
            return 0

    def glauber(energy_sampled, energy_new, T):
        return 0.5*(1-math.tanh((energy_sampled-energy_new)/(2*T)))

    length         = 32
    total_mc_steps = 500
    #aprox. igual a 330º, usado no artigo
    theta          = 2*math.pi*(330/360)

    matrix_of_spins = np.zeros((length, length))

    for row in range(length):
        for column in range(length):
            matrix_of_spins[row][column] = np.random.uniform(0, 2*math.pi)

    print(calculate_magnetization(matrix_of_spins, length))

    #sempre verdadeiro pra rede quadrada
    angle_between_spin_and_neighbour_below = (3*math.pi)/2
    angle_between_spin_and_neighbour_left  = math.pi #(2*math.pi)/2
    angle_between_spin_and_neighbour_above = (math.pi)/2
    angle_between_spin_and_neighbour_right = 0 #(0*math.pi)/2

    for step in range(total_mc_steps):
        for sample in range(length**2):
            sampled_row    = np.random.randint(0, length)
            sampled_column = np.random.randint(0, length)

            sampled_spin = matrix_of_spins[sampled_row][sampled_column]

            #vizinho abaixo
            spin_neighbour_below = matrix_of_spins[(sampled_row + 1) % length][sampled_column]
            #vizinho acima
            spin_neighbour_above = matrix_of_spins[sampled_row - 1][sampled_column]
            #vizinho a direita
            spin_neighbour_right = matrix_of_spins[sampled_row][(sampled_column + 1) % length]
            #vizinho a esquerda
            spin_neighbour_left  = matrix_of_spins[sampled_row][sampled_column - 1]

            total_energy_sampled = coupling_constant(sampled_spin, angle_between_spin_and_neighbour_below, theta) * math.cos(sampled_spin - spin_neighbour_below) \
                                 + coupling_constant(sampled_spin, angle_between_spin_and_neighbour_right, theta) * math.cos(sampled_spin - spin_neighbour_right) \
                                 + coupling_constant(sampled_spin, angle_between_spin_and_neighbour_left, theta) * math.cos(sampled_spin - spin_neighbour_left) \
                                 + coupling_constant(sampled_spin, angle_between_spin_and_neighbour_above, theta) * math.cos(sampled_spin - spin_neighbour_above)

            new_spin_candidate = np.random.uniform(0, 2*math.pi)

            total_energy_new    = coupling_constant(new_spin_candidate, angle_between_spin_and_neighbour_below, theta) * math.cos(sampled_spin - spin_neighbour_below) \
                                + coupling_constant(new_spin_candidate, angle_between_spin_and_neighbour_right, theta) * math.cos(sampled_spin - spin_neighbour_right) \
                                + coupling_constant(new_spin_candidate, angle_between_spin_and_neighbour_left, theta) * math.cos(sampled_spin - spin_neighbour_left) \
                                + coupling_constant(new_spin_candidate, angle_between_spin_and_neighbour_above, theta) * math.cos(sampled_spin - spin_neighbour_above)

            #glauber transition rate
            if np.random.uniform(0, 1) < glauber(total_energy_sampled, total_energy_new, temperature):
                matrix_of_spins[sampled_row][sampled_column] = new_spin_candidate

    plot_snapshot(matrix_of_spins)
    print('Simulation for temperature T={:.2f} finished.'.format(temperature))
    return calculate_magnetization(matrix_of_spins, length)

list_of_magnetizations = []
range_of_temps = np.arange(1e-3, 1.4, 0.2)

for temperature in range_of_temps:
    magnetization = simulation_for_temperature(temperature)
    list_of_magnetizations.append(magnetization)
    print(magnetization)

print(range_of_temps, list_of_magnetizations)

plt.scatter(range_of_temps, list_of_magnetizations)
plt.savefig('magnetization.png', dpi=400)
