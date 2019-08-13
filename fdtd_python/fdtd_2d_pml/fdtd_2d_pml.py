import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from timeit import default_timer as timer
from fdtd_helper import (light_speed,
                         vacuum_permittivity,
                         vacuum_permeability,
                         get_grid_width,
                         get_time_resolution_2d,
                         get_electrical_update_coefficient_1d,
                         get_magnetical_update_coefficient_1d,
                         get_gaussian_tau_from_frequency,
                         rectangular_impulse,
                         gaussian_impulse_1d,
                         ramped_sinus,
                         constant_sinus)


figure, axis = plt.subplots()
x_indices = np.linspace(0, 99, 100, dtype=np.float64)
y_indices = np.linspace(0, 99, 100, dtype=np.float64)
z_indices = np.ndarray((len(x_indices), len(y_indices)), dtype=np.float64)

old_ez_universe = z_indices.copy()
old_hx_universe = z_indices.copy()
old_hy_universe = z_indices.copy()
current_ez_universe = z_indices.copy()
ez_integral_universe = z_indices.copy()
current_hx_universe = z_indices.copy()
hx_integral_universe = z_indices.copy()
current_hy_universe = z_indices.copy()
hy_integral_universe = z_indices.copy()

sigmax_universe = np.zeros((len(x_indices), len(y_indices)), dtype=np.float64)
sigmay_universe = np.zeros((len(x_indices), len(y_indices)), dtype=np.float64)
permittivity_universe = np.ones((len(x_indices), len(y_indices)), dtype=np.float64)
permeability_universe = np.ones((len(x_indices), len(y_indices)), dtype=np.float64)

e_factor1_universe = 0 * z_indices
e_factor2_universe = 0 * z_indices
e_factor3_universe = 0 * z_indices
hx_factor1_universe = 0 * z_indices
hx_factor2_universe = 0 * z_indices
hx_factor3_universe = 0 * z_indices
hy_factor1_universe = 0 * z_indices
hy_factor2_universe = 0 * z_indices
hy_factor3_universe = 0 * z_indices


max_frequency = 1.5e6
gausz_width = get_gaussian_tau_from_frequency(max_frequency)
max_fraction = min_fraction = 1
samplerate_time = 10
samplerate_space = 40
grid_width_x = get_grid_width(max_frequency, max_fraction, samplerate_space)
grid_width_y = get_grid_width(max_frequency, max_fraction, samplerate_space)
time_step = get_time_resolution_2d(min_fraction, grid_width_x, grid_width_y) / samplerate_time
current_time = 0
max_time = 100 * time_step

current_time_artist = plt.text(0, 7, '', fontsize=8)
image = plt.imshow(z_indices, interpolation="bilinear", animated=True, cmap="jet", vmin=-1, vmax=1)
# image = plt.imshow(z_indices, interpolation="bilinear", animated=True, cmap="jet", vmin=0, vmax=10)


def set_pml_geometry():
    """
    """
    pml_width = 20
    pml_height = 20
    f = vacuum_permittivity / (2.0 * time_step)
    for y in range(1, len(y_indices) - 1):
        for x in range(1, pml_width):
            sigmax_universe[pml_width - x - 1][y] = f * (x / pml_width)**3
            sigmax_universe[len(x_indices) - pml_width + x][y] = f * (x / pml_width) ** 3
    for x in range(1, len(x_indices) - 1):
        for y in range(1, pml_height):
            sigmay_universe[x][pml_height - y - 1] = f * (y / pml_height)**3
            sigmay_universe[x][len(y_indices) - pml_height + y] = f * (y / pml_height) ** 3


def set_geometry():
    """
    Sets the permettivity and permeability of the cells with geometry.
    Free space cells dont have to be set as the arrays are initialised with
    ones
    """
    pass


def calculate_factors():
    """
    Uses the permettivity and permeability universe to calculate both
    factor universes
    """
    for x in range(0, len(x_indices)):
        for y in range(0, len(y_indices)):
            e_factor1_universe[x][y] = 1 - ((time_step * (sigmax_universe[x][y] + sigmay_universe[x][y])) / vacuum_permittivity)
            e_factor2_universe[x][y] = (sigmax_universe[x][y] * sigmay_universe[x][y] * (time_step**2)) / (vacuum_permittivity**2)
            e_factor3_universe[x][y] = (light_speed * time_step) / permittivity_universe[x][y]

            hx_factor1_universe[x][y] = (1 / time_step - (sigmay_universe[x][y] / (2 * vacuum_permittivity))) / (1 / time_step + (sigmay_universe[x][y] / (2 * vacuum_permittivity)))
            hx_factor2_universe[x][y] = (light_speed / permittivity_universe[x][y]) / ((1 / time_step) + (sigmay_universe[x][y] / (2 * vacuum_permittivity)))
            nominator = (light_speed * sigmax_universe[x][y] * time_step) / (permeability_universe[x][y] * vacuum_permittivity)
            denominator = 1 / time_step + (sigmay_universe[x][y] / (2 * vacuum_permittivity))
            hx_factor3_universe[x][y] = nominator / denominator

            hy_factor1_universe[x][y] = (1 / time_step - (sigmax_universe[x][y] / (2 * vacuum_permittivity))) / (1 / time_step + (sigmax_universe[x][y] / (2 * vacuum_permittivity)))
            hy_factor2_universe[x][y] = (light_speed / permittivity_universe[x][y]) / ((1 / time_step) + (sigmax_universe[x][y] / (2 * vacuum_permittivity)))
            nominator = (light_speed * sigmay_universe[x][y] * time_step) / (permeability_universe[x][y] * vacuum_permittivity)
            denominator = 1 / time_step + (sigmax_universe[x][y] / (2 * vacuum_permittivity))
            hy_factor3_universe[x][y] = nominator / denominator


def handle_perfect_electrical_conductors():
    """
    As perfect electrical conductors can not be expressed by the permettivity
    they have to be treated specially after the e-field is calculated
    by setting the e-field inside a perfect electrical conductor to zero
    """
    for x in range(0, len(x_indices)):
        for y in range(0, len(y_indices)):
            if permittivity_universe[x][y] == math.inf:
                current_ez_universe[x][y] = 0.0


def handle_perfect_magnetical_conductors():
    """
    As perfect electrical conductors can not be expressed by the permettivity
    they have to be treated specially after the e-field is calculated
    by setting the e-field inside a perfect electrical conductor to zero
    """
    for x in range(0, len(x_indices)):
        for y in range(0, len(y_indices)):
            if permeability_universe[x][y] == math.inf:
                current_hx_universe[x][y] = 0.0
                current_hy_universe[x][y] = 0.0


def update_fields():
    global current_time
    global old_hx_universe
    global old_hy_universe
    global current_hx_universe
    global current_hy_universe
    global old_ez_universe
    global current_ez_universe

    if current_time > max_time:
        exit(0)

    start_time = timer()

    # Update the magnetic field
    for y in range(1, len(y_indices) - 1):
        for x in range(1, len(x_indices) - 1):
            hx_curl_term = ((old_ez_universe[x][y + 1] - old_ez_universe[x][y]) / grid_width_y)
            hy_curl_term = ((old_ez_universe[x + 1][y] - old_ez_universe[x][y]) / grid_width_x)
            hx_integral_universe[x][y] += hx_curl_term
            hy_integral_universe[x][y] += hy_curl_term
            current_hx_universe[x][y] = hx_factor1_universe[x][y] * old_hx_universe[x][y] \
                - hx_factor2_universe[x][y] * hx_curl_term - hx_factor3_universe[x][y] * hx_integral_universe[x][y]
            current_hy_universe[x][y] = hy_factor1_universe[x][y] * old_hy_universe[x][y] \
                + hy_factor2_universe[x][y] * hy_curl_term + hy_factor3_universe[x][y] * hy_integral_universe[x][y]
    handle_perfect_magnetical_conductors()

    # Update the electric field
    for y in range(1, len(y_indices) - 1):
        for x in range(1, len(x_indices) - 1):
            ez_integral_universe[x][y] += e_factor2_universe[x][y] * old_ez_universe[x][y]
            curl_term = ((current_hy_universe[x][y] - current_hy_universe[x - 1][y]) / grid_width_x) - ((current_hx_universe[x][y] - current_hx_universe[x][y - 1]) / grid_width_y)
            # current_ez_universe[x][y] = e_factor1_universe[x][y] * old_ez_universe[x][y] - e_factor2_universe[x][y] * ez_integral_universe[x][y] + e_factor3_universe[x][y] * curl_term
            current_ez_universe[x][y] = e_factor1_universe[x][y] * old_ez_universe[x][y] - ez_integral_universe[x][y] + e_factor3_universe[x][y] * curl_term
    handle_perfect_electrical_conductors()

    # Add the source

    middle_x = len(x_indices) // 2
    middle_y = len(y_indices) // 2
    # current_ez_universe[middle_x][middle_y] += 5 * gaussian_impulse_1d(current_time, 3 * gausz_width, gausz_width)
    # current_ez_universe[middle_x][middle_y] += 10 * constant_sinus(current_time, 1, max_frequency)
    current_ez_universe[middle_x][middle_y] += 10 * ramped_sinus(current_time, 50 * time_step, 1, max_frequency)
    # current_ez_universe[middle_x][middle_y] += 10 * rectangular_impulse(current_time, 6 * time_step, 10 * time_step)

    old_ez_universe = current_ez_universe.copy()
    old_hx_universe = current_hx_universe.copy()
    old_hy_universe = current_hy_universe.copy()
    current_time += time_step
    end_time = timer()

    print(end_time - start_time)


def update_image(*args):
    """
    The callback function of an FuncAnimation instance.
    Used to update the fields and yield the image to draw
    as a result
    """
    update_fields()
    current_time_artist.set_text("current time: " + str(current_time))
    image.set_data(current_ez_universe)
    # image.set_data(permittivity_universe)
    return current_time_artist, image


if __name__ == "__main__":
    set_geometry()
    set_pml_geometry()
    calculate_factors()
    anim = animation.FuncAnimation(figure, update_image, save_count=100,
            interval=1, blit=True)
    plt.show()
