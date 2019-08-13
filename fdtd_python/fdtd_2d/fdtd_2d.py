import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from fdtd_helper import (light_speed,
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
x_indices = np.linspace(0, 100, 101, dtype=np.float64)
y_indices = np.linspace(0, 100, 101, dtype=np.float64)
z_indices = np.ndarray((len(x_indices), len(y_indices)), dtype=np.float64)

old_ez_universe = z_indices.copy()
old_hx_universe = z_indices.copy()
old_hy_universe = z_indices.copy()
current_ez_universe = z_indices.copy()
current_hx_universe = z_indices.copy()
current_hy_universe = z_indices.copy()
current_p_universe = z_indices.copy()

permittivity_universe = np.ones((len(x_indices), len(y_indices)), dtype=np.float64)
permeability_universe = np.ones((len(x_indices), len(y_indices)), dtype=np.float64)

e_factor_universe = 0 * z_indices
h_factor_universe = 0 * z_indices

max_frequency = 1.5e6
gausz_width = get_gaussian_tau_from_frequency(max_frequency)
max_fraction = min_fraction = 1
samplerate_space = 20
grid_width_x = get_grid_width(max_frequency, max_fraction, samplerate_space)
grid_width_y = get_grid_width(max_frequency, max_fraction, samplerate_space)
time_step = get_time_resolution_2d(min_fraction, grid_width_x, grid_width_y)
current_time = 0
max_time = 1024 * time_step

image = plt.imshow(z_indices, interpolation="bilinear", animated=True, cmap="jet", vmin=-1, vmax=1)


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
            e_factor_universe[x][y] = get_electrical_update_coefficient_1d(time_step,
                    permittivity_universe[x][y])
            h_factor_universe[x][y] = get_magnetical_update_coefficient_1d(time_step,
                    permeability_universe[x][y])


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
    global current_p_universe
    # Update the magnetic field
    for y in range(1, len(y_indices) - 1):
        for x in range(1, len(x_indices) - 1):
            current_hx_universe[x][y] = old_hx_universe[x][y] - h_factor_universe[x][y] * ((old_ez_universe[x][y + 1] - old_ez_universe[x][y]) / grid_width_y)
            current_hy_universe[x][y] = old_hy_universe[x][y] + h_factor_universe[x][y] * ((old_ez_universe[x + 1][y] - old_ez_universe[x][y]) / grid_width_x)
    handle_perfect_magnetical_conductors()

    # Update the electric field
    for y in range(1, len(y_indices) - 1):
        for x in range(1, len(x_indices) - 1):
            curl_term = ((current_hy_universe[x][y] - current_hy_universe[x - 1][y]) / grid_width_x) - ((current_hx_universe[x][y] - current_hx_universe[x][y - 1]) / grid_width_y)
            current_ez_universe[x][y] = old_ez_universe[x][y] + e_factor_universe[x][y] * curl_term
    handle_perfect_electrical_conductors()

    # Add the source
    middle_x = len(x_indices) // 2
    middle_y = len(y_indices) // 2
    current_ez_universe[middle_x][middle_y] += 50 * gaussian_impulse_1d(current_time, 6 * gausz_width, gausz_width)
    # current_ez_universe[middle_x][middle_y] += 10 * constant_sinus(current_time, 1, max_frequency)
    # current_ez_universe[middle_x][middle_y] += 10 * ramped_sinus(current_time, 50 * time_step, 1, max_frequency)
    # current_ez_universe[middle_x][middle_y] += 10 * rectangular_impulse(current_time, 6 * time_step, 10 * time_step)


    old_ez_universe = current_ez_universe.copy()
    old_hx_universe = current_hx_universe.copy()
    old_hy_universe = current_hy_universe.copy()
    current_time += time_step


def update_image(*args):
    """
    The callback function of an FuncAnimation instance.
    Used to update the fields and yield the image to draw
    as a result
    """
    update_fields()
    image.set_data(current_ez_universe)
    return image,


if __name__ == "__main__":
    set_geometry()
    calculate_factors()
    anim = animation.FuncAnimation(figure, update_image, save_count=100,
                                        interval=5, blit=True)
    plt.show()
