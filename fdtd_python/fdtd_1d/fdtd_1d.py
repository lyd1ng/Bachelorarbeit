import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from fdtd_helper import (get_grid_width,
                         get_time_resolution_1d_pabc,
                         get_electrical_update_coefficient_1d,
                         get_magnetical_update_coefficient_1d,
                         get_gaussian_tau_from_frequency,
                         rectangular_impulse,
                         gaussian_impulse_1d,
                         ramped_sinus,
                         constant_sinus)


figure, axis = plt.subplots()
x_indices = np.linspace(0, 200, 201, dtype=np.float64)
old_electrical_universe = 0 * x_indices
old_magnetical_universe = 0 * x_indices
current_electrical_universe = 0 * x_indices
current_magnetical_universe = 0 * x_indices
permittivity_universe = np.ones(len(x_indices), dtype=np.float64)
permeability_universe = np.ones(len(x_indices), dtype=np.float64)
e_factor_universe = 0 * x_indices
h_factor_universe = 0 * x_indices
current_time_artist = plt.text(0, 7, '', fontsize=8)
electrical_line, = axis.plot(x_indices, current_electrical_universe)
magnetical_line, = axis.plot(x_indices, current_magnetical_universe)
poynting_line, = axis.plot(x_indices, current_magnetical_universe)


max_frequency = 2.45e9
gausz_width = get_gaussian_tau_from_frequency(max_frequency)
max_fraction = min_fraction = 1
samplerate_space = 20
current_time = 0
grid_width = get_grid_width(max_frequency, max_fraction, samplerate_space)
time_step = get_time_resolution_1d_pabc(max_fraction, grid_width)
last_lower_boundary_value = 0
last_last_lower_boundary_value = 0
last_upper_boundary_value = 0
last_last_upper_boundary_value = 0


def set_geometry():
    """
    Sets the permettivity and permeability of the cells with geometry.
    Free space cells dont have to be set as the arrays are initialised with
    ones
    """
    middle_x = len(x_indices) // 2 + 20
    for x in range(middle_x - 10, middle_x + 10):
        permittivity_universe[x] = 80
    pass


def calculate_factors():
    """
    Uses the permettivity and permeability universe to calculate both
    factor universes
    """
    for x in range(0, len(x_indices)):
        e_factor_universe[x] = get_electrical_update_coefficient_1d(time_step,
                               permittivity_universe[x])
        h_factor_universe[x] = get_magnetical_update_coefficient_1d(time_step,
                               permeability_universe[x])


def handle_perfect_electrical_conductors():
    """
    As perfect electrical conductors can not be expressed by the permettivity
    they have to be treated specially after the e-field is calculated
    by setting the e-field inside a perfect electrical conductor to zero
    """
    for x in range(0, len(x_indices)):
        if permittivity_universe[x] == math.inf:
            current_electrical_universe[x] = 0.0


def handle_perfect_magnetical_conductors():
    """
    As perfect electrical conductors can not be expressed by the permettivity
    they have to be treated specially after the e-field is calculated
    by setting the e-field inside a perfect electrical conductor to zero
    """
    for x in range(0, len(x_indices)):
        if permeability_universe[x] == math.inf:
            current_magnetical_universe[x] = 0.0


def update_fields():
    """
    Update e- and h-field staggered in time and space
    """
    global current_time
    global old_magnetical_universe
    global old_electrical_universe
    global current_magnetical_universe
    global current_electrical_universe
    # Handle the numerical boundary condition at the upper side of the
    # universe using the dirichlet boundary condition
    current_magnetical_universe[-1] = old_magnetical_universe[-1] + h_factor_universe[-1] * ((last_last_upper_boundary_value - old_electrical_universe[-1]) / grid_width)
    # Calculate the magnetic field from the old electrical field
    for x in range(0, len(x_indices) - 1):
        current_magnetical_universe[x] = old_magnetical_universe[x] + h_factor_universe[x] * ((old_electrical_universe[x + 1] - old_electrical_universe[x]) / grid_width)
    abc_upper_boundary()
    handle_perfect_magnetical_conductors()

    # Handle the numerical boundary condition at the lower side of the
    # universe using the dirichlet boundary condition
    current_electrical_universe[0] = old_electrical_universe[0] + e_factor_universe[0] * ((current_magnetical_universe[0] - last_last_lower_boundary_value) / grid_width)
    # Update the electrical field from the CURRENT magnetic field and
    # not the old magnetic field as both fields are staggered in time by
    # half the time step. The simulation will be unstable if you do otherwise
    for x in range(1, len(x_indices)):
        current_electrical_universe[x] = old_electrical_universe[x] + e_factor_universe[x] * ((current_magnetical_universe[x] - current_magnetical_universe[x - 1]) / grid_width)
    abc_lower_boundary()
    # handle_perfect_electrical_conductors()

    # Add the source as a soft electrical source
    # current_electrical_universe[len(x_indices) // 2] += rectangular_impulse(current_time, 60 * time_step, 3 * time_step) * 5
    current_electrical_universe[len(x_indices) // 2] += gaussian_impulse_1d(current_time, 6 * gausz_width, gausz_width)
    # current_electrical_universe[len(x_indices) // 2] += ramped_sinus(current_time, time_step * 100, 1, max_frequency)
    # current_electrical_universe[len(x_indices) // 2] += constant_sinus(current_time, 1, max_frequency)
    # Of course the current fields will be the old field in the next step
    old_magnetical_universe = current_magnetical_universe.copy()
    old_electrical_universe = current_electrical_universe.copy()
    # This is only needed for the simulated source as the numerical solution
    # of the maxwell equations does not include the absolute point in time
    current_time += time_step


def init_electrical_animator():
    """
    Set the lines to zero
    """
    electrical_line.set_ydata([0] * len(x_indices))
    magnetical_line.set_ydata([0] * len(x_indices))
    return electrical_line, magnetical_line


def animate(frame):
    """
    Calculate the e- and h-field and display the fields
    """
    current_time_artist.set_text("current time: " + str(current_time))
    update_fields()
    electrical_line.set_xdata(grid_width * x_indices)
    magnetical_line.set_xdata(grid_width * x_indices)
    electrical_line.set_ydata(current_electrical_universe)
    magnetical_line.set_ydata(current_magnetical_universe)
    poynting_line.set_xdata(grid_width * x_indices)
    poynting_line.set_ydata(current_electrical_universe * current_magnetical_universe)
    return current_time_artist, electrical_line, magnetical_line, poynting_line


def abc_lower_boundary():
    """
    Gets the field values from the lower boundary
    """
    global last_last_lower_boundary_value
    global last_lower_boundary_value
    last_last_lower_boundary_value = last_lower_boundary_value
    last_lower_boundary_value = current_magnetical_universe[0]


def abc_upper_boundary():
    """
    Gets the field values from the upper boundary
    """
    global last_last_upper_boundary_value
    global last_upper_boundary_value
    last_last_upper_boundary_value = last_upper_boundary_value
    last_upper_boundary_value = current_electrical_universe[-1]


if __name__ == "__main__":
    text_box = dict(boxstyle='round', facecolor='grey', alpha=0.5)
    plt.xlabel("z dimension")
    labels = np.linspace(0, len(x_indices), len(x_indices) // 10 + 1) * grid_width
    axis.set_xticks(labels)
    axis.set_xticklabels(labels, rotation=45, ha='right')
    figure.subplots_adjust(bottom=0.2)
    axis.text(0, 9, "time_step: " + str(time_step), fontsize=8, verticalalignment='top', bbox=text_box)
    axis.text(0, 8, "grid_width: " + str(grid_width), fontsize=8, verticalalignment='top', bbox=text_box)

    plt.ylim(-10, 10)
    plt.xlim(0, len(x_indices) * grid_width)
    set_geometry()
    calculate_factors()
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=24, metadata=dict(artist='Me'), bitrate=1800)
    anim = animation.FuncAnimation(figure, animate,
                                   init_func=init_electrical_animator,
                                   interval=1, blit=True,
                                   save_count=1000)
    anim.save('material_blob.mp4', writer=writer)
    # plt.show()
