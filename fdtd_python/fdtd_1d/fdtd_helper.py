import math
from numpy import float64

light_speed = 299792458
vacuum_permittivity = 8.854187817e-12
vacuum_permeability = 1.2566370614e-6


def get_lambda_min(max_frequency: float64, max_fraction_index: float64)\
        -> float64:
    """
    Get the minimum occuring wavelength in meters
    """
    return light_speed / (max_frequency * max_fraction_index)


def get_electrical_update_coefficient_1d(time_step: float64,
                                         permittivity: float64) -> float64:
    """
    Gets the electrical update coefficient for a cell
    """
    nominator = light_speed * time_step
    denominator = permittivity
    return nominator / denominator


def get_magnetical_update_coefficient_1d(time_step: float64,
                                         permeability: float64) -> float64:
    """
    Gets the magnetical update coefficients for a cell
    """
    nominator = light_speed * time_step
    denominator = permeability
    return nominator / denominator


def get_grid_width(max_frequency: float64,
                   max_fraction_index: float64,
                   sample_rate: int) -> float64:
    """
    Get the grid resolution in meters which is the minimum wavelength divided
    by the sample rate
    sample_rate=10 => 10 samples in one wavelength
    """
    return get_lambda_min(max_frequency, max_fraction_index) / sample_rate


def get_time_resolution_1d(min_fraction_index: float64,
                           grid_width_x: float64) -> float64:
    """
    Gets the minimum timestep for a 1d universe following
    the courant stabiliy condition.
    THIS IS NOT FOR 1D PERFECT ABSORBING BOUNDARY CONDITION
    """
    nominator = min_fraction_index * grid_width_x
    denominator = light_speed
    return nominator / denominator


def get_time_resolution_1d_pabc(boundary_fraction_index: float64,
                                grid_width_x: float64) -> float64:
    """
    Gets the minimum timestep for a 1d universe NOT following
    the courant stability condtition.
    Instead it is assures that the physical wave will travel one cell
    in exact two timesteps (at the boundaries).
    This way a perfect absorbing boundary condition is fullfilled very
    easily as the force field of the grid outside the universe is
    the force field of the boundary grid two timesteps ago.
    """
    nominator = boundary_fraction_index * grid_width_x
    denominator = 2 * light_speed
    return nominator / denominator


def get_time_resolution_2d(min_fraction_index: float64,
                           grid_width_x: float64,
                           grid_width_y: float64) -> float64:
    """
    Gets the minimum timestep for a 2d universe following
    the courant stability condition
    """
    nominator = min_fraction_index
    denominator = light_speed * \
        math.sqrt(grid_width_x**(-2) + grid_width_y**(-2))
    return nominator / denominator


def get_time_resolution_3d(min_fraction_index: float64,
                           grid_width_x: float64,
                           grid_width_y: float64,
                           grid_width_z: float64) -> float64:
    """
    Gets the minimum timestep for a 3d universe following
    the courant stability condition
    """
    nominator = min_fraction_index
    denominator = light_speed * \
        math.sqrt(grid_width_x**(-1) + grid_width_y**(-2) + grid_width_z**(-2))
    return nominator / denominator


def get_gaussian_tau_from_frequency(max_frequency: float64) -> float64:
    """
    Gets the tau for the gaussian impulse depending on the max
    frequency which should be included in the source
    """
    nominator = 1.0
    denominator = math.pi * max_frequency
    return nominator / denominator


def get_time_resolution_from_gaussian_impulse(tau: float64,
                                              sample_rate: int) -> float64:
    """
    Gets the time resolution depending of the duration of the gaussian
    impulse when used as a source.
    The sample_rate should be at least 10
    """
    return tau / sample_rate


def rectangular_impulse(t: float64, t0: float64, width: float64) -> float64:
    """
    A rectangular impulse
    """
    if t >= t0 - width / 2 and t <= t0 + width / 2:
        return 1.0
    return 0.0


def gaussian_impulse_1d(t: float64, t0: float64, tau: float64) -> float64:
    """
    A gassian impulse in one dimension.
    This could be used as a soft or hard source
    t:   The current time
    t0:  The position of the gauss impulse
    tau: The width of the gauss impulse.
         It specifies the max observed frequency
    """
    nominator = t - t0
    denominator = tau
    return math.exp(-(nominator / denominator)**2)


def ramped_sinus(t: float64,
                 ramp_len: float64,
                 amplitude: float64,
                 frequency: float64) -> float64:
    """
    A ramped sin wave.
    This could be used as a soft or hard soure
    """
    ramp_factor = 1.0
    if t <= ramp_len:
        ramp_factor = t / ramp_len
    return ramp_factor * amplitude * math.sin(t * 2.0 * math.pi * frequency)


def constant_sinus(t: float64,
                   amplitude: float64,
                   frequency: float64) -> float64:
    """
    A constant sin wave.
    This could be used as a soft or hard soure
    """
    return amplitude * math.sin(t * 2.0 * math.pi * frequency)
