#include <math.h>

const int light_speed = 299792458;
const double vacuum_permittivity = 8.854187817 * pow(10, -12);
const double vacuum_permeability = 1.2566370614 * pow(10, -6);

double get_lambda_min(double max_frequency, double max_fraction)
{
	return light_speed / (max_frequency * max_fraction);
}

double get_grid_width(double max_frequency, double max_fraction, double samplerate)
{
	return get_lambda_min(max_frequency, max_fraction) / samplerate;
}

double get_time_resolution_2d(double min_fraction, double grid_width_x, double grid_width_y)
{
	double denominator = light_speed * sqrt(pow(grid_width_x, -2) + pow(grid_width_y, -2));
	return min_fraction / denominator;
}

double ramped_sinus(double t, double ramp_len, double amplitude, double frequency)
{
	double ramp_factor = 1.0;
	if (t <= ramp_len) { ramp_factor = t / ramp_len; } 
	return ramp_factor * amplitude * sin(t * 2.0 * M_PI * frequency);
}

