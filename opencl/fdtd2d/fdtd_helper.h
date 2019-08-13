#ifndef FDTD_HELPER
#define FDTD_HELPER

const int light_speed;
const double vacuum_permittivity;
const double vacuum_permeability;

double get_lambda_min(double max_frequency, double max_fraction);
double get_grid_width(double max_frequency, double max_fraction, double samplerate);
double get_time_resolution_2d(double min_fraction, double grid_width_x, double grid_width_y);
double ramped_sinus(double t, double ramp_len, double amplitude, double frequency);

#endif

