#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "fdtd_helper.h"

#define UWIDTH 100
#define UHEIGHT 100

double current_ez_universe[UWIDTH][UHEIGHT];
double current_hx_universe[UWIDTH][UHEIGHT];
double current_hy_universe[UWIDTH][UHEIGHT];
double old_ez_universe[UWIDTH][UHEIGHT];
double old_hx_universe[UWIDTH][UHEIGHT];
double old_hy_universe[UWIDTH][UHEIGHT];
double current_ez_universe[UWIDTH][UHEIGHT]; 
double current_hx_universe[UWIDTH][UHEIGHT];
double current_hy_universe[UWIDTH][UHEIGHT];
double hx_integral_universe[UWIDTH][UHEIGHT]; 
double ez_integral_universe[UWIDTH][UHEIGHT];
double hy_integral_universe[UWIDTH][UHEIGHT];
double sigmax_universe[UWIDTH][UHEIGHT];
double sigmay_universe[UWIDTH][UHEIGHT];
double permittivity_universe[UWIDTH][UHEIGHT];
double permeability_universe[UWIDTH][UHEIGHT];

double e_factor1_universe[UWIDTH][UHEIGHT];             
double e_factor2_universe[UWIDTH][UHEIGHT];             
double e_factor3_universe[UWIDTH][UHEIGHT];             
double hx_factor1_universe[UWIDTH][UHEIGHT]; 
double hx_factor2_universe[UWIDTH][UHEIGHT];
double hx_factor3_universe[UWIDTH][UHEIGHT];
double hy_factor1_universe[UWIDTH][UHEIGHT];
double hy_factor2_universe[UWIDTH][UHEIGHT];
double hy_factor3_universe[UWIDTH][UHEIGHT];

double max_frequency;
double fraction;
int samplerate_space;
int samplerate_time;
double grid_width_x;
double grid_width_y;
double time_step;
double current_time;
double max_time;

void init_params()
{
	max_frequency = 1.5 * pow(10, 6);
	fraction = 1;
	samplerate_space = 40;
	samplerate_time = 10;
	grid_width_x = get_grid_width(max_frequency, fraction, samplerate_space);
	grid_width_y = get_grid_width(max_frequency, fraction, samplerate_space);
	time_step = get_time_resolution_2d(fraction, grid_width_x, grid_width_y) / (double)samplerate_time;
	current_time = 0;
	max_time = 100 * time_step;
}


void set_pml_geometry(const int uwidth, const int uheight, float* sigmax, float* sigmay)
{
    int pml_width = 20;
    int pml_height = 20;
    double f = vacuum_permittivity / (2.0 * time_step);
    for (int y=1; y < UHEIGHT - 1; y++)
    {
        for (int x=1; x < pml_width; x++)
	{
            sigmax_universe[pml_width - x - 1][y] = f * pow(x / pml_width, 3);
            sigmax_universe[UWIDTH - pml_width + x][y] = f * pow(x / pml_width, 3);
	}
    }
    for (int x=1; x < UWIDTH - 1; x++)
    {
        for (int y=1; y < pml_height; y++)
	{
            sigmay_universe[x][pml_height - y - 1] = f * pow(y / pml_height, 3);
            sigmay_universe[x][UHEIGHT - pml_height + y] = f * pow(y / pml_height, 3);
	}
    }
}

void calculate_factors()
{
    for (int x=0; x < UWIDTH; x++)
    {
        for (int y=0; y < UHEIGHT; y++)
	{
		double nominator = 0;
		double denominator = 0;
            	e_factor1_universe[x][y] = 1 - ((time_step * (sigmax_universe[x][y] + sigmay_universe[x][y])) / vacuum_permittivity);
            	e_factor2_universe[x][y] = (sigmax_universe[x][y] * sigmay_universe[x][y] * pow(time_step, 2)) / pow(vacuum_permittivity, 2);
            	e_factor3_universe[x][y] = (light_speed * time_step) / permittivity_universe[x][y];

            	hx_factor1_universe[x][y] = (1 / time_step - (sigmay_universe[x][y] / (2 * vacuum_permittivity))) / (1 / time_step + (sigmay_universe[x][y] / (2 * vacuum_permittivity)));
            	hx_factor2_universe[x][y] = (light_speed / permittivity_universe[x][y]) / ((1 / time_step) + (sigmay_universe[x][y] / (2 * vacuum_permittivity)));
            	nominator = (light_speed * sigmax_universe[x][y] * time_step) / (permeability_universe[x][y] * vacuum_permittivity);
            	denominator = 1 / time_step + (sigmay_universe[x][y] / (2 * vacuum_permittivity));
            	hx_factor3_universe[x][y] = nominator / denominator;

            	hy_factor1_universe[x][y] = (1 / time_step - (sigmax_universe[x][y] / (2 * vacuum_permittivity))) / (1 / time_step + (sigmax_universe[x][y] / (2 * vacuum_permittivity)));
            	hy_factor2_universe[x][y] = (light_speed / permittivity_universe[x][y]) / ((1 / time_step) + (sigmax_universe[x][y] / (2 * vacuum_permittivity)));
            	nominator = (light_speed * sigmay_universe[x][y] * time_step) / (permeability_universe[x][y] * vacuum_permittivity);
            	denominator = 1 / time_step + (sigmax_universe[x][y] / (2 * vacuum_permittivity));
            	hy_factor3_universe[x][y] = nominator / denominator;
	}
    }
}


void handle_perfect_electrical_conductors()
{
	for (int x=0; x < UWIDTH; x++)
	{
		for (int y=0; y < UHEIGHT; y++)
		{
			if (permittivity_universe[x][y] == 12345)
			{
				current_ez_universe[x][y] = 0.0;
			}
		}
	}
}

void handle_perfect_magnetical_conductors()
{
	for (int x=0; x < UWIDTH; x++)
	{
		for (int y=0; y < UHEIGHT; y++)
		{
			if (permeability_universe[x][y] == 12345)
			{
				current_hx_universe[x][y] = 0.0;
				current_hy_universe[x][y] = 0.0;
			}
		}
	}
}

void update_fields()
{
	if (current_time > max_time) { exit(-1); }

	clock_t start = clock();

	// Update the magnetic field
	for (int y=1; y < UHEIGHT - 1; y++)
	{
		for (int x=1; x < UWIDTH - 1; x++)
		{
			double hx_curl_term = ((old_ez_universe[x][y + 1] - old_ez_universe[x][y]) / grid_width_y);
			double hy_curl_term = ((old_ez_universe[x + 1][y] - old_ez_universe[x][y]) / grid_width_x);
			hx_integral_universe[x][y] += hx_curl_term;
			hy_integral_universe[x][y] += hy_curl_term;
			current_hx_universe[x][y] = hx_factor1_universe[x][y] * old_hx_universe[x][y] - hx_factor2_universe[x][y] * hx_curl_term - hx_factor3_universe[x][y] * hx_integral_universe[x][y];
			current_hy_universe[x][y] = hy_factor1_universe[x][y] * old_hy_universe[x][y] + hy_factor2_universe[x][y] * hy_curl_term + hy_factor3_universe[x][y] * hy_integral_universe[x][y];
		}
	}
	handle_perfect_magnetical_conductors();


	// Update the electrical field
	for (int y=1; y < UHEIGHT - 1; y++)
	{
		for (int x=1; x < UWIDTH - 1; x++)
		{
			ez_integral_universe[x][y] += e_factor2_universe[x][y] * old_ez_universe[x][y];
			double curl_term = ((current_hy_universe[x][y] - current_hy_universe[x - 1][y]) / grid_width_x) - ((current_hx_universe[x][y] - current_hx_universe[x][y - 1]) / grid_width_y);
            		current_ez_universe[x][y] = e_factor1_universe[x][y] * old_ez_universe[x][y] - ez_integral_universe[x][y] + e_factor3_universe[x][y] * curl_term;
		}
	}
	handle_perfect_electrical_conductors();



	int middle_x = (int)UWIDTH / 2;
	int middle_y = (int)UHEIGHT / 2;
	current_ez_universe[middle_x][middle_y] += 10 * ramped_sinus(current_time, 50 * time_step, 1, max_frequency);

	memcpy(old_ez_universe, current_ez_universe, sizeof(double) * UWIDTH * UHEIGHT);
	memcpy(old_hx_universe, current_hx_universe, sizeof(double) * UWIDTH * UHEIGHT);
	memcpy(old_hy_universe, current_hy_universe, sizeof(double) * UWIDTH * UHEIGHT);
	current_time += time_step;
	clock_t end = clock();

	printf("%.16f\n", (double)(end - start) / CLOCKS_PER_SEC);
}

int main()
{
	init_params();
	set_pml_geometry();
	calculate_factors();
	while (1)
	{
		update_fields();
	}
	return 0;
}
