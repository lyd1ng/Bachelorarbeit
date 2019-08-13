for (float current_time = 0; current_time < max_time; current_time += time_step)
{
	for (int z=1; z < UDEPTH - 1; z++)
	{
		for (int y=1; y < UHEIGHT - 1; y++)
		{
			for (int x=1; x < UWIDTH - 1; x++)
			{
				fp hx_curl_term = ((old_ez_universe[x][y + 1][z] - old_ez_universe[x][y][z]) / grid_width_y) - ((old_ey_universe[x][y][z + 1] - old_ey_universe[x][y][z]) / grid_width_z);
				fp hy_curl_term = ((old_ex_universe[x][y][z + 1] - old_ex_universe[x][y][z]) / grid_width_z) - ((old_ez_universe[x + 1][y][z] - old_ez_universe[x][y][z]) / grid_width_x);
				fp hz_curl_term = ((old_ey_universe[x + 1][y][z] - old_ey_universe[x][y][z]) / grid_width_x) - ((old_ex_universe[x][y + 1][z] - old_ex_universe[x][y][z]) / grid_width_y);
				hx_curl_integral_universe[x][y][z] += hx_curl_term;
				hy_curl_integral_universe[x][y][z] += hy_curl_term;
				hz_curl_integral_universe[x][y][z] += hz_curl_term;
				hx_field_integral_universe[x][y][z] += old_hx_universe[x][y][z];
				hy_field_integral_universe[x][y][z] += old_hy_universe[x][y][z];
				hz_field_integral_universe[x][y][z] += old_hz_universe[x][y][z];
				current_hx_universe[x][y][z] = + hx_factor1_universe[x][y][z] * old_hx_universe[x][y][z] - hx_factor2_universe[x][y][z] * hx_curl_term
											   - hx_factor3_universe[x][y][z] * hx_curl_integral_universe[x][y][z] - hx_factor4_universe[x][y][z] * hx_field_integral_universe[x][y][z];
				current_hy_universe[x][y][z] = + hy_factor1_universe[x][y][z] * old_hy_universe[x][y][z] - hy_factor2_universe[x][y][z] * hy_curl_term
											   - hy_factor3_universe[x][y][z] * hy_curl_integral_universe[x][y][z] - hy_factor4_universe[x][y][z] * hy_field_integral_universe[x][y][z];
				current_hz_universe[x][y][z] = + hz_factor1_universe[x][y][z] * old_hz_universe[x][y][z] - hz_factor2_universe[x][y][z] * hz_curl_term
											   - hz_factor3_universe[x][y][z] * hz_curl_integral_universe[x][y][z] - hz_factor4_universe[x][y][z] * hz_field_integral_universe[x][y][z];
			}
		}
	}
	handle_perfect_magnetical_conductors();


	// Update the electrical field
	for (int z=1; z < UDEPTH - 1; z++)
	{
		for (int y=1; y < UHEIGHT - 1; y++)
		{
			for (int x=1; x < UWIDTH - 1; x++)
			{
				fp ex_curl_term = ((current_hz_universe[x][y][z] - current_hz_universe[x][y - 1][z]) / grid_width_y) - ((current_hy_universe[x][y][z] - current_hy_universe[x][y][z - 1]) / grid_width_z);
				fp ey_curl_term = ((current_hx_universe[x][y][z] - current_hx_universe[x][y][z - 1]) / grid_width_z) - ((current_hz_universe[x][y][z] - current_hz_universe[x - 1][y][z]) / grid_width_x);
				fp ez_curl_term = ((current_hy_universe[x][y][z] - current_hy_universe[x - 1][y][z]) / grid_width_x) - ((current_hx_universe[x][y][z] - current_hx_universe[x][y - 1][z]) / grid_width_y);
				ex_curl_integral_universe[x][y][z] += ex_curl_term;
				ey_curl_integral_universe[x][y][z] += ey_curl_term;
				ez_curl_integral_universe[x][y][z] += ez_curl_term;
				ex_field_integral_universe[x][y][z] += old_ex_universe[x][y][z];
				ey_field_integral_universe[x][y][z] += old_ey_universe[x][y][z];
				ez_field_integral_universe[x][y][z] += old_ez_universe[x][y][z];
				current_ex_universe[x][y][z] = ex_factor1_universe[x][y][z] * old_ex_universe[x][y][z] - ex_factor2_universe[x][y][z] * ex_field_integral_universe[x][y][z]
											   + ex_factor3_universe[x][y][z] * ex_curl_term + ex_factor4_universe[x][y][z] * ex_curl_integral_universe[x][y][z];
				current_ey_universe[x][y][z] = ey_factor1_universe[x][y][z] * old_ey_universe[x][y][z] - ey_factor2_universe[x][y][z] * ey_field_integral_universe[x][y][z]
											   + ey_factor3_universe[x][y][z] * ey_curl_term + ey_factor4_universe[x][y][z] * ey_curl_integral_universe[x][y][z];
				current_ez_universe[x][y][z] = ez_factor1_universe[x][y][z] * old_ez_universe[x][y][z] - ez_factor2_universe[x][y][z] * ez_field_integral_universe[x][y][z]
											   + ez_factor3_universe[x][y][z] * ez_curl_term + ez_factor4_universe[x][y][z] * ez_curl_integral_universe[x][y][z];
			}
		}
	}

	handle_perfect_electrical_conductors();

	int middle_x = (int)UWIDTH / 2;
	int middle_y = (int)UHEIGHT / 2;
	int middle_z = (int)UDEPTH / 2;
	current_ez_universe[middle_x][middle_y][middle_z] = ramped_sinus(current_time, ramped_sin_length * time_step, ramped_sin_amplitude, ramped_sin_frequency);

	memcpy(old_ex_universe, current_ex_universe, sizeof(fp) * UWIDTH * UHEIGHT * UDEPTH);
	memcpy(old_ey_universe, current_ey_universe, sizeof(fp) * UWIDTH * UHEIGHT * UDEPTH);
	memcpy(old_ez_universe, current_ez_universe, sizeof(fp) * UWIDTH * UHEIGHT * UDEPTH);

	memcpy(old_hx_universe, current_hx_universe, sizeof(fp) * UWIDTH * UHEIGHT * UDEPTH);
	memcpy(old_hy_universe, current_hy_universe, sizeof(fp) * UWIDTH * UHEIGHT * UDEPTH);
	memcpy(old_hz_universe, current_hz_universe, sizeof(fp) * UWIDTH * UHEIGHT * UDEPTH);
}
