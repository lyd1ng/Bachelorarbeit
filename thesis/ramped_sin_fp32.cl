inline float ramped_sin_fp32(float t, float ramp_len, float frequency)
{
	return smoothstep(0, ramp_len, t) * sin(2.0 * M_PI * frequency * t);
}
