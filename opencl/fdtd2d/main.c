#include <math.h>
#include <stdio.h>
#include <alloca.h>
#include <CL/opencl.h>
#include "c_helper.h"
#include "fdtd_helper.h"

#define NX 4
#define NY 4
#define N NX * NY
#define BPP 4
#define IMG_SIZE N * BPP * sizeof(float)


void verbose(int error, char* msg)
{
	if (error) { printf("%s\n", msg); }
}



void set_pml_geometry(const int uwidth, const int uheight, const float time_step, float sigmax[NY][NX], float sigmay[NY][NX])
{
    int pml_width = 20;
    int pml_height = 20;
    double f = vacuum_permittivity / (2.0 * time_step);
    for (int y=1; y < uheight - 1; y++)
    {
        for (int x=1; x < pml_width; x++)
	{
            sigmax[y][pml_width - x - 1] = f * pow(x / pml_width, 3);
            sigmax[y][uwidth - pml_width + x] = f * pow(x / pml_width, 3);
	}
    }
    for (int x=1; x < uwidth - 1; x++)
    {
        for (int y=1; y < pml_height; y++)
	{
            sigmay[pml_height - y - 1][x] = f * pow(y / pml_height, 3);
            sigmay[uheight - pml_height + y][x] = f * pow(y / pml_height, 3);
	}
    }
}

void calculate_factors(const int uwidth,
		       const int uheight,
		       float time_step,
		       float permittivity[NY][NX],
		       float permeability[NY][NX],
		       float sigmax[NY][NX],
		       float sigmay[NY][NX],
		       float ez_factor1[NY][NX],
		       float ez_factor2[NY][NX],
		       float ez_factor3[NY][NX],
		       float hx_factor1[NY][NX],
		       float hx_factor2[NY][NX],
		       float hx_factor3[NY][NX],
		       float hy_factor1[NY][NX],
		       float hy_factor2[NY][NX],
		       float hy_factor3[NY][NX])
{
    for (int y=0; y < uheight; y++)
    {
        for (int x=0; x < uwidth; x++)
	{
		double nominator = 0;
		double denominator = 0;
            	ez_factor1[y][x] = 1 - ((time_step * (sigmax[y][x] + sigmay[y][x]))
				   / vacuum_permittivity);
            	ez_factor2[y][x] = (sigmax[y][x] * sigmay[y][x] * pow(time_step, 2))
				   / pow(vacuum_permittivity, 2);
            	ez_factor3[y][x] = (light_speed * time_step) / permittivity[y][x];

            	hx_factor1[y][x] = (1 / time_step - (sigmay[y][x]
				   / (2 * vacuum_permittivity))) / (1 / time_step
				   + (sigmay[y][x] / (2 * vacuum_permittivity)));
            	hx_factor2[y][x] = (light_speed / permittivity[y][x])
				   / ((1 / time_step) + (sigmay[y][x]
				   / (2 * vacuum_permittivity)));
            	nominator = (light_speed * sigmax[y][x] * time_step)
			    / (permeability[y][x] * vacuum_permittivity);
            	denominator = 1 / time_step + (sigmay[y][x]
			      / (2 * vacuum_permittivity));
            	hx_factor3[y][x] = nominator / denominator;

            	hy_factor1[y][x] = (1 / time_step - (sigmax[y][x]
				   / (2 * vacuum_permittivity)))
				   / (1 / time_step + (sigmax[y][x]
				   / (2 * vacuum_permittivity)));
            	hy_factor2[y][x] = (light_speed / permittivity[y][x])
				   / ((1 / time_step) + (sigmax[y][x]
				   / (2 * vacuum_permittivity)));
            	nominator = (light_speed * sigmay[y][x] * time_step)
			    / (permeability[y][x] * vacuum_permittivity);
            	denominator = 1 / time_step + (sigmax[y][x]
			      / (2 * vacuum_permittivity));
            	hy_factor3[y][x] = nominator / denominator;
	}
    }
}



int main()
{
	cl_int errNum;
	cl_kernel kernel;
	cl_program program;
	cl_context context;
	cl_device_id device;
	cl_command_queue queue;
	cl_platform_id platform;

	// The necessary arrays
	float ez_old[NY][NX] = {0};
	float hx_old[NY][NX] = {0};
	float hy_old[NY][NX] = {0};
	float ez_current[NY][NX] = {0};
	float hx_current[NY][NX] = {0};
	float hy_current[NY][NX] = {0};
	float ez_integral[NY][NX] = {0};
	float hx_integral[NY][NX] = {0};
	float hy_integral[NY][NX] = {0};
	float ez_factor1[NY][NX] = {0};
	float ez_factor2[NY][NX] = {0};
	float ez_factor3[NY][NX] = {0};
	float hx_factor1[NY][NX] = {0};
	float hx_factor2[NY][NX] = {0};
	float hx_factor3[NY][NX] = {0};
	float hy_factor1[NY][NX] = {0};
	float hy_factor2[NY][NX] = {0};
	float hy_factor3[NY][NX] = {0};
	float sigmax[NY][NX] = {0};
	float sigmay[NY][NX] = {0};
	float permittivity[NY][NX] = {0};
	float permeability[NY][NX] = {0};

	// The permeability and permittivity are relativ
	// to the free space values. So 1 is the correct value
	// to simulate the free space.
	for (int y = 0; y < NY; y++)
	{
		for (int x = 0; x < NX; x++)
		{
			permittivity[y][x] = 1;
			permeability[y][x] = 1;
		}
	}
	printf("Array Initialisation passed\n");

	// The simulation parameters
	double max_frequency = 1.5 * pow(10, 6);
	double fraction = 1;
	int samplerate_space = 40;
	int samplerate_time = 10;
	double grid_width_x = get_grid_width(max_frequency, fraction, samplerate_space);
	double grid_width_y = get_grid_width(max_frequency, fraction, samplerate_space);
	double time_step = get_time_resolution_2d(fraction, grid_width_x, grid_width_y)
			  / (double)samplerate_time;
	double current_time = 0;
	double max_time = 1 * time_step;
	printf("Parameter Initialisation passed\n");


	// Set the pml geometry
	// set_pml_geometry(NX, NY, time_step, sigmax, sigmay);
	printf("Set PML passed\n");

	// Calculate factors
	calculate_factors(NX, NY, time_step, permittivity, permeability, sigmax, sigmay,
			ez_factor1, ez_factor2, ez_factor3,
			hx_factor1, hx_factor2, hx_factor3,
			hy_factor1, hy_factor2, hy_factor3);
	
	for (int y=0; y < NY; y++)
	{
		for(int x = 0; x < NX; x++)
		{
			printf("%f ", ez_factor1[y][x]);
			/*
			printf("%f ", ez_factor2[y][x]);
			printf("%f ", ez_factor3[y][x]);
			printf("%f ", hx_factor1[y][x]);
			printf("%f ", hx_factor2[y][x]);
			printf("%f ", hx_factor3[y][x]);
			printf("%f ", hy_factor1[y][x]);
			printf("%f ", hy_factor2[y][x]);
			printf("%f ", hy_factor3[y][x]);
			*/
		}
		printf("\n");
	}


	// Get the first platfrom available
	errNum = clGetPlatformIDs(1, &platform, NULL);
	verbose(errNum, "Error getting the platform");

	// Get a GPU device
	errNum = clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 1, &device, NULL);
	verbose(errNum, "Error getting the device");

	// Create a Context
	cl_context_properties context_properties[] =
	{
			CL_CONTEXT_PLATFORM, (cl_context_properties)platform, 0
	};
	context = clCreateContext(context_properties, 1, &device, NULL, NULL, &errNum);
	verbose(errNum, "Error creating the context");

	// Create the command queue
	queue = clCreateCommandQueue(context, device, 0, &errNum);
	verbose(errNum, "Error creating the command queue");

	// Create and build the program
	FILE* code_fd = fopen("fdtd2d.c", "r");
	char *code;
	size_t code_length;
	read_fileh(code_fd, &code, &code_length);
	program = clCreateProgramWithSource(context, 1,
			(const char**)&code, &code_length, &errNum);
	verbose(errNum, "Error creating the program");
	errNum = clBuildProgram(program, 1, &device, NULL, NULL, NULL);
	free(code);

	if (errNum)
	{
		size_t log_size;
		char* log;
		clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, 0,
				NULL, &log_size);
		log = (char*)alloca(log_size * sizeof(char));
		clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, log_size,
				log, NULL);
		printf("%s\n", log);
		return 0;
	}

	// Oneshot source
	ez_old[NY / 2][NX / 2] = 0.00001;

	// Create the buffers
	// Old Buffers (read only)
	int n = (NX - 1) * (NY - 1);
	cl_mem ez_oldb = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,
				 n * sizeof(float), &ez_old[1][1], &errNum);
	verbose(errNum, "Error creating ez_oldb");
	cl_mem hx_oldb = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,
				 n * sizeof(float), &hx_old[1][1], &errNum);
	verbose(errNum, "Error creating hx_oldb");
	cl_mem hy_oldb = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,
				 n * sizeof(float), &hy_old[1][1], &errNum);
	verbose(errNum, "Error creating hy_oldb");

	// Current Buffers (write only)
	cl_mem ez_currentb = clCreateBuffer(context, CL_MEM_WRITE_ONLY
			| CL_MEM_USE_HOST_PTR, n * sizeof(float), &ez_current[1][1],
			&errNum);
	verbose(errNum, "Error creating ez_currentb");
	cl_mem hx_currentb = clCreateBuffer(context, CL_MEM_WRITE_ONLY
			| CL_MEM_USE_HOST_PTR, n * sizeof(float), &hx_current[1][1],
			&errNum);
	verbose(errNum, "Error creating hx_currentb");
	cl_mem hy_currentb = clCreateBuffer(context, CL_MEM_WRITE_ONLY
			| CL_MEM_USE_HOST_PTR, n * sizeof(float), &hy_current[1][1],
			&errNum);
	verbose(errNum, "Error creating hy_currentb");

	// Integral (read only)
	cl_mem ez_integralb = clCreateBuffer(context, CL_MEM_READ_ONLY
			| CL_MEM_USE_HOST_PTR, n * sizeof(float), &ez_integral[1][1],
			&errNum);
	verbose(errNum, "Error creating ez_integralb");
	cl_mem hx_integralb = clCreateBuffer(context, CL_MEM_READ_ONLY
			| CL_MEM_USE_HOST_PTR, n * sizeof(float), &hx_integral[1][1],
			&errNum);
	verbose(errNum, "Error creating hx_integralb");
	cl_mem hy_integralb = clCreateBuffer(context, CL_MEM_READ_ONLY
			| CL_MEM_USE_HOST_PTR, n * sizeof(float), &hy_integral[1][1],
			&errNum);
	verbose(errNum, "Error creating hy_integralb");

	// Factors (read only)
	cl_mem ez_factor1b = clCreateBuffer(context, CL_MEM_READ_ONLY
			| CL_MEM_USE_HOST_PTR, n * sizeof(float), &ez_factor1[1][1],
			&errNum);
	verbose(errNum, "Error creating ez_factor1b");
	cl_mem ez_factor2b = clCreateBuffer(context, CL_MEM_READ_ONLY
			| CL_MEM_USE_HOST_PTR, n * sizeof(float), &ez_factor2[1][1],
			&errNum);
	verbose(errNum, "Error creating ez_factor2b");
	cl_mem ez_factor3b = clCreateBuffer(context, CL_MEM_READ_ONLY
			| CL_MEM_USE_HOST_PTR, n * sizeof(float), &ez_factor3[1][1],
			&errNum);
	verbose(errNum, "Error creating ez_factor3b");
	cl_mem hx_factor1b = clCreateBuffer(context, CL_MEM_READ_ONLY
			| CL_MEM_USE_HOST_PTR, n * sizeof(float), &hx_factor1[1][1],
			&errNum);
	verbose(errNum, "Error creating hx_factor1b");
	cl_mem hx_factor2b = clCreateBuffer(context, CL_MEM_READ_ONLY
			| CL_MEM_USE_HOST_PTR, n * sizeof(float), &hx_factor2[1][1],
			&errNum);
	verbose(errNum, "Error creating hx_factor2b");
	cl_mem hx_factor3b = clCreateBuffer(context, CL_MEM_READ_ONLY
			| CL_MEM_USE_HOST_PTR, n * sizeof(float), &hx_factor3[1][1],
			&errNum);
	verbose(errNum, "Error creating hx_factor3b");
	cl_mem hy_factor1b = clCreateBuffer(context, CL_MEM_READ_ONLY
			| CL_MEM_USE_HOST_PTR, n * sizeof(float), &hy_factor1[1][1],
			&errNum);
	verbose(errNum, "Error creating hy_factor1b");
	cl_mem hy_factor2b = clCreateBuffer(context, CL_MEM_READ_ONLY
			| CL_MEM_USE_HOST_PTR, n * sizeof(float), &hy_factor2[1][1],
			&errNum);
	verbose(errNum, "Error creating hy_factor2b");
	cl_mem hy_factor3b = clCreateBuffer(context, CL_MEM_READ_ONLY
			| CL_MEM_USE_HOST_PTR, n * sizeof(float), &hy_factor3[1][1],
			&errNum);
	verbose(errNum, "Error creating hy_factor3b");
	cl_mem sigmaxb = clCreateBuffer(context, CL_MEM_READ_ONLY
			| CL_MEM_USE_HOST_PTR, n * sizeof(float), &sigmax[1][1],
			&errNum);
	verbose(errNum, "Error creating sigmaxb");
	cl_mem sigmayb = clCreateBuffer(context, CL_MEM_READ_ONLY
			| CL_MEM_USE_HOST_PTR, n * sizeof(float), &sigmay[1][1],
			&errNum);
	verbose(errNum, "Error creating sigmayb");

	// Pass the arrays to the device
	errNum = clEnqueueReadBuffer(queue, ez_oldb, CL_TRUE, 0, 

	// Create the kernel object
	kernel = clCreateKernel(program, "fdtd2d_noiter", &errNum);
	verbose(errNum, "Error creating kernel");
	
	// Set the kernel args
	int uwidth = NX;
	errNum = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void*)&ez_oldb);
	verbose(errNum, "Error setting ez_oldb arg");
	errNum = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void*)&hx_oldb);
	verbose(errNum, "Error setting hx_oldb arg");
	errNum = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void*)&hy_oldb);
	verbose(errNum, "Error setting hy_oldb arg");
	errNum = clSetKernelArg(kernel, 3, sizeof(cl_mem), (void*)&ez_currentb);
	verbose(errNum, "Error setting ez_currentb arg");
	errNum = clSetKernelArg(kernel, 4, sizeof(cl_mem), (void*)&hx_currentb);
	verbose(errNum, "Error setting hx_currentb arg");
	errNum = clSetKernelArg(kernel, 5, sizeof(cl_mem), (void*)&hy_currentb);
	verbose(errNum, "Error setting hy_currentb arg");
	errNum = clSetKernelArg(kernel, 6, sizeof(cl_mem), (void*)&ez_integralb);
	verbose(errNum, "Error setting ez_integralb arg");
	errNum = clSetKernelArg(kernel, 7, sizeof(cl_mem), (void*)&hx_integralb);
	verbose(errNum, "Error setting hx_integralb arg");
	errNum = clSetKernelArg(kernel, 8, sizeof(cl_mem), (void*)&hy_integralb);
	verbose(errNum, "Error setting hy_integralb arg");
	errNum = clSetKernelArg(kernel, 9, sizeof(cl_mem), (void*)&ez_factor1b);
	verbose(errNum, "Error setting ez_factor1b arg");
	errNum = clSetKernelArg(kernel, 10, sizeof(cl_mem), (void*)&ez_factor2b);
	verbose(errNum, "Error setting ez_factor2b arg");
	errNum = clSetKernelArg(kernel, 11, sizeof(cl_mem), (void*)&ez_factor3b);
	verbose(errNum, "Error setting ez_factor3b arg");
	errNum = clSetKernelArg(kernel, 12, sizeof(cl_mem), (void*)&hx_factor1b);
	verbose(errNum, "Error setting hx_factor1b arg");
	errNum = clSetKernelArg(kernel, 13, sizeof(cl_mem), (void*)&hx_factor2b);
	verbose(errNum, "Error setting hx_factor2b arg");
	errNum = clSetKernelArg(kernel, 14, sizeof(cl_mem), (void*)&hx_factor3b);
	verbose(errNum, "Error setting hx_factor3b arg");
	errNum = clSetKernelArg(kernel, 15, sizeof(cl_mem), (void*)&hy_factor1b);
	verbose(errNum, "Error setting hy_factor1b arg");
	errNum = clSetKernelArg(kernel, 16, sizeof(cl_mem), (void*)&hy_factor2b);
	verbose(errNum, "Error setting hy_factor2b arg");
	errNum = clSetKernelArg(kernel, 17, sizeof(cl_mem), (void*)&hy_factor3b);
	verbose(errNum, "Error setting hy_factor3b arg");
	errNum = clSetKernelArg(kernel, 18, sizeof(float), (void*)&grid_width_x);
	verbose(errNum, "Error setting grid_width_x");
	errNum = clSetKernelArg(kernel, 19, sizeof(float), (void*)&grid_width_y);
	verbose(errNum, "Error setting grid_width_y");
	errNum = clSetKernelArg(kernel, 20, sizeof(int), (void*)&uwidth);
	verbose(errNum, "Error setting uwidth");

	// Invoke the kernel
	size_t gws[2] = {NX - 2, NY - 2};
	while (current_time < max_time)
	{
		errNum = clEnqueueNDRangeKernel(queue, kernel, 2, NULL,
				(size_t*)&gws, NULL, 0, 0, 0);
		errNum = clEnqueueCopyBuffer(queue, ez_currentb, ez_oldb, 0, 0,
				N * sizeof(float), 0, NULL, NULL);
		errNum = clEnqueueCopyBuffer(queue, hx_currentb, hx_oldb, 0, 0,
				N * sizeof(float), 0, NULL, NULL);
		errNum = clEnqueueCopyBuffer(queue, hy_currentb, hy_oldb, 0, 0,
				N * sizeof(float), 0, NULL, NULL);
		clFinish(queue);
		current_time += time_step;
	}
	
	errNum = clEnqueueReadBuffer(queue, ez_currentb, CL_TRUE, 0, N * sizeof(float),
		(void*)(&ez_current[1][1]), 0, NULL, NULL);
	
	// Print the results in an ugly list
	printf("\n\nThe results are:\n");
	for (int y = 0; y < NY; y++)
	{
		for (int x = 0; x < NX; x++)
		{
			printf("%f ", ez_current[y][x]);
		}
		printf("\n");
	}
	return 0;
}
