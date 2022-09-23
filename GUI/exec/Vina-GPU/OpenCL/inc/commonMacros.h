#ifndef COMMON_MACROS_H
#define COMMON_MACROS_H
#define _CRT_SECURE_NO_WARNINGS
//#define CL_USE_DEPRECATED_OPENCL_3_0_APIS

//#define OPENCL_PART_1
#define OPENCL_PART_2

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <fstream>
#ifdef __APPLE__
	#include "OpenCL/opencl.h"
#else
	#include "CL/opencl.h"
#endif

#ifdef WIN32
#include <Windows.h>
#else
#include <unistd.h>
#endif

#include "kernel2.h"
#define KERNEL_1_NUM_OF_KERNELS 1 //user define (FIX IT!)
#define KERNEL_2_NUM_OF_KERNELS 2
#define NUM_OF_FILES 5
#endif

