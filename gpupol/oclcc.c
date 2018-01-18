#include <stdio.h>
#if __APPLE__
	#include <OpenCL/cl.h>
#else
	#include <CL/cl.h>
#endif
#include <string.h>
#include "util.h"

int main(int argc, char** argv){
	if((argc < 4 && argc>5) || (strcmp(argv[1],"cpu") &&  strcmp(argv[1],"gpu"))){
		printf("Error: please specify device, source file and output file: [cpu|gpu] srcFile outFile [options]\n");
		return 192;
	}
	
	char* targetDeviceType=argv[1];
	char* sourceFile=argv[2];
	char* outputFile=argv[3];
	char compOptions[1000];
// 	if(argc==5)
// 		compOptions=argv[4];
// 	else
// 		compOptions="";
	
	const char* src= CStringFromFile(sourceFile);
	cl_platform_id platformId;
	cl_uint nPlatforms;
	char tString[1000];
	cl_device_id deviceIds[3];
	cl_device_id device;
	cl_device_type deviceType;
	cl_uint nDevices;
	cl_int ret;
	cl_context context;
	cl_program program;
	int gpu=-1;
	int cpu=-1;
	int i;
	char deviceVendor[1000];
	
	clGetPlatformIDs(1, &platformId, &nPlatforms);
	printf("Number of platforms: %i\n", nPlatforms);
	clGetPlatformInfo(platformId, CL_PLATFORM_VENDOR, 1000, (void*) tString, NULL);
	printf("Platform vendor name: %s\n", tString);
	clGetDeviceIDs(platformId, CL_DEVICE_TYPE_ALL, 3, deviceIds, &nDevices);
	printf("Number of devices: %i\n", nDevices);
	for(i=0; i<nDevices; i++){
		ret=clGetDeviceInfo(deviceIds[i], CL_DEVICE_TYPE, sizeof(deviceType), (void*) &deviceType, NULL);
		if(deviceType==CL_DEVICE_TYPE_CPU) {
			cpu=i;
		}
		else if(deviceType==CL_DEVICE_TYPE_GPU) {
			gpu=i;
		}
		else sprintf(tString, "other");
		ret=clGetDeviceInfo(deviceIds[i], CL_DEVICE_NAME, 1000, (void*) tString, NULL);
		ret=clGetDeviceInfo(deviceIds[i], CL_DEVICE_VENDOR, 1000, (void*) tString, NULL);
		printf("Vendor: %s\n", tString);
	}
	if(!strcmp(targetDeviceType,"gpu")){
		if(gpu==-1){
			printf("No GPU device found! Aborting...\n");
			return 0;
		}
		ret=clGetDeviceInfo(deviceIds[gpu], CL_DEVICE_VENDOR, 1000, (void*) deviceVendor, NULL);
		device = deviceIds[gpu];
		if(!strcmp(deviceVendor, "Advanced Micro Devices, Inc."))
			sprintf(compOptions,"-D ARCH=__ATI_GPU__");
		else if(!strcmp(deviceVendor, "NVIDIA Corporation"))
			sprintf(compOptions,"-D ARCH=__NVIDIA_GPU__");
		else
			compOptions[0]='\0';
	}
	else{
		if(cpu==-1){
			printf("No CPU device found! Aborting...\n");
			return 0;
		}
		sprintf(compOptions,"-D ARCH=__CPU__");
		device = deviceIds[cpu];
	}
	context = clCreateContext(NULL, 1, &device, NULL, NULL, &ret);
	if(ret==CL_SUCCESS)
		printf("Succesfully created context\n");
	else
		return 193;
	
	program=clCreateProgramWithSource(context, 1, &src, NULL, &ret);
	if(ret==CL_SUCCESS)
		printf("Succesfully created program object\n");
	else{
		char buf[0x10000];
		clGetProgramBuildInfo( program, device, CL_PROGRAM_BUILD_LOG, 0x10000, buf, NULL);
		printf("\n%s\n", buf);
		return 194;
	}
	
	sprintf(compOptions, "%s -DIS_MAIN", compOptions);
	ret = clBuildProgram(program, 1, &device, compOptions, NULL, NULL);
	if(ret==CL_SUCCESS)
		printf("Succesfully built program\n");
	else{
		char buf[0x100000];
		clGetProgramBuildInfo( program, device, CL_PROGRAM_BUILD_LOG, 0x100000, buf, NULL);
		printf("\n%s\n", buf);
// 		printf("\n%s\n", src);
		return 195;	
	}
	size_t progSize;
	cl_uint nDev;
	char* progBin;
	clGetProgramInfo(program, CL_PROGRAM_NUM_DEVICES, sizeof(cl_uint), &nDev, NULL);
	clGetProgramInfo(program, CL_PROGRAM_BINARY_SIZES, sizeof(size_t), (void*)&progSize, NULL);
	progBin = (char*) malloc(2*progSize);
	clGetProgramInfo(program, CL_PROGRAM_BINARIES, sizeof(char*), &progBin, NULL);
	FILE* pFile;
	
	pFile = fopen(outputFile, "wb");
	fwrite(progBin, sizeof(char), progSize, pFile);
	fclose(pFile);
	clReleaseContext(context);
	return 0;
}
	
