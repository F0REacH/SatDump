#include "opencl.h"
#include "logger.h"
#include "core/config.h"
#include <fstream>

#ifdef USE_OPENCL
namespace satdump
{
    namespace opencl
    {
        std::vector<OCLDevice> getAllDevices()
        {
            std::vector<OCLDevice> devs;
            cl_platform_id platforms_ids[100];
            cl_uint platforms_cnt = 0;
            cl_device_id devices_ids[100];
            cl_uint devices_cnt = 0;
            char device_name[200];
            size_t device_name_len = 0;

            if (clGetPlatformIDs(100, platforms_ids, &platforms_cnt) != CL_SUCCESS)
                return devs;

            for (int p = 0; p < (int)platforms_cnt; p++)
            {
                if (clGetDeviceIDs(platforms_ids[p], CL_DEVICE_TYPE_ALL, 100, devices_ids, &devices_cnt) != CL_SUCCESS)
                    continue;

                for (int d = 0; d < (int)devices_cnt; d++)
                    if (clGetDeviceInfo(devices_ids[d], CL_DEVICE_NAME, 200, device_name, &device_name_len) == CL_SUCCESS)
                        devs.push_back({p, d, std::string(&device_name[0], &device_name[device_name_len])});
            }

            return devs;
        }

        cl_context_properties opencl_ctx_props[7] = {0, 0, 0, 0, 0, 0, 0};

        void initOpenCL()
        {
            std::vector<OCLDevice> devices = getAllDevices();
            logger->info("Found OpenCL Devices ({:d}) :", devices.size());
            for (OCLDevice &d : devices)
                logger->debug(" - " + d.name);
        }

        bool context_is_init = false;
        cl_context ocl_context;
        cl_device_id ocl_device;

        void setupOCLContext()
        {
            int platform_id = satdump::config::main_cfg["satdump_general"]["opencl_device"]["platform"].get<int>();
            int device_id = satdump::config::main_cfg["satdump_general"]["opencl_device"]["device"].get<int>();

            cl_platform_id platforms_ids[100];
            cl_uint platforms_cnt = 0;
            cl_device_id devices_ids[100];
            cl_uint devices_cnt = 0;
            char device_platform_name[200];
            size_t device_platform_name_len = 0;
            cl_int err = 0;

            if (!context_is_init)
            {
                logger->trace("First OpenCL context request. Initializing...");

                if (clGetPlatformIDs(100, platforms_ids, &platforms_cnt) != CL_SUCCESS)
                {
                    logger->error("Could not get OpenCL platform IDs!");
                    return;
                }

                if (platforms_cnt == 0)
                    std::runtime_error("No platforms found. Check OpenCL installation!");

                cl_platform_id platform = platforms_ids[platform_id];
                if (clGetPlatformInfo(platform, CL_PLATFORM_NAME, 200, device_platform_name, &device_platform_name_len) == CL_SUCCESS)
                    logger->info("Using platform: {:s}", std::string(&device_platform_name[0], &device_platform_name[device_platform_name_len]));
                else
                    logger->error("Could not get platform name!");

                if (clGetDeviceIDs(platform, CL_DEVICE_TYPE_ALL, 100, devices_ids, &devices_cnt) != CL_SUCCESS)
                {
                    logger->error("Could not get OpenCL devices IDs!");
                    return;
                }

                if (devices_cnt == 0)
                    std::runtime_error("No devices found. Check OpenCL installation!");

                ocl_device = devices_ids[device_id];
                if (clGetDeviceInfo(ocl_device, CL_DEVICE_NAME, 200, device_platform_name, &device_platform_name_len) == CL_SUCCESS)
                    logger->info("Using device: {:s}", std::string(&device_platform_name[0], &device_platform_name[device_platform_name_len]));

                ocl_context = clCreateContext(opencl_ctx_props, 1, &ocl_device, NULL, NULL, &err);
                if (err != CL_SUCCESS)
                {
                    logger->error("Could not init OpenCL context!");
                    return;
                }

                context_is_init = true;
            }
            else
            {
                logger->trace("OpenCL context already initilized.");
            }

            if (!context_is_init)
                throw std::runtime_error("OpenCL context not initialized!");
        }

        std::map<std::string, cl_program> cached_kernels;

        cl_program buildCLKernel(std::string path, bool use_cache)
        {
            if (use_cache)                          // If cache enabled...
                if (cached_kernels.count(path) > 0) // ...check if we already have this kernel
                    return cached_kernels[path];

            std::ifstream isf(path);
            std::string kernel_src(std::istreambuf_iterator<char>{isf}, {});

            const char *srcs[1] = {kernel_src.c_str()};
            const size_t lens[1] = {kernel_src.length()};
            cl_int err = 0;
            char error_msg[10000];
            size_t error_len = 0;

            cl_program prg = clCreateProgramWithSource(ocl_context, 1, srcs, lens, &err);
            err = clBuildProgram(prg, 1, &ocl_device, NULL, NULL, NULL);

            if (err != CL_SUCCESS)
            {
                if (clGetProgramBuildInfo(prg, ocl_device, CL_PROGRAM_BUILD_LOG, 10000, error_msg, &error_len) == CL_SUCCESS)
                    throw std::runtime_error("Error building: " + std::string(&error_msg[0], &error_msg[error_len]));
                else
                    throw std::runtime_error("Error building, and could not read error log!");
            }

            if (use_cache)                              // If cache enabled...
                if (cached_kernels.count(path) == 0)    // ...and we don't already have the kernel...
                    cached_kernels.insert({path, prg}); // ...return it

            return prg;
        }
    }
}
#endif