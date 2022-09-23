# Vina-GPU+
Vina-GPU+ further accelerates Vina-GPU and facilitates single receptor-multi-ligand docking.
## Compiling and Running
**Note**: at least one GPU card is required and make sure the version of GPU driver is up to date
### Windows
#### Run on the executable file
1. For the first time to use Vina-GPU+, please run `Vina-GPU+_K.exe` with command `./Vina-GPU+_K.exe --config=./input_file_example/2bm2_config.txt`
You are supposed to have the docking results `test_out` file ,`Kernel1_Opt.bin` file and `Kernel2_Opt.bin` file
2. Once you have the `Kernel1_Opt.bin` file and `Kernel2_Opt.bin` file, you can run `Vina-GPU+.exe` without compiling the kernel files (thus to save more runtime)
>When you run `Vina-GPU+.exe`, please make sure `Kernel1_Opt.bin` file and `Kernel2_Opt.bin` file are in the same directory

For the usage and limitaiton of Vina-GPU+, please check [Usage](#Usage) and [Limitation](#Limitation).
A graphic user interface (GUI) is also provided for Windows users, please check [GUI](#GUI)
#### Build from source file
>Visual Studio 2019 is recommended for build Vina-GPU from source
1. install [boost library](https://www.boost.org/) (current version is 1.77.0)
2. install [CUDA Toolkit](https://developer.nvidia.com/zh-cn/cuda-toolkit) (current version is v11.5) if you are using NVIDIA GPU cards

    Note: the OpenCL library can be found in CUDA installation path for NVIDIA or in the driver installation path for AMD

3. add `./lib` `./OpenCL/inc` `$(YOUR_BOOST_LIBRARY_PATH)` `$(YOUR_BOOST_LIBRARY_PATH)/boost` `$(YOUR_CUDA_TOOLKIT_LIBRARY_PATH)/CUDA/v11.5/include` in the include directories
4. add `$(YOUR_BOOST_LIBRARY_PATH)/stage/lib` `$(YOUR_CUDA_TOOLKIT_PATH)/CUDA/lib/x64`in the addtional library 
5. add `OpenCL.lib` in the additional dependencies 
6. add `--config=./input_file_example/2bm2_config.txt` in the command arguments
7.  add `NVIDIA_PLATFORM` `OPENCL_3_0` `WINDOWS` in the preprocessor definitions if necessary
8. if you want to compile the binary kernel file on the fly, add `BUILD_KERNEL_FROM_SOURCE` in the preprocessor definitions
9. build & run

Note: ensure the line ending are CLRF
## Usage
|Arguments| Description|Default value
|--|--|--|
|--config | the config file (in .txt format) that contains all the following arguments for the convenience of use| no default
| --receptor | the recrptor file (in .pdbqt format)| no default
|--ligand_directory| this path contains all the ligand files，the ligand file (in .pdbqt fotmat)| no default
|--thread| the scale of parallelism (docking lanes)|1000
|--search_depth| the number of searching iterations in each docking lane| heuristically determined
|--center_x/y/z|the center of searching box in the receptor|no default
|--size_x/y/z|the volume of the searching box|no default 

## Limitation
|Arguments| Description|Limitation
|--|--|--|
|--thread| the scale of parallelism (docking lanes)| preferably less than 10000
|--size_x/y/z|the volume of the searching box |less than 30/30/30