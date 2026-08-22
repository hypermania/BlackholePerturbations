##################################################################
##################################################################
# User settings: modify these for your use case
##################################################################
##################################################################
# Host compiler choice (needs support for C++20)
HOST_COMPILER ?= g++


##################################################################
# CUDA related settings

# Set if CUDA should be disabled or not (CUDA is enabled by default)
# You can also disable cuda by calling "make disable-cuda=true"
disable-cuda := false

# Location of the CUDA Toolkit
CUDA_PATH ?= /usr/local/cuda
# CUDA include path
CUDA_INCLUDE_DIR := "/usr/local/cuda/include/"
# CUDA library path
CUDA_LIBRARY_DIR := "/usr/local/cuda/lib64"

# Command for NVCC
NVCC := $(CUDA_PATH)/bin/nvcc -ccbin $(HOST_COMPILER)
NVCCFLAGS   := -m64 --threads 2

# Gencode arguments
# Use 86 for RTX 3060 Ti and 89 for RTX 4080. Change this for other GPUs / CUDA Toolkit version.
SMS ?= 89 # 50 52 60 61 70 75 80 86

ifeq ($(SMS),)
	$(info >>> WARNING - no SM architectures have been specified - waiving sample <<<)
	SAMPLE_ENABLED := 0
endif

ifeq ($(GENCODE_FLAGS),)
# Generate SASS code for each SM architecture listed in $(SMS)
$(foreach sm,$(SMS),$(eval GENCODE_FLAGS += -gencode arch=compute_$(sm),code=sm_$(sm)))
endif



##################################################################
##################################################################
# Non-user-settings: You probably won't need to change these.
##################################################################
##################################################################
# File names and file paths for the program
program_NAME := main
src_DIR := src
program_C_SRCS := $(wildcard $(src_DIR)/*.c)
program_CXX_SRCS := $(wildcard $(src_DIR)/*.cpp) $(wildcard $(src_DIR)/teukolsky_generated/*.cpp)
program_H_SRCS := $(wildcard $(src_DIR)/*.h)
program_HPP_SRCS := $(wildcard $(src_DIR)/*.hpp)
program_GEN_SRCS := $(wildcard $(src_DIR)/*.gen)
program_C_OBJS := ${program_C_SRCS:.c=.o}
program_CXX_OBJS := ${program_CXX_SRCS:.cpp=.o}
program_CXX_ASMS := ${program_CXX_SRCS:.cpp=.s}

program_OBJS := $(program_C_OBJS) $(program_CXX_OBJS)
program_INCLUDE_DIRS := "external"
program_LIBRARY_DIRS :=
program_LIBRARIES := m dl quadmath # fftw3 


# Names for CUDA C++ files
program_CU_SRCS := $(wildcard $(src_DIR)/*.cu) $(wildcard $(src_DIR)/*/*.cu)
program_CUH_SRCS := $(wildcard $(src_DIR)/*.cuh) $(wildcard $(src_DIR)/*/*.cuh)
program_CU_OBJS := ${program_CU_SRCS:.cu=.o}
device_link_OBJ := $(src_DIR)/device_link.o


# Option to disable CUDA
ifeq ($(disable-cuda),false)
	program_INCLUDE_DIRS += $(CUDA_INCLUDE_DIR)
	program_LIBRARY_DIRS += $(CUDA_LIBRARY_DIR)
	program_LIBRARIES += cudart_static cufft_static culibos cufile rt
	program_OBJS += $(program_CU_OBJS) $(device_link_OBJ)
else
	CXXFLAGS += -D DISABLE_CUDA
endif



# Compiler flags
CXXFLAGS += $(foreach includedir,$(program_INCLUDE_DIRS),-I$(includedir))
CXXFLAGS += -std=c++20 -Wall -DEIGEN_DONT_PARALLELIZE -DEIGEN_NO_CUDA -ftemplate-depth=20000
#-fext-numeric-literals  	#-DEIGEN_HAS_CONSTEXPR=1 #-DEIGEN_NO_DEBUG
CXXFLAGS += -march=alderlake -pthread -fopenmp
#CXXFLAGS += -march=native -pthread
CXXFLAGS += -O3 -ffast-math
#CXXFLAGS += -g -fno-omit-frame-pointer -fext-numeric-literals
CXXFLAGS += -DNDEBUG

NVCC_OPTIMIZE_FLAGS := --ftz=false
NVCC_INCLUDE_DIR_FLAGS += $(foreach includedir,$(program_INCLUDE_DIRS),-I$(includedir))
NVCCFLAGS += -std=c++20 -DCUDA_API_PER_THREAD_DEFAULT_STREAM -DEIGEN_NO_CUDA
NVCCFLAGS += $(foreach library,$(program_LIBRARIES),-l$(library))


# Add linker flags
LDFLAGS += $(foreach librarydir,$(program_LIBRARY_DIRS),-L$(librarydir)) 
LDLIBS += $(foreach library,$(program_LIBRARIES),-l$(library))


.PHONY: all clean distclean benchmark-precise benchmark-sds-precise \
	check-precise-performance check-precise-correctness \
	check-precise-sanitizers check-sds-precise-performance \
	check-sds-precise-correctness check-sds-precise-sanitizers

all: $(program_NAME)

benchmark-precise: test/benchmark_teukolsky_precise

benchmark-sds-precise: test/benchmark_sds_precise

check-precise-performance: benchmark-precise
	OMP_PROC_BIND=close OMP_PLACES=cores OMP_WAIT_POLICY=active \
		./test/benchmark_teukolsky_precise 50000 60 6

check-sds-precise-performance: benchmark-sds-precise
	OMP_PROC_BIND=close OMP_PLACES=cores OMP_WAIT_POLICY=active \
		./test/benchmark_sds_precise 50000 180 6

test/benchmark_teukolsky_precise: test/benchmark_teukolsky_precise.cpp src/teukolsky_precise.hpp
	$(HOST_COMPILER) -Iexternal -Isrc -std=c++20 -O3 -DNDEBUG -march=native \
		-fopenmp $< -lquadmath -o $@

test/benchmark_sds_precise: test/benchmark_sds_precise.cpp src/sds_precise.hpp \
		src/odeint_eigen/eigen_operations.hpp
	$(HOST_COMPILER) -Iexternal -Isrc -std=c++20 -O3 -DNDEBUG -march=native \
		-fopenmp $< -lquadmath -o $@

test/test_teukolsky_precise_correctness: test/test_teukolsky_precise_correctness.cpp \
		src/teukolsky_precise.hpp src/odeint_eigen/eigen_operations.hpp
	$(HOST_COMPILER) -Iexternal -Isrc -std=c++20 -O2 -Wall -Wextra \
		-fopenmp $< -lquadmath -o $@

test/test_eigen_scale_sums: test/test_eigen_scale_sums.cpp \
		src/odeint_eigen/eigen_operations.hpp
	$(HOST_COMPILER) -Iexternal -Isrc -std=c++20 -O2 -Wall -Wextra \
		-fopenmp $< -lquadmath -o $@

test/test_sds_precise_correctness: test/test_sds_precise_correctness.cpp \
		src/sds_precise.hpp src/odeint_eigen/eigen_operations.hpp
	$(HOST_COMPILER) -Iexternal -Isrc -std=c++20 -O2 -Wall -Wextra \
		-fopenmp $< -lquadmath -o $@

test/test_teukolsky_precise_correctness_prod: test/test_teukolsky_precise_correctness.cpp \
		src/teukolsky_precise.hpp src/odeint_eigen/eigen_operations.hpp
	$(HOST_COMPILER) -Iexternal -Isrc -std=c++20 -O3 -ffast-math -DNDEBUG \
		-Wall -Wextra -fopenmp $< -lquadmath -o $@

test/test_eigen_scale_sums_prod: test/test_eigen_scale_sums.cpp \
		src/odeint_eigen/eigen_operations.hpp
	$(HOST_COMPILER) -Iexternal -Isrc -std=c++20 -O3 -ffast-math -DNDEBUG \
		-Wall -Wextra -fopenmp $< -lquadmath -o $@

test/test_sds_precise_correctness_prod: test/test_sds_precise_correctness.cpp \
		src/sds_precise.hpp src/odeint_eigen/eigen_operations.hpp
	$(HOST_COMPILER) -Iexternal -Isrc -std=c++20 -O3 -ffast-math -DNDEBUG \
		-Wall -Wextra -fopenmp $< -lquadmath -o $@

check-precise-correctness: test/test_teukolsky_precise_correctness \
		test/test_eigen_scale_sums test/test_teukolsky_precise_correctness_prod \
		test/test_eigen_scale_sums_prod
	OMP_PROC_BIND=close OMP_PLACES=cores ./test/test_teukolsky_precise_correctness
	OMP_PROC_BIND=close OMP_PLACES=cores ./test/test_eigen_scale_sums
	OMP_PROC_BIND=close OMP_PLACES=cores ./test/test_teukolsky_precise_correctness_prod
	OMP_PROC_BIND=close OMP_PLACES=cores ./test/test_eigen_scale_sums_prod

check-sds-precise-correctness: test/test_sds_precise_correctness \
		test/test_sds_precise_correctness_prod
	OMP_PROC_BIND=close OMP_PLACES=cores ./test/test_sds_precise_correctness
	OMP_PROC_BIND=close OMP_PLACES=cores ./test/test_sds_precise_correctness_prod

check-precise-sanitizers:
	$(HOST_COMPILER) -Iexternal -Isrc -std=c++20 -O1 -g -Wall -Wextra \
		-fopenmp -fsanitize=address,undefined -fno-omit-frame-pointer \
		test/test_teukolsky_precise_correctness.cpp -lquadmath \
		-o test/test_teukolsky_precise_correctness_san
	$(HOST_COMPILER) -Iexternal -Isrc -std=c++20 -O1 -g -Wall -Wextra \
		-fopenmp -fsanitize=address,undefined -fno-omit-frame-pointer \
		test/test_eigen_scale_sums.cpp -lquadmath -o test/test_eigen_scale_sums_san
	OMP_NUM_THREADS=2 OMP_PROC_BIND=close OMP_PLACES=cores \
		./test/test_teukolsky_precise_correctness_san
	OMP_NUM_THREADS=2 OMP_PROC_BIND=close OMP_PLACES=cores \
		./test/test_eigen_scale_sums_san

check-sds-precise-sanitizers:
	$(HOST_COMPILER) -Iexternal -Isrc -std=c++20 -O1 -g -Wall -Wextra \
		-fopenmp -DSDS_SKIP_2000_DIGIT_REFERENCE \
		-fsanitize=address,undefined -fno-omit-frame-pointer \
		test/test_sds_precise_correctness.cpp -lquadmath \
		-o test/test_sds_precise_correctness_san
	OMP_NUM_THREADS=2 OMP_PROC_BIND=close OMP_PLACES=cores \
		./test/test_sds_precise_correctness_san

$(program_NAME): $(program_OBJS)
	$(LINK.cc) $(program_OBJS) -o $(program_NAME) $(LDLIBS)

$(program_OBJS): $(program_H_SRCS) $(program_HPP_SRCS) $(program_CUH_SRCS) $(program_GEN_SRCS)

%.o: %.cpp %.gen
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

%.o: %.cu
	$(NVCC) $(NVCC_INCLUDE_DIR_FLAGS) $(NVCCFLAGS) $(GENCODE_FLAGS) $(NVCC_OPTIMIZE_FLAGS) --diag-suppress 20012,20014 -o $@ -dc $<

$(device_link_OBJ): $(program_CU_OBJS)
	$(NVCC) $(NVCC_INCLUDE_DIR_FLAGS) $(NVCCFLAGS) $(GENCODE_FLAGS)  -o $@ --device-link $(program_CU_OBJS)

%.s: %.cpp
	$(CXX) $(CXXFLAGS) -S -fverbose-asm $< -o $@

asm: $(program_CXX_ASMS)

clean:
	$(RM) $(program_NAME)
	$(RM) test/benchmark_teukolsky_precise
	$(RM) test/benchmark_sds_precise
	$(RM) test/test_teukolsky_precise_correctness test/test_eigen_scale_sums
	$(RM) test/test_sds_precise_correctness test/test_sds_precise_correctness_prod
	$(RM) test/test_teukolsky_precise_correctness_prod test/test_eigen_scale_sums_prod
	$(RM) test/test_teukolsky_precise_correctness_san test/test_eigen_scale_sums_san
	$(RM) test/test_sds_precise_correctness_san
	$(RM) $(program_OBJS)
	$(RM) $(program_CXX_ASMS)
	$(RM) $(wildcard *~)
	$(RM) -r html latex

cudaclean:
	$(RM) $(program_CU_OBJS)
	$(RM) $(device_link_OBJ)

distclean: clean

show:
	echo $(CXX)
	echo $(GXX)
	echo $(GCC)
	echo $(LINK.cc)
	echo $(CC)
	echo $(CPP)
	echo $(RM)
	echo $(CXXFLAGS)
	echo $(NVCC)
	echo $(program_CXX_SRCS) "\n"
	echo $(program_HPP_SRCS) "\n"
	echo $(program_CXX_OBJS) "\n"
	echo $(program_OBJS) "\n"
	echo $(program_CU_SRCS) "\n"
	echo $(program_CUH_SRCS) "\n"
	echo $(program_CU_OBJS) "\n"
	echo $(device_link_OBJ) "\n"
	echo $(program_CXX_ASMS) "\n"
