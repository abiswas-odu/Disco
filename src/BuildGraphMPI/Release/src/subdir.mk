################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS_BUILDG_MPI += \
../../src/Dataset.cpp \
../../src/Edge.cpp \
../../src/HashTable.cpp \
../../src/OverlapGraph.cpp \
../../src/Read.cpp \
../../src/main.cpp 

OBJS_BUILDG_MPI += \
./src/Dataset.o \
./src/Edge.o \
./src/HashTable.o \
./src/OverlapGraph.o \
./src/Read.o \
./src/main.o 

CPP_DEPS_BUILDG_MPI += \
./src/Dataset.d \
./src/Edge.d \
./src/HashTable.d \
./src/OverlapGraph.d \
./src/Read.d \
./src/main.d 


# Each subdirectory must supply rules for building sources it contributes
src/BuildGraphMPI/src/%.o: ../../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	-$(CC) -O3 -g3 -Wall -fopenmp -lgomp -std=c++11 -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


