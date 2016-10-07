################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/Config.cpp \
../src/DataSet.cpp \
../src/Edge.cpp \
../src/OverlapGraph.cpp \
../src/Read.cpp \
../src/Utils.cpp \
../src/dna.cpp \
../src/main.cpp 

OBJS += \
./src/Config.o \
./src/DataSet.o \
./src/Edge.o \
./src/OverlapGraph.o \
./src/Read.o \
./src/Utils.o \
./src/dna.o \
./src/main.o 

CPP_DEPS += \
./src/Config.d \
./src/DataSet.d \
./src/Edge.d \
./src/OverlapGraph.d \
./src/Read.d \
./src/Utils.d \
./src/dna.d \
./src/main.d 

ifeq ($(READGZ), 1)  #at this point, the makefile checks if FEATURE is enabled
OPTS = -DINCLUDE_READGZ #variable passed to g++
GZLIB = -lz 
endif

# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ $(OPTS) -O3 -Wall -c -fmessage-length=0 -fopenmp -Wno-sign-compare $(GZLIB) -lgomp -std=c++11 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


