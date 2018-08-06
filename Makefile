# C++ compiler to use
CXX = g++

# Compiler flags
CXXFLAGS = -std=c++11

# Add options to generate dependencies
OUTPUT_OPTION=-MMD -MP -o $@

# Check for "debug=true" or "debug=1" flag, modify CFLAGS accordingly
ifneq (,$(filter $(debug),true 1))
   CXXFLAGS += -Wall -g -DDEBUG
else
   CXXFLAGS += -O3
endif

# Check for "omp=true" or "omp=1" flag, modify CFLAGS accordingly
ifneq (,$(filter $(omp),true 1))
   CXXFLAGS += -DOPENMP
endif

# Check for "exceptions=true" or "exceptions=1" flag, modify CFLAGS accordingly
ifneq (,$(filter $(exceptions),true 1))
   CXXFLAGS += -DUSE_THREAD_EXCEPTIONS
endif

# Libraries
LIBS = 

# Sources
SRCS = main.cpp\
       conversions.cpp\
       div_b_fix.cpp\
       file_access.cpp\
       vector_utilities.cpp\
       enums.cpp\
       scheme_central_fd.cpp\
       scheme_eno.cpp\
       time_steppers.cpp\
       spatialintegrationmethod.cpp\
       spatialmethodcentralfd2.cpp\
       spatialmethodeno.cpp\
       spatialmethodenoroe.cpp\
       spatialmethodenolf.cpp\
       timeintegrationmethod.cpp\
       timeintegrationeuler.cpp\
       timeintegrationrk3.cpp

# Object files
OBJS = $(SRCS:.cpp=.o)

# Dependencies
DEPS = $(SRCS:.cpp=.d)

# Default values, to be overwritten from Makefile.local
INCLUDES = 
LFLAGS = 
MAIN = mhd2d

# Read definitions from Makefile.local, if it exists
-include Makefile.local

.PHONY: clean

default: all

# Main
all: $(MAIN)
	@echo Compilation finished.

-include $(DEPS)

# Link object files into executable
$(MAIN): $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $(MAIN) $(OBJS) $(LFLAGS) $(LIBS)

# Compile into object files
.cpp.o:
	$(CXX) $(CXXFLAGS) $(INCLUDES) -MMD -c $< -o $@

# Clean executable, object files and temporary files
clean:
	$(RM) $(OBJS) $(DEPS) *~ $(MAIN)

