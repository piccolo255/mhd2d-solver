# Makefile modeled after the example at
#   http://www.cs.swarthmore.edu/~newhall/unixhelp/howto_makefiles.html

# C++ compiler to use
CPP = g++

# Compiler flags
CPPFLAGS = -std=c++11

# Check for "debug=true" flag, modify CFLAGS accordingly
ifeq ($(debug),true)
	CPPFLAGS += -Wall -pedantic -g -DDEBUG
else
	CPPFLAGS += -Wall -O3
endif

# Libraries
LIBS = 

# Sources
SRCS =  main.cpp conversions.cpp file_access.cpp riemann.cpp
SRCS += scheme_central_fd.cpp scheme_eno.cpp time_steppers.cpp vector_utilities.cpp

# Object files
OBJS = $(SRCS:.cpp=.o)

# Default values, to be overwritten from Makefile.local
INCLUDES = 
LFLAGS = 
MAIN = euler1d

# Read definitions from Makefile.local, if it exists
-include Makefile.local

# .PHONY: clean

# Main
all: $(MAIN)
	@echo Compilation finished.

# Link object files into executable
$(MAIN): $(OBJS) 
	$(CPP) $(CPPFLAGS) $(INCLUDES) -o $(MAIN) $(OBJS) $(LFLAGS) $(LIBS)

# Compile into object files
.cpp.o:
	$(CPP) $(CPPFLAGS) $(INCLUDES) -c $<  -o $@

# Clean executable, object files and temporary files
clean:
	$(RM) *.o *~ $(MAIN)
