################################# Explicit ####################################


# # Compiler
# CXX = g++

# # Compiler flags
# CXXFLAGS = -Wall -g

# # Source files
# SRC = main.cpp fibo.cpp

# # Object files
# OBJ = main.o fibo.o

# # Executable file
# TARGET = fibo

# # Default target
# all: $(TARGET)

# # Rule to build the target
# $(TARGET): $(OBJ)
# 	$(CXX) $(OBJ) -o $(TARGET)

# # Rule to build object files from source files
# main.o: main.cpp fibo.h
# 	$(CXX) $(CXXFLAGS) -c main.cpp

# fibo.o: fibo.cpp fibo.h
# 	$(CXX) $(CXXFLAGS) -c fibo.cpp

# # Clean up object files and executable
# clean:
# 	rm -f $(OBJ) $(TARGET)

################################# Implicit ####################################

# Makefile with implicit rules

# Compiler
CXX = g++

# Compiler flags
CXXFLAGS = -Wall -g

# Source files
SRC = main.cpp fibo.cpp

# Object files
OBJ = main.o fibo.o

# Executable file
TARGET = fibo

# Default target
all: $(TARGET)

# Link object files to create executable
$(TARGET): $(OBJ)
	$(CXX) $(OBJ) -o $(TARGET)

# Clean up object files and executable
clean:
	rm -f $(OBJ) $(TARGET)
