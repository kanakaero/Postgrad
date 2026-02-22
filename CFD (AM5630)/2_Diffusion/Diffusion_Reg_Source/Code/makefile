# ----------------------------------------------------------------------
# FORTRAN Makefile for the 2D Diffusion Solver
# ----------------------------------------------------------------------
# Author: Kanak Agarwal (kanakaero.github.io)
# Last Modified: 25th September 2025
# ----------------------------------------------------------------------

# Compiler
FC = gfortran
FCFLAGS = -O2 -Wall

# Executable
EXEC = 2d_diff

# Source files
SRC = main.f90 \
      var_dec.f90 \
      mesh.f90 \
      fparser.f90 \
      discret.f90 \
      solver.f90 \
      bc_s_k.f90 \
      dealloc.f90

# Object files
OBJ = $(SRC:.f90=.o)

# Module files
MOD = $(SRC:.f90=.mod)

# Default target
all: $(EXEC)

# Link the executable
$(EXEC): $(OBJ)
	$(FC) $(FCFLAGS) -o $(EXEC) $(OBJ)

# Compile source files into objects
%.o: %.f90
	$(FC) $(FCFLAGS) -c $<

# Clean files
clean:
	rm -f $(OBJ) $(EXEC)
	rm -f *.dat
