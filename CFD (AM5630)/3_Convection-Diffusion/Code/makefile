# ----------------------------------------------------------------------
# FORTRAN Makefile for the 2D Convection Diffusion Solver
# ----------------------------------------------------------------------
# Author: Kanak Agarwal (kanakaero.github.io)
# Last Modified: 12th November 2025
# ----------------------------------------------------------------------

# Compiler
FC = gfortran
FCFLAGS = -O2 -Wall

# Executable
EXEC = conv_diff

# Source files
SRC = main.f90 \
      var_dec.f90 \
      gridcon.f90 \
	  discret.f90 \
	  bc.f90 \
	  gauss_seidel.f90 \
	  thomas.f90 \
	  tdma.f90 \
	  postproc.f90 \
	  debug.f90 \
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
	rm -f *.tec
