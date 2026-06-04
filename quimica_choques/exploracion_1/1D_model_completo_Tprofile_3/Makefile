# ==========================
# Compiler
# ==========================
FC      = gfortran
FFLAGS  = -O2 -Wall
MODDIR  = mod

# ==========================
# Directories
# ==========================
SRC_DIR = src

# ==========================
# Files
# ==========================
SRC = $(SRC_DIR)/network.f90 \
      $(SRC_DIR)/solver.f90 \
      $(SRC_DIR)/test_1D_1.f90

EXE = 1D_test.exe

# ==========================
# Default target
# ==========================
all: $(EXE)

# ==========================
# Build executable
# ==========================
$(EXE): $(SRC)
	mkdir -p $(MODDIR)
	$(FC) -J$(MODDIR) $(FFLAGS) $(SRC) -o $(EXE)

# ==========================
# Clean
# ==========================
clean:
	rm -rf $(MODDIR) *.o *.mod $(EXE)

# ==========================
# Rebuild from scratch
# ==========================
rebuild: clean all