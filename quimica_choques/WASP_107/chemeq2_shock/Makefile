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

SRC_OMP = $(SRC_DIR)/network.f90 \
          $(SRC_DIR)/solver.f90 \
          $(SRC_DIR)/test_1D_multimodel_omp.f90

EXE     = 1D_test.exe
EXE_OMP = 1D_test_multi_omp.exe

# ==========================
# Default target
# ==========================
all: $(EXE)

# ==========================
# Build executable (original)
# ==========================
$(EXE): $(SRC)
	mkdir -p $(MODDIR)
	$(FC) -J$(MODDIR) $(FFLAGS) $(SRC) -o $(EXE)

# ==========================
# Build executable (OpenMP)
# ==========================
$(EXE_OMP): $(SRC_OMP)
	mkdir -p $(MODDIR)
	$(FC) -J$(MODDIR) $(FFLAGS) -fopenmp $(SRC_OMP) -o $(EXE_OMP)

omp: $(EXE_OMP)

# ==========================
# Clean
# ==========================
clean:
	rm -rf $(MODDIR) *.o *.mod $(EXE) $(EXE_OMP)

# ==========================
# Rebuild from scratch
# ==========================
rebuild: clean all
rebuild_omp: clean omp