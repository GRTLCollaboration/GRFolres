# GRTeclyn tests for the GRFolres modified-gravity source.
#
# Layout (grteclyn branch): the GRTeclyn framework lives in the sibling
#   ../../../GRTeclyn   checkout, the modified-gravity source in ../Source .
# This makefile mirrors GRTeclyn/Tests/GNUmakefile, adding the GRFolres source
# directories and the local test/Common directories.

GRTECLYN_HOME ?= $(realpath ../../../GRTeclyn)

include $(GRTECLYN_HOME)/Tools/GNUMake/Make.defaults

# We don't want to use MPI by default for the tests
USE_MPI       = FALSE
USE_PARTICLES = TRUE

EBASE = Tests

AMREX_HOME    ?= $(realpath ../../../amrex)
GRFOLRES_HOME ?= $(realpath ../..)

include $(GRTECLYN_HOME)/Tools/GNUMake/Make.defs

include $(AMREX_HOME)/Src/Amr/Make.package
include $(AMREX_HOME)/Src/AmrCore/Make.package
include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/Boundary/Make.package
include $(AMREX_HOME)/Src/Particle/Make.package

GRTECLYN_SOURCE     = $(GRTECLYN_HOME)/Source
GRFOLRES_SOURCE     = $(GRFOLRES_HOME)/GRTeclyn/Source
GRTECLYN_TESTS_HOME = $(realpath .)

src_dirs := $(GRTECLYN_SOURCE)/ParticleInterpolator \
            $(GRTECLYN_SOURCE)/BlackHoles \
            $(GRTECLYN_SOURCE)/CCZ4 \
            $(GRTECLYN_SOURCE)/Grids \
            $(GRTECLYN_SOURCE)/GRTeclynCore \
            $(GRTECLYN_SOURCE)/IO \
            $(GRTECLYN_SOURCE)/Maths \
            $(GRTECLYN_SOURCE)/Matter \
            $(GRTECLYN_SOURCE)/Tagging \
            $(GRFOLRES_SOURCE)/FourDerivScalarTensor \
            $(GRFOLRES_SOURCE)/ModifiedCCZ4

# Find all directories that end in "Test" or "Tests"
# (hidden directories are not included)
test_dirs_match := $(wildcard $(GRTECLYN_TESTS_HOME)/*Test/ \
                              $(GRTECLYN_TESTS_HOME)/*Tests/)
test_dirs       := $(patsubst %/, %, $(test_dirs_match))

# Finally add in the Common directory
test_dirs += $(GRTECLYN_TESTS_HOME)/Common

search_dirs := $(realpath .) $(src_dirs) $(test_dirs)

doctest_dir     := $(GRTECLYN_HOME)/external/doctest
doctest_headers := doctest.h

# Find all source and header files
CEXE_sources += $(notdir $(foreach dir, $(search_dirs), $(wildcard $(dir)/*.cpp)))
CEXE_headers += $(notdir $(foreach dir, $(search_dirs), $(wildcard $(dir)/*.hpp)))
CEXE_headers += $(notdir $(foreach dir, $(search_dirs), $(wildcard $(dir)/*.impl.hpp)))

# Add doctest header
CEXE_headers += $(doctest_headers)

VPATH_LOCATIONS          += $(src_dirs) $(test_dirs)
INCLUDE_LOCATIONS        += $(realpath .) $(src_dirs) $(test_dirs)
SYSTEM_INCLUDE_LOCATIONS += $(doctest_dir)

include $(GRTECLYN_HOME)/Tools/GNUMake/Make.rules

ifeq ($(USE_MPI), TRUE)
  ifeq ($(origin SLURM_NTASKS),undefined)
    LAUNCH = mpiexec -n 2 ./
  else
    LAUNCH = mpiexec -n ${SLURM_NTASKS} ./
  endif
else
  LAUNCH = ./
endif

DOCTEST_OPTIONS := -dt-d=yes

run: $(executable)
	$(info Running tests for configuration: $(optionsSuffix))
	$(SILENT) $(LAUNCH)$(executable) $(DOCTEST_OPTIONS)
