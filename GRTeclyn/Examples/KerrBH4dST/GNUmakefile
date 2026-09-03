GRTECLYN_HOME ?= $(realpath ../../../../GRTeclyn)
GRFOLRES_HOME ?= $(realpath ../../../../GRFolres)

include $(GRTECLYN_HOME)/Tools/GNUMake/Make.defaults

EBASE = KerrBH4dST

# Line extraction uses the particle interpolator.
USE_PARTICLES = TRUE

AMREX_HOME ?= $(realpath ../../../../amrex)

include $(GRTECLYN_HOME)/Tools/GNUMake/Make.defs

include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/Boundary/Make.package
include $(AMREX_HOME)/Src/AmrCore/Make.package
include $(AMREX_HOME)/Src/Amr/Make.package
include $(AMREX_HOME)/Src/Particle/Make.package

GRTECLYN_SOURCE = $(GRTECLYN_HOME)/Source
GRFOLRES_SOURCE = $(GRFOLRES_HOME)/GRTeclyn/Source

src_dirs := $(GRTECLYN_SOURCE)/BlackHoles \
            $(GRTECLYN_SOURCE)/CCZ4 \
            $(GRTECLYN_SOURCE)/Grids \
            $(GRTECLYN_SOURCE)/GRTeclynCore \
            $(GRTECLYN_SOURCE)/IO \
            $(GRTECLYN_SOURCE)/Maths \
            $(GRTECLYN_SOURCE)/Matter \
            $(GRTECLYN_SOURCE)/ParticleInterpolator \
            $(GRTECLYN_SOURCE)/Tagging \
            $(GRFOLRES_SOURCE)/FourDerivScalarTensor \
            $(GRFOLRES_SOURCE)/ModifiedCCZ4 \
            $(realpath .)

# Find all source and header files in the src_dirs
CEXE_sources += $(notdir $(foreach dir, $(src_dirs), $(wildcard $(dir)/*.cpp)))
CEXE_headers += $(notdir $(foreach dir, $(src_dirs), $(wildcard $(dir)/*.hpp)))
CEXE_headers += $(notdir $(foreach dir, $(src_dirs), $(wildcard $(dir)/*.impl.hpp)))

VPATH_LOCATIONS += $(src_dirs)
INCLUDE_LOCATIONS += $(src_dirs)

include $(GRTECLYN_HOME)/Tools/GNUMake/Make.rules
