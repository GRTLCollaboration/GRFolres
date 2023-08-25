# -*- Mode: Makefile -*-

### This makefile produces an executable for each name in the `ebase'
###  variable using the libraries named in the `LibNames' variable.

# Included makefiles need an absolute path to the Chombo installation
# CHOMBO_HOME := Please set the CHOMBO_HOME locally (e.g. export CHOMBO_HOME=... in bash)

# User needs to specify GRCHOMBO_SOURCE location, e.g. using
# export GRCHOMBO_SOURCE = ~/GRChombo/Source
export MODIFIEDGRAVITY_SOURCE = $(shell pwd)/../../Source

ebase := Main_BinaryBHCubicHorndeski

LibNames := AMRTimeDependent AMRTools BoxTools

src_dirs := $(GRCHOMBO_SOURCE)/utils \
            $(GRCHOMBO_SOURCE)/simd  \
            $(MODIFIEDGRAVITY_SOURCE)/ModifiedCCZ4  \
	    $(MODIFIEDGRAVITY_SOURCE)/CubicHorndeski  \
            $(MODIFIEDGRAVITY_SOURCE)/TaggingCriteria  \
	    $(GRCHOMBO_SOURCE)/CCZ4  \
            $(GRCHOMBO_SOURCE)/BoxUtils  \
            $(GRCHOMBO_SOURCE)/GRChomboCore  \
            $(GRCHOMBO_SOURCE)/AMRInterpolator  \
            $(GRCHOMBO_SOURCE)/TaggingCriteria  \
	    $(GRCHOMBO_SOURCE)/ApparentHorizonFinder \
            $(GRCHOMBO_SOURCE)/InitialConditions/BlackHoles \
            $(GRCHOMBO_SOURCE)/BlackHoles

include $(CHOMBO_HOME)/mk/Make.test
