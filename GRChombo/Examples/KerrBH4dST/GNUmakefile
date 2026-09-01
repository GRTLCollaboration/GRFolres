# -*- Mode: Makefile -*-

### This makefile produces an executable for each name in the `ebase'
###  variable using the libraries named in the `LibNames' variable.

# Included makefiles need an absolute path to the Chombo installation
# CHOMBO_HOME := Please set the CHOMBO_HOME locally (e.g. export CHOMBO_HOME=... in bash)

# User needs to specify GRCHOMBO_SOURCE location, e.g. using
# export GRCHOMBO_SOURCE = ~/GRChombo/Source
export GRFOLRES_SOURCE = $(shell pwd)/../../Source

ebase := Main_KerrBH4dST

LibNames := AMRTimeDependent AMRTools BoxTools

src_dirs := $(GRCHOMBO_SOURCE)/utils \
            $(GRCHOMBO_SOURCE)/simd  \
            $(GRCHOMBO_SOURCE)/CCZ4  \
            $(GRFOLRES_SOURCE)/ModifiedCCZ4 \
	    $(GRFOLRES_SOURCE)/FourDerivScalarTensor  \
            $(GRCHOMBO_SOURCE)/BoxUtils  \
            $(GRCHOMBO_SOURCE)/GRChomboCore  \
            $(GRCHOMBO_SOURCE)/TaggingCriteria  \
            $(GRCHOMBO_SOURCE)/InitialConditions/BlackHoles \
            $(GRCHOMBO_SOURCE)/BlackHoles \
            $(GRCHOMBO_SOURCE)/AMRInterpolator \
	    $(GRCHOMBO_SOURCE)/ApparentHorizonFinder	    

include $(CHOMBO_HOME)/mk/Make.test
