# -*- Mode: Makefile -*-

### This makefile produces an executable for each name in the `ebase'
###  variable using the libraries named in the `LibNames' variable.

# Included makefiles need an absolute path to the Chombo installation
# CHOMBO_HOME := Please set the CHOMBO_HOME locally(e.g. export CHOMBO_HOME=... in bash)

# User needs to specify GRCHOMBO_SOURCE location, e.g. using
# export GRCHOMBO_SOURCE = ~/GRChombo/Source
export MODIFIEDGRAVITY_SOURCE = $(shell pwd)/../../Source

ebase := TestField4dST_Test

LibNames := BoxTools

src_dirs := $(GRCHOMBO_SOURCE)/utils  \
	    $(GRCHOMBO_SOURCE)/simd  \
	    $(MODIFIEDGRAVITY_SOURCE)/ModifiedCCZ4  \
	    $(MODIFIEDGRAVITY_SOURCE)/TestField4dST  \
	    $(GRCHOMBO_SOURCE)/CCZ4  \
	    $(GRCHOMBO_SOURCE)/BoxUtils

include $(CHOMBO_HOME)/mk/Make.test
