#!/bin/bash
# Script to make all the MINEOS executables and libraries

#======================  SET PATHS ======================#
source MINEOS_paths

# cd to source directory
cd src

#======================  COMPILE LIBRARIES  ======================#
# These are written out in full for a bit of extra transparency

# CIP
cd libcip
make -f MAKE_ciplib.mk
cd ..

# TAU
cd libtau
make
cd ..

# TIM
cd libtim
make
cd ..

# UTIL
cd libutil
make
cd ..

#======================  COMPILE MINEOS executables  ======================#

cd mineos
for imk in MAKE*
do 
	# make with an -i flag to ignore errors
	echo $imk
	rm *.o
	make -i -f $imk
done
cd ..

#plot_wk
cd plot_wk
make -f MAKE_plot_wk.mk 
cd ..

# cd back to parent directory
cd ..