#include /save/alegarra/f90-05/Makeinit
#include /save/alegarra/f90-05/Makeinit_debug
#include ./Makeinit_debug
include ./Makeinit
prog=vrr

dir=./libs

$(prog):	$(prog).o  $(dir)/sparsem.a 
	$(f90) $(optf90) $(prog).o   \
	$(dir)/sparsem.a  -o $(prog)
	
$(prog).o:	$(prog).f90  $(dir)/sparsem.a 
	$(f90) $(optdir)$(dir) -c $(optf90) $(prog).f90 

	
