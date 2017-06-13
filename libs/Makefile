# Makes libraries

#include ../Makeinit_debug
include ../Makeinit

all:  sparsem.a 

sparsem.a: kind.o hash.o sparse_reduced.o sparse2.o \
	 sparssub.o second.o 
	ar cr sparsem.a kind.o hash.o sparse_reduced.o sparse2.o  \
            sparssub.o second.o

kind.o:	kind.f90
	$(f90) -c  $(optf90) kind.f90

hash.o:	hash.f90
	$(f90) -c $(optf90)  hash.f90

sparse_reduced.o:	sparse_reduced.f90 
	$(f90) -c  $(optf90) sparse_reduced.f90

sparse2.o:	sparse2.f 
	$(f90) -c  $(optf77) sparse2.f

sparssub.o:	sparssub.f 
	$(f90) $(optf77) -c sparssub.f

second.o:	second.f
	$(f90) -c  $(optf77) second.f

clean:
	rm *.o *.a *.mod

	
windows:
	cat kind.f90  lapack90r.f90 denseop.f90 sparse.f90 fspak90.f90 \
	    pcg.f90 model.f90 blups1.f90 prob.f90 gibbs.f90 >../windows/bglibs90.f90
	cat sparse2.f fspak.f fspaksub.f ginv.f sparssub.f ranlib.f \
			blups2.f second.f > ../windows/bglibs77.f	
