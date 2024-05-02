#
# Makefile for DOM examples
#
default: all
 all:  phis2e
#
#---------------------------
#MK=$(FLIB_ROOT)#/fortran.mk
#echo '$(MK)'
#include $(MK)
MK=/home/nico/Workspace/Progs/xmlf90-1.2g/macros/fortran.mk
#MK=$(FLIB_ROOT)/fortran.mk
include $(MK)
#---------------------------
#
# Uncomment the following line for debugging support
#
FFLAGS=-ffree-line-length-none
#FFLAGS=$(FFLAGS_CHECK)
#
#LIBS=$(LIB_PREFIX)$(LIB_STD)     -lflib  -I/opt/intel/mklinclude  -L/opt/intel/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core
LIBS=$(LIB_PREFIX)$(LIB_STD)     -lflib -mkl
#LIBS=$(LIB_PREFIX)$(LIB_STD) -lflib -L/usr/lib64/atlas/ -llapack -lblas -Bstatic_pgi -Bstatic 
##LIBS=$(LIB_PREFIX)$(LIB_STD) -lflib -L/opt/intel/mkl/lib/intel64 -I/opt/intel/mkl/lib/intel64  -llapack

#LDFLAGS = -p -g -pg
#LDFLAGS = -O3 -parallel
#LDFLAGS =  -qopenmp -O3 #-xHost  -no-prec-div 
#FC = fortran-intel
#testov: general.o overlap.o testOv.o       
#	$(FC) $(LDFLAGS) -o testOv general.o overlap.o testOv.o $(LIBS)

#testpot: general.o f11.o nuclear.o special_functions.o  testNuclear.o
#	$(FC) $(LDFLAGS) -o testPot general.o f11.o nuclear.o special_functions.o  testNuclear.o $(LIBS)
#	$(FC) $(LDFLAGS) -o testPot general.o f11.o nuclear.o special_functions.o  testNuclear.o $(LIBS)

phis2e:  general.o setup.o newtypes.o centerlib.o linearinterp.o  splineinterp.o tools.o  matrices.o inputlib.o fdn.o colllib_2e.o cgto.o  overlap.o special_functions.o f11.o nuclear.o repulsion.o matrep.o  onee_int_cgto.o twoe_int_cgto.o phis_twoe_int.o
	$(FC) $(LDFLAGS) -o   phis2e general.o setup.o newtypes.o centerlib.o linearinterp.o  splineinterp.o  tools.o  matrices.o inputlib.o fdn.o colllib_2e.o cgto.o  overlap.o special_functions.o f11.o nuclear.o repulsion.o matrep.o  onee_int_cgto.o  twoe_int_cgto.o phis_twoe_int.o $(LIBS) 
#phis2e:  general.o setup.o newtypes.o centerlib.o linearinterp.o  splineinterp.o tools.o  matrices.o inputlib.o fdn.o colllib_2e.o cgto.o  overlap.o special_functions.o f11.o nuclear.o repulsion.o matrep.o collint.o collint_2e.o onee_int_cgto.o twoe_int_cgto.o phis_twoe_int.o
#	$(FC) $(LDFLAGS) -o   phis2e general.o setup.o newtypes.o centerlib.o linearinterp.o  splineinterp.o  tools.o  matrices.o inputlib.o fdn.o colllib_2e.o cgto.o  overlap.o special_functions.o f11.o nuclear.o repulsion.o matrep.o collint.o collint_2e.o onee_int_cgto.o  twoe_int_cgto.o phis_twoe_int.o $(LIBS) 
### to remove the useless modules for collint.o collint_2e.o.....

#
clean: 
	rm -f  *.o *.$(MOD_EXT)

