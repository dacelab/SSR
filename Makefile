## FORTRAN REGULAR COMPILER OPTIONS
FC = gfortran
FCFLAGS = -O3 -fopenmp  -fPIC -fmax-stack-var-size=2048 -fbounds-check# -O3

# Linking libraries 
LLFLAGS = -lblas -llapack -lm -lgomp
# LLFLAGS = -L/usr/lib/ -lblas -llapack 

PINC = -I/usr/include

## FORTRAN DEBUG COMPILER OPTIONS
D_FCFLAGS = -fopenmp -fPIC -fbounds-check -g  
D_FCFLAGS += -fmax-stack-var-size=2048

OBJ_SYNTH_LIST1 = lpp_mod.o Sobol_mod.o Utilities_mod.o minimiz_routines.o \
       lbfgs_routines.o glmnet_routines.o coordDesc_mod.o ALS_mod.o \
       crossValid_mod.o BSR_mod.o SRCD_mod.o sigmatune_mod.o main_SSRCD_synthetic.f90
       
OBJ_SYNTH_LIST2 = lpp_mod.o Sobol_mod.o Utilities_mod.o minimiz_routines.o \
       lbfgs_routines.o glmnet_routines.o coordDesc_mod.o ALS_mod.o \
       crossValid_mod.o BSR_mod.o SRCD_mod.o sigmatune_mod.o main_BSSR_synthetic.f90

OBJ_SYNTH_LIST3 = lpp_mod.o Sobol_mod.o Utilities_mod.o minimiz_routines.o \
       lbfgs_routines.o glmnet_routines.o coordDesc_mod.o ALS_mod.o \
       crossValid_mod.o BSR_mod.o SRCD_mod.o sigmatune_mod.o main_ALS_synthetic.f90

OBJ_REAL_LIST1 = lpp_mod.o Sobol_mod.o Utilities_mod.o minimiz_routines.o \
       lbfgs_routines.o glmnet_routines.o coordDesc_mod.o ALS_mod.o \
       crossValid_mod.o BSR_mod.o SRCD_mod.o sigmatune_mod.o main_SSRCD_real_splits.f90

OBJ_REAL_LIST2 = lpp_mod.o Sobol_mod.o Utilities_mod.o minimiz_routines.o \
       lbfgs_routines.o glmnet_routines.o coordDesc_mod.o ALS_mod.o \
       crossValid_mod.o BSR_mod.o SRCD_mod.o sigmatune_mod.o main_SSRCD_real_single.f90

lpp_mod.o: ./Common/lpp_mod.f90
	$(FC) $(FCFLAGS) -c ./Common/lpp_mod.f90 $(LLFLAGS)

Sobol_mod.o: ./Common/Sobol_mod.f90
	$(FC) $(FCFLAGS) -c ./Common/Sobol_mod.f90 $(LLFLAGS)

Utilities_mod.o: ./Common/Utilities_mod.f90
	$(FC) $(FCFLAGS) -c ./Common/Utilities_mod.f90 $(LLFLAGS)

# These routines MUST be compiled with optimization flags OFF:
minimiz_routines.o: ./Algos/minimiz_routines.f
	$(FC) -fPIC -c ./Algos/minimiz_routines.f

lbfgs_routines.o: ./Algos/lbfgs_routines.f
	$(FC) -fPIC -c ./Algos/lbfgs_routines.f

glmnet_routines.o: ./Algos/glmnet_routines.f
	$(FC) $(FCFLAGS) -c ./Algos/glmnet_routines.f $(LLFLAGS)

coordDesc_mod.o: ./Algos/coordDesc_mod.f90
	$(FC) $(FCFLAGS) -c ./Algos/coordDesc_mod.f90 $(LLFLAGS)

ALS_mod.o: ./Algos/ALS_mod.f90
	$(FC) $(FCFLAGS) -c ./Algos/ALS_mod.f90 $(LLFLAGS)

crossValid_mod.o: ./Algos/crossValid_mod.f90
	$(FC) $(FCFLAGS) -c ./Algos/crossValid_mod.f90 $(LLFLAGS)

BSR_mod.o: ./Algos/BSR_mod.f90
	$(FC) $(FCFLAGS) -c ./Algos/BSR_mod.f90 $(LLFLAGS)

SRCD_mod.o: ./Algos/SRCD_mod.f90
	$(FC) $(FCFLAGS) -c ./Algos/SRCD_mod.f90 $(LLFLAGS)

sigmatune_mod.o: ./Algos/sigmatune_mod.f90
	$(FC) $(FCFLAGS) -c ./Algos/sigmatune_mod.f90 $(LLFLAGS)

main_SSRCD_synthetic.o: main_SSRCD_synthetic.f90
	$(FC) $(FCFLAGS) -c main_SSRCD_synthetic.f90 $(LLFLAGS)

main_BSSR_synthetic.o: main_BSSR_synthetic.f90
	$(FC) $(FCFLAGS) -c main_BSSR_synthetic.f90 $(LLFLAGS)

main_ALS_synthetic.o: main_ALS_synthetic.f90
	$(FC) $(FCFLAGS) -c main_ALS_synthetic.f90 $(LLFLAGS)

main_SSRCD_real_splits.o: main_SSRCD_real_splits.f90
	$(FC) $(FCFLAGS) -c main_SSRCD_real_splits $(LLFLAGS)

main_SSRCD_real_single.o: main_SSRCD_real_single.f90
	$(FC) $(FCFLAGS) -c main_SSRCD_real_single $(LLFLAGS)

synth_ssrcd: $(OBJ_SYNTH_LIST1)
	$(FC) $(FCFLAGS) -o synth_ssrcd $(OBJ_SYNTH_LIST1) $(LLFLAGS)

synth_bssr: $(OBJ_SYNTH_LIST2)
	$(FC) $(FCFLAGS) -o synth_bssr $(OBJ_SYNTH_LIST2) $(LLFLAGS)

synth_als: $(OBJ_SYNTH_LIST3)
	$(FC) $(FCFLAGS) -o synth_als $(OBJ_SYNTH_LIST3) $(LLFLAGS)

real_ssrcd: $(OBJ_REAL_LIST1)
	$(FC) $(FCFLAGS) -o real_ssrcd $(OBJ_REAL_LIST1) $(LLFLAGS)

real_single_ssrcd: $(OBJ_REAL_LIST2)
	$(FC) $(FCFLAGS) -o real_single_ssrcd $(OBJ_REAL_LIST2) $(LLFLAGS)

clean:
	rm *.o *.mod