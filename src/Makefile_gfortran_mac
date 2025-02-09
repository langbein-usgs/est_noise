
FC = gfortran

# useful for debugging
FFLAGS = -Ofast  -w -std=f2008 -lm  -g -fcheck=all  #  -mcmodel=medium

FFLAGS1 = -Ofast -w -std=legacy
#  fft --- https://www.netlib.org/fftpack/
#          https://fortran-lang.github.io/fftpack/ 
FFT = dfft.o

##  using apple'version of blas --- seems to work best with apple cpu (M1)
BLAS = -fexternal-blas -framework accelerate

# homebrew   similar run time as apple's but seems to use fewer resources
#BLAS = -I/usr/local/opt/openblas/include -lopenblas -L/usr/local/opt/openblas/lib
#BLAS = -I/opt/homebrew/opt/openblas/include -lopenblas -L/opt/homebrew/opt/openblas/lib

# macports  runtime is about 30% slower than either apple's version or homebrew's version
#BLAS= -I/opt/local/include -lopenblas -L/opt/local/lib

SRC = time_mod.f90  alloc_mod.f90 est_noise_mod.f90 filterfunc_mod.f90  psd_calc8.f90 bust_mod.f90  GetStuff_mod.f90  bust8.f90 \
  adjust8.f90 genNoise_mod.f90  randnor.f gen_noise8.f90  mle_mod.f90   NedlerMeadSimplex_mod.f90  model_fit_A_mod.f90 est_noise8.01.f90 \
  compare_wander_mod.f90 compare_wander8.01.f90


%.o: %.f90  #wildcard rule, creation of *.o depends on *.f90
	$(FC) $(FFLAGS)  -c $<

dfft:
	gfortran -c dfft.f -std=legacy -w
	
OBJ1 = time_mod.o alloc_mod.o filterfunc_mod.o  psd_calc8.o

#  need to include the fft library 	
psd_calc8: $(OBJ1)
	$(FC)  $(FFT) $(BLAS) -o $@ $(OBJ1) 

OBJ2 = bust_mod.o time_mod.o GetStuff_mod.o  bust8.o    
bust8: $(OBJ2)
	$(FC) -o $@ $(OBJ2)
	
OBJ3 = GetStuff_mod.o time_mod.o adjust8.o
adjust8: $(OBJ3)
	$(FC) -o $@ $(OBJ3)
	
OBJ4 = alloc_mod.o genNoise_mod.o time_mod.o randnor.o filterfunc_mod.o  gen_noise8.o
gen_noise8: $(OBJ4)
	$(FC) $(FFT) -o $@ $(OBJ4)

OBJ5 =   alloc_mod.o  est_noise_mod.o model_fit_A_mod.o  mle_mod.o filterfunc_mod.o  genNoise_mod.o randnor.o GetStuff_mod.o \
        time_mod.o  NedlerMeadSimplex_mod.o est_noise8.01.o	
est_noise8: $(OBJ5)
	$(FC) $(BLAS) $(FFT) -o $@ $(OBJ5)
	
OBJ6= time_mod.o alloc_mod.o GetStuff_mod.o compare_wander_mod.o randnor.o filterfunc_mod.o genNoise_mod.o compare_wander8.01.o
compare_wander8: $(OBJ6)
	$(FC) $(FFT) -o $@ $(OBJ6)
	
all:  dfft psd_calc8 bust8 adjust8 gen_noise8 est_noise8 compare_wander8
	@echo "Ready to run make install"
	
install:
	mv psd_calc8 bust8 adjust8 gen_noise8 est_noise8 compare_wander8 ../bin
	@echo "Ready to run make clean"
clean:
	rm -fr *.o *.mod *.dSYM core psd_calc8 bust8 adjust8 gen_noise8 est_noise8 compare_wander8 *.out *.jrn


###

