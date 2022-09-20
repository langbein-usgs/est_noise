#!/bin/sh

#  script used to compile est_noise
#  modify the gfortran call for mac to include "legacy" to compile the dfft stuff  (-std=legacy)
#   May need to do this for linux, too
if [ "$#" -ne 2 ]
then
  echo Script used to compile est_noise
  echo "Usage:  compile.sh mac/linux i/g "
  echo "  where argument one designates the OS, mac or linux"
  echo "    and argument two designates the compiler"
  echo "    i -- intel ifort"
  echo "    g -- gfortran"
  exit
fi

ver=7.30

file=`echo "-o bin/est_noise"$ver" ./src/est_noise"$ver".f ./src/time.f  ./src/funmin"$ver".f ./src/dfft.f ./src/NedlerMeadSimplex"$ver".f  ./src/modify.f ./src/NedlerMeadSimplexEXP"$ver".f" ./src/filterfunc.f ./src/genNoise.f src/randnor.f ./src/GetData.f  `
#file=`echo "-o bin/est_noise"$ver" ./src/est_noise"$ver".f ./src/time.f  ./src/funmin"$ver"_test2.f ./src/dfft.f ./src/NedlerMeadSimplex"$ver".f  ./src/modify.f ./src/NedlerMeadSimplexEXP"$ver".f" ./src/filterfunc.f ./src/genNoise.f src/randnor.f ./src/GetData.f  `

echo $file

filepsd=`echo -o ./bin/psd_calc7 ./src/psd_calc7.f ./src/filterfunc.f ./src/dfft.f`
 
filegen=`echo -o ./bin/gen_noise7.02 ./src/gen_noise7.02.f ./src/genNoise.f ./src/filterfunc.f ./src/dfft.f ./src/randnor.f ./src/time.f `

filebust=`echo -o ./bin/bust_5 ./src/bust_5.f ./src/time.f ./src/GetData.f `

filewand=`echo -o bin/compare_wander7.01 src/compare_wander7.01.f src/time.f src/GetData.f src/filterfunc.f src/dfft.f src/randnor.f `

fileadj=`echo -o bin/adjust_1 src/adjust_1.f src/time.f src/GetData.f `
rm -f bin/est_noise"$ver" bin/psd_calc7 bin/gen_noise7.02 bin/adjust* bin/compare_wander*

if [ "$1" = "linux" ]
then
  echo OS is $1
  if [ "$2" = "g" ]
  then
    echo compiler is gfortran

     echo compile est_noise7
###########  UNCOMMENT ONE OF THE TW NEXT LINES
#     gfortran $file -L/usr/lib -lopenblas -mcmodel=medium -O3
#  this seems to work for CentOS
     gfortran $file -L/usr/lib64 -lopenblas -mcmodel=medium -O3 -w -std=legacy   ## CentOS
#  This seems to work for ubuntu
#     gfortran $file -static -L/usr/lib -lopenblas -mcmodel=medium -O3 -lpthread   ## Ubuntu
################################
#     gfortran $file lib/libopenblas.so -mcmodel=medium -O3   ####### Not used!
     echo compile psd_calc7
    gfortran $filepsd -O3 -mcmodel=medium -w -std=legacy
     echo compile gen_noise
    gfortran $filegen -O3 -mcmodel=medium -w -std=legacy
     echo compile compare_wander
    gfortran $filewand -O3 -mcmodel=medium -w -std=legacy
     echo compile bust_5
    gfortran $filebust -O3 -w -std=legacy
     echo compile adjust
    gfortran $fileadj -O3 -w -std=legacy
  elif [ "$2" = "i" ]
  then
    echo compiler is intel ifort
     echo compile est_noise7
    MKLPATH=/opt/intel/mkl/lib/intel64
    MKLINCLUDE=/opt/intel/mkl/include/intel64
    ifort $file -mcmodel medium -shared-intel   -warn none -heap-arrays 2048  -L$MKLPATH -I$MKLINCLUDE -I$MKLINCLUDE/lp64 -lmkl_lapack95_lp64 \
  -Wl,--start-group $MKLPATH/libmkl_intel_lp64.a $MKLPATH/libmkl_intel_thread.a $MKLPATH/libmkl_core.a \
   /opt/intel/lib/intel64/libiomp5.a   -lpthread -Wl,--end-group
     echo compile psd_calc7
    ifort $filepsd -O3 -mcmodel medium 
     echo compile gen_noise
    ifort $filegen -O3 -mcmodel medium 
     echo compile compare_wander
    ifort $filewand -O3 -mcmodel medium 
     echo compile bust_5
    ifort $filebust -O3
     echo compile adjust
    ifort $fileadj -O3
  else
    echo "Invalid compiler option"
    exit
  fi

elif  [ "$1" = "mac" ]
then
  echo OS is $1
  if [ "$2" = "g" ]
  then
    echo compiler is gfortran
     echo compile est_noise
#     gfortran $file -O3 -framework Accelerate -mcmodel=medium  -Wl,-stack_size,0x80000000
#  currently works 
#     gfortran $file -O3 -framework Accelerate -mcmodel=medium  -Wa,-q -w -std=legacy
#  use with openblas from homebrew
     gfortran $file -I/usr/include/openblas/ -I/usr/local/Cellar/openblas/0.3.21/include/  -lopenblas -mcmodel=medium -O3 -w -std=legacy  -L/usr//local/Cellar/openblas/0.3.21/lib/  
#  use with valgrind
#     gfortran $file -O0 -g -framework Accelerate  -mcmodel=medium -Wl,-stack_size,0x80000000   #  -Wa,-q
#     exit
#    gfortran $file -O3 -framework Accelerate -mdynamic-no-pic
     echo compile psd_calc7
    gfortran $filepsd -O3  -Wa,-q -w -std=legacy
     echo compile gen_noise
    gfortran $filegen -O3  -Wa,-q -w -std=legacy
     echo compile compare_wander
    gfortran $filewand -O3  -Wa,-q -w -std=legacy
     echo compile bust_5  
    gfortran $filebust -O3  -Wa,-q -w -std=legacy
     echo compile adjust
    gfortran $fileadj -O3  -Wa,-q -w -std=legacy
  elif [ "$2" = "i" ]
  then
    echo compiler is intel ifort
    echo compile est_noise
    /opt/intel/bin/compilervars.sh intel64
    /opt/intel/mkl/bin/mklvars.sh intel64
 
    ifort $file -w -i8 -I${MKLROOT}/include/intel64/ilp64 -I${MKLROOT}/include -O3 \
    ${MKLROOT}/lib/libmkl_blas95_ilp64.a ${MKLROOT}/lib/libmkl_intel_ilp64.a ${MKLROOT}/lib/libmkl_sequential.a ${MKLROOT}/lib/libmkl_core.a -lpthread -lm -ldl \
         -Wl,-stack_size,0x40000000 

#    ifort $file -i8 -I${MKLROOT}/include/intel64/ilp64 -I${MKLROOT}/include  ${MKLROOT}/lib/libmkl_blas95_ilp64.a ${MKLROOT}/lib/libmkl_lapack95_ilp64.a -L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl
#    ifort $file  -Wl,-stack_size,0x40000000 -O3 -i8 -I${MKLROOT}/include/intel64/ilp64 -I${MKLROOT}/include  ${MKLROOT}/lib/libmkl_blas95_ilp64.a ${MKLROOT}/lib/libmkl_lapack95_ilp64.a -L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl
#    ifort $file -Wl,-stack_size,0x40000000 -O3  -L/Library/Frameworks/Intel_MKL.framework/Versions/Current/lib/universal -lmkl_lapack -lmkl_core -lmkl_intel_lp64 -lmkl_intel_ilp64 -lmkl_intel_thread
    echo compile psd_calc7
    ifort -O3 $filepsd
    echo compile gen_noise
    ifort -O3 $filegen
     echo compile compare_wander
    ifort $filewand -O3
     echo compile bust_5
    ifort $filebust -O3
     echo compile adjust
    ifort $fileadj -O3
  else
    echo "Invalid compiler option"
    exit
  fi

else

  echo "invalid OS"
  exit
fi

ls -l bin/est_noise"$ver"
ls -l bin/psd_calc7
ls -l bin/gen_noise7.02
ls -l bin/compare_wander7.01
ls -l bin/bust_5
ls -l bin/adjust_1
exit
