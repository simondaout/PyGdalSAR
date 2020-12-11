# Compilation.. 
# WIP
ifort ACP_unw.f90  -L/usr/lib64 -llapack -lblas -lm -o ACP_unw
gfortran -w -o flatten_stack flatten_stack.f
gfortran   -cpp -O3 -fcray-pointer -fbackslash -ffixed-line-length-none -ffree-line-length-none -fno-underscoring -D_FILE_OFFSET_BITS=64 -o my_deroul_interf_filt my_deroul_interf_filt.f90

gfortran -cpp -O3 -fcray-pointer -fbackslash -ffixed-line-length-none -ffree-line-length-none -fno-underscoring -D_FILE_OFFSET_BITS=64 -I/home/cometsoft/Ubuntu/nsbas/contrib/fftw/api -I../../src/common/base -I/home/cometsoft/Ubuntu/nsbas/src/common/blas -I/home/cometsoft/Ubuntu/nsbas/src/common/math -I/home/cometsoft/Ubuntu/nsbas/include -I.  -c -o flatten_stack.o flatten_stack.f90
gfortran -cpp -O3 -fcray-pointer -fbackslash -ffixed-line-length-none -ffree-line-length-none -fno-underscoring -D_FILE_OFFSET_BITS=64 -I/home/cometsoft/Ubuntu/nsbas/contrib/fftw/api -I/home/cometsoft/Ubuntu/nsbas/src/common/base -I/home/cometsoft/Ubuntu/nsbas/src/common/blas -I/home/cometsoft/Ubuntu/nsbas/src/common/math -I/home/cometsoft/Ubuntu/nsbas/include -I.  -o ../bin/flatten_stack flatten_stack.o /home/cometsoft/Ubuntu/nsbas/src/common/roipac/libnsbas-roipac.a /home/cometsoft/Ubuntu/nsbas/src/common/fftw/libnsbas-fftw3f.a /home/cometsoft/Ubuntu/nsbas/src/common/math/libnsbas-math.a /home/cometsoft/Ubuntu/nsbas/src/common/base/libnsbas-base.a /home/cometsoft/Ubuntu/nsbas/src/common/lapack/libnsbas-lapack.a /home/cometsoft/Ubuntu/nsbas/src/common/blas/libnsbas-blas.a -lm
