# Compilation.. 
# WIP
ifort ACP_unw.f90  -L/usr/lib64 -llapack -lblas -lm -o ACP_unw
gfortran -w -o flatten_model flatten_model.f
gfortran   -cpp -O3 -fcray-pointer -fbackslash -ffixed-line-length-none -ffree-line-length-none -fno-underscoring -D_FILE_OFFSET_BITS=64 -o my_deroul_interf_filt my_deroul_interf_filt.f90
