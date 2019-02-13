# Compilation.. 
# WIP
ifort ACP_unw.f90  -L/usr/lib64 -llapack -lblas -lm -o ACP_unw
gfortran -w -o flatten_model flatten_model.f
