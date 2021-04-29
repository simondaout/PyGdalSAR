# Compilation.. 
# WIP
ifort ACP_unw.f90  -L/usr/lib64 -llapack -lblas -lm -o ACP_unw
gfortran -w -o ../bin/flatten_stack flatten_stack.f
gfortran -w -o ../bin/unflatten_stack unflatten_stack.f90

