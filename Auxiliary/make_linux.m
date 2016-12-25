% This make.m is used under Linux
mex CFLAGS="\$CFLAGS -std=c99" -O -lblas -largeArrayDims updateSval.c
% mex CFLAGS="\$CFLAGS -std=c99" -O -lblas -largeArrayDims updateSval_blas.c

mex CFLAGS="\$CFLAGS -std=c99" -O -largeArrayDims partXY.c



% mex -O -largeArrayDims updateSval.c
mex CFLAGS="\$CFLAGS -std=c99" -O -largeArrayDims -lblas partXY_blas.c
% mex -O -largeArrayDims -lblas partXY_blas.c
% mex -O -largeArrayDims -lblas updateSval.c