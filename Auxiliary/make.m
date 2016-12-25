% This make.m is used under Windows
cd D:\tmingkui\matrix_completion\RiemannianMatrixCompletion_17Aug2012\Auxiliary
mex -lmwlapack -lmwblas -largeArrayDims updateSval_blas.c 
mex -lmwlapack -lmwblas -largeArrayDims partXY_blas.c 
% mex -O -largeArrayDims -c partXY_blas.c
mex -O -largeArrayDims -c partXY.c
% mex -O -largeArrayDims -c partXY_blas.c
% mex -O -largeArrayDims -c updateSval.c
% mex -O -largeArrayDims -c updateSval_blas.c
mex -O -largeArrayDims  partXY.obj
% mex -O -largeArrayDims  partXY_blas.obj
% mex -O -largeArrayDims  partXY.c
% mex -O -largeArrayDims  partXY_blas.obj
% mex -O -largeArrayDims  updateSval.obj
% mex -O -largeArrayDims  updateSval_blas.obj

%mex -O -largeArrayDims mex_sparse_coding.cpp matrix_op.obj conjugate.obj cover_tree.obj libomp.obj myblas.obj omputils.obj mexvertification.obj
%mex -O -largeArrayDims train.c -I..\ tron.obj linear.obj linear_model_matlab.obj ..\blas\*.obj
%mex -O -largeArrayDims predict.c -I..\ tron.obj linear.obj linear_model_matlab.obj ..\blas\*.obj
% mex -O -largeArrayDims libsvmread.c
% mex -O -largeArrayDims libsvmwrite.c  -lmwlapack -lmwblas
