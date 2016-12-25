% This make.m is used under Windows

mex -O -largeArrayDims -c matrix_op.cpp
mex -O -largeArrayDims -c conjugate.cpp
mex -O -largeArrayDims -c cover_tree.cpp
mex -O -largeArrayDims -c libomp.cpp
mex -O -largeArrayDims -c myblas.cpp
mex -O -largeArrayDims -c omputils.cpp
mex -O -largeArrayDims -c mexvertification.cpp

mex -O -largeArrayDims mex_sparse_coding.cpp matrix_op.obj conjugate.obj cover_tree.obj libomp.obj myblas.obj omputils.obj mexvertification.obj
%mex -O -largeArrayDims train.c -I..\ tron.obj linear.obj linear_model_matlab.obj ..\blas\*.obj
%mex -O -largeArrayDims predict.c -I..\ tron.obj linear.obj linear_model_matlab.obj ..\blas\*.obj
% mex -O -largeArrayDims libsvmread.c
% mex -O -largeArrayDims libsvmwrite.c
