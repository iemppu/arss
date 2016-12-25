/* -------------------------------------------------------------------------- */
/* partXY_mex mexFunction */
/* -------------------------------------------------------------------------- */

/* Compile with mex -lmwlapack -lmwblas -largeArrayDims partXY_blas.c */

#include "mex.h"
#include "blas.h"


// double ddot_(int *n, double *sx, int *incx, double *sy, int *incy)
// {
//   long int i, m, nn, iincx, iincy;
//   double stemp;
//   long int ix, iy;
// 
//   /* forms the dot product of two vectors.   
//      uses unrolled loops for increments equal to one.   
//      jack dongarra, linpack, 3/11/78.   
//      modified 12/3/93, array(1) declarations changed to array(*) */
// 
//   /* Dereference inputs */
//   nn = *n;
//   iincx = *incx;
//   iincy = *incy;
// 
//   stemp = 0.0;
//   if (nn > 0)
//   {
//     if (iincx == 1 && iincy == 1) /* code for both increments equal to 1 */
//     {
//       m = nn-4;
//       for (i = 0; i < m; i += 5)
//         stemp += sx[i] * sy[i] + sx[i+1] * sy[i+1] + sx[i+2] * sy[i+2] +
//                  sx[i+3] * sy[i+3] + sx[i+4] * sy[i+4];
// 
//       for ( ; i < nn; i++)        /* clean-up loop */
//         stemp += sx[i] * sy[i];
//     }
//     else /* code for unequal increments or equal increments not equal to 1 */
//     {
//       ix = 0;
//       iy = 0;
//       if (iincx < 0)
//         ix = (1 - nn) * iincx;
//       if (iincy < 0)
//         iy = (1 - nn) * iincy;
//       for (i = 0; i < nn; i++)
//       {
//         stemp += sx[ix] * sy[iy];
//         ix += iincx;
//         iy += iincy;
//       }
//     }
//   }
// 
//   return stemp;
// } /* ddot_ */



/* compute a part of X*Y */
void mexFunction
(
    int nargout,
    mxArray *pargout [ ],
    int nargin,
    const mxArray *pargin  [ ]
)
{
    double *Xt, *Y, *Z, *I, *J, *v, LL;
    ptrdiff_t m, n, r, L, p, ir, jr, k;
    ptrdiff_t inc = 1;

    if (nargin != 5 || nargout > 1)
        mexErrMsgTxt ("Usage: v = partXY (Xt, Y, I, J, L)") ;

    /* ---------------------------------------------------------------- */
    /* inputs */
    /* ---------------------------------------------------------------- */
    
    Xt = mxGetPr( pargin [0] );
    Y  = mxGetPr( pargin [1] );
    I  = mxGetPr( pargin [2] );
    J  = mxGetPr( pargin [3] );
    LL = mxGetScalar( pargin [4] ); L = (ptrdiff_t) LL;
    m  = mxGetN( pargin [0] );
    n  = mxGetN( pargin [1] );
    r  = mxGetM( pargin [0] ); 
    if ( r != mxGetM( pargin [1] ))
        mexErrMsgTxt ("rows of Xt must be equal to rows of Y") ;
    if ( r > m || r > n )
        mexErrMsgTxt ("rank must be r <= min(m,n)") ;
    
    /* ---------------------------------------------------------------- */
    /* output */
    /* ---------------------------------------------------------------- */

    pargout [0] = mxCreateDoubleMatrix(1, L, mxREAL);
    v = mxGetPr( pargout [0] );
    
    /* C array indices start from 0 */
    for (p = 0; p < L; p++) {
        ir = ( I[p] - 1 ) * r;
        jr = ( J[p] - 1 ) * r;
         v[p] = 0.0;
       /* for (k = 0; k < r; k++)
            v[p] += Xt[ ir + k ] * Y[ jr + k ]; */
        v[p] = ddot(&r, Xt+ir, &inc, Y+jr, &inc);
    }
    
    return;
}

