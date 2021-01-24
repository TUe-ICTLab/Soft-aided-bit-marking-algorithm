/*
 * bchenc_mex.cpp
 * The MEX interface function bchenc_mex provide interface to the
 * the BCH encoder function implemented in C.
 *
 * Syntax
 *  codeword = bchenc_mex(m, n, t, message)
 * Arguments
 *  codeword - BCH codeword, binary row vector
 *  m - Order of the Galois Field, scalar, 2~20
 *  n - Length of the codeword, codeword, scalar
 *  t - Error correction capability, scalar
 *  message - Information bits, binary row vector
 * Note
 *  Not all m, n, t combinations are valid.
 *
 * Author: Alan Bao Jian ZHOU
 * Created: 20141021
 * Last-modified: 20141021
 */

/*
 * Note
 * 0. Make sure that you have installed the C/C++ compiler
 * (comes with gcc/g++, Windows SDK, or MSVC++) that is supported
 * by your version of MATLAB. Also, remember to use
 *      mex -setup
 * to assign a C/C++ compiler to MATLAB MEX if this is the first
 * time you compile a mex function.
 *
 * 1. The source files
 *      bchenc_mex.cpp
 *      bchenc.h
 *      bchenc.c
 * are required for compiling the MEX function bchenc_mex.
 *
 * 2. Compiling command
 *      mex -O bchenc_mex.cpp bchenc.c
 * You can control optimization options, e.g., in Windows,
 *      mex -v -largeArrayDims -O OPTIMFLAGS="/Ox" bchenc_mex.cpp bchenc.c
 *
 * 3. For debug purpose (e.g., MSVC++|Tools|Attache to Process), use
 *      mex -g bchenc_mex.cpp bchenc.c
 * to build the file with debugging symbols included.
 */


#include "mex.h"
#include "matrix.h" // MEX functions

// The header file that contains the C function this MEX function is interfacing to
#include "bchenc.h"

#include <stdlib.h> // malloc


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // Verify MEX-File Input and Output Parameters.
    // 4 inputs
    if(nrhs != 4)
    {
        mexErrMsgIdAndTxt("bchenc_mex:nrhs", "4 inputs required.");
    }
    // 1 output
    if(nlhs != 1)
    {
        mexErrMsgIdAndTxt("bchenc_mex:nlhs", "1 output required.");
    }
    // m is a scalar.
    if( !mxIsDouble(prhs[0]) ||
         mxIsComplex(prhs[0]) ||
         mxGetNumberOfElements(prhs[0]) != 1 )
    {
        mexErrMsgIdAndTxt("bchenc_mex:mNotScalar",
                "Input m must be a scalar.");
    }
    // n is a scalar.
    if( !mxIsDouble(prhs[1]) ||
         mxIsComplex(prhs[1]) ||
         mxGetNumberOfElements(prhs[1]) != 1 )
    {
        mexErrMsgIdAndTxt("bchenc_mex:nNotScalar",
                "Input n must be a scalar.");
    }
    // t is a scalar.
    if( !mxIsDouble(prhs[2]) ||
         mxIsComplex(prhs[2]) ||
         mxGetNumberOfElements(prhs[2]) != 1 )
    {
        mexErrMsgIdAndTxt("bchenc_mex:tNotScalar",
                "Input t must be a scalar.");
    }
    // message is type double.
    if( !mxIsDouble(prhs[3]) ||
         mxIsComplex(prhs[3]))
    {
        mexErrMsgIdAndTxt("bchenc_mex:messageNotDouble",
                "Input message must be type double.");
    }
    // message is a row vector.
    if(mxGetM(prhs[3])!=1)
    {
        mexErrMsgIdAndTxt("bchenc_mex:messageNotRowVector",
                "Input message must be a row vector.");
    }
    
    // Read Input Data.
    int m = (int)mxGetScalar(prhs[0]); // Order of the Galois Field, input
    int n = (int)mxGetScalar(prhs[1]); // Length of the codeword, input
    int t = (int)mxGetScalar(prhs[2]); // Error correction capacity, input
    double *message = mxGetPr(prhs[3]); // Information bits, input
    int k = (int)mxGetN(prhs[3]); // Length of information bits

    // Prepare Output Data.
    plhs[0] = mxCreateDoubleMatrix(1,n,mxREAL);
    double *codeword = mxGetPr(plhs[0]); // Codeword bits, output
    
    // Prepare encoding.
    // Allocate memory, type unsigned char for compatibility with bchenc.c
    // Remember to free.
    unsigned char *ucCodeword = (unsigned char*)malloc(n * sizeof(unsigned char));
    // Copy message into codeword.
    // Caution: Adjust the codeword s.t. it's compatible with the MATLAB built-in function bchenc.
    // In MATLAB, message: ->, redundancy: ===>, codeword: ->===>
    // In bchenc.c, the LSB and MSB are reversed w.r.t. that in MATLAB.
    // To obtain the same encoding result, the input should be message: <-, and the codeword is <-<===.    
    int idxMessage = 0;
    int idxUcCodeword = k - 1;
    for (idxMessage = 0; idxMessage < k; idxMessage++) // Flip the message.
    {
        ucCodeword[idxUcCodeword] = (unsigned char)message[idxMessage];
        idxUcCodeword--;
    }
    
    // Encode
    bchenc(m, n, t, ucCodeword);
    
    // Post-processing after encoding
    // Copy unsigned char codeword into double codeword.
    // Caution: Adjust the codeword s.t. it's compatible with the MATLAB built-in function bchenc.
    // In MATLAB, message: ->, redundancy: ===>, codeword: ->===>
    // In bchenc.c, the LSB and MSB are reversed w.r.t. that in MATLAB.
    // To obtain the same encoding result, the input should be message: <-, and the codeword is <-<===.
    int idxCodeword = k - 1;
    for (idxUcCodeword = 0; idxUcCodeword < k; idxUcCodeword++) // Flip the message part.
    {
        codeword[idxCodeword] = ucCodeword[idxUcCodeword];
        idxCodeword--;
    }
    idxCodeword = n - 1;
    for (idxUcCodeword = k; idxUcCodeword < n; idxUcCodeword++) // Flip the redundancy part.
    {
        codeword[idxCodeword] = ucCodeword[idxUcCodeword];
        idxCodeword--;
    }
    
    // Release memory.
    free(ucCodeword);
} // fxn mexFunction
