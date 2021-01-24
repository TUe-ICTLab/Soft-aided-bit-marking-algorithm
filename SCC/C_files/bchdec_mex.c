/*
 * bchdec_mex.cpp
 * The MEX interface function bchdec_mex provide interface to the
 * the BCH decoder function implemented in C.
 *
 * Syntax
 *  message = bchdec_mex(m, k, t, rx)
 * Arguments
 *  message - Decoded message, binary row vector
 *  m - Order of the Galois Field, scalar, 2~20
 *  k - Length of the information bits, scalar
 *  t - Error correction capability, scalar
 *  rx - Received corrupted codeword, binary row vector
 * Note
 *  Not all m, n, t combinations are valid.
 *
 * Author: Alan Bao Jian ZHOU
 * Created: 20141022
 * Last-modified: 20141022
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
 *      bchdec_mex.cpp
 *      bchdec.h
 *      bchdec.c
 * are required for compiling the MEX function bchdec_mex.
 *
 * 2. Compiling command
 *      mex -O bchdec_mex.cpp bchdec.c
 * You can control optimization options, e.g., in Windows,
 *      mex -v -largeArrayDims -O OPTIMFLAGS="/Ox" bchdec_mex.cpp bchdec.c
 *
 * 3. For debug purpose (e.g., MSVC++|Tools|Attache to Process), use
 *      mex -g bchdec_mex.cpp bchdec.c
 * to build the file with debugging symbols included.
 */


#include "mex.h"
#include "matrix.h" // MEX functions

// The header file that contains the C function this MEX function is interfacing to
#include "bchdec.h"

#include <stdlib.h> // malloc


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // Verify MEX-File Input and Output Parameters.
    // 4 inputs
    if(nrhs != 4)
    {
        mexErrMsgIdAndTxt("bchdec_mex:nrhs",
                "4 inputs required.");
    }
    // 1 output
    if(nlhs != 1)
    {
        mexErrMsgIdAndTxt("bchdec_mex:nlhs",
                "1 output required.");
    }
    // m is a scalar.
    if( !mxIsDouble(prhs[0]) ||
         mxIsComplex(prhs[0]) ||
         mxGetNumberOfElements(prhs[0]) != 1 )
    {
        mexErrMsgIdAndTxt("bchdec_mex:mNotScalar",
                "Input m must be a scalar.");
    }
    // k is a scalar.
    if( !mxIsDouble(prhs[1]) ||
         mxIsComplex(prhs[1]) ||
         mxGetNumberOfElements(prhs[1]) != 1 )
    {
        mexErrMsgIdAndTxt("bchdec_mex:nNotScalar",
                "Input k must be a scalar.");
    }
    // t is a scalar.
    if( !mxIsDouble(prhs[2]) ||
         mxIsComplex(prhs[2]) ||
         mxGetNumberOfElements(prhs[2]) != 1 )
    {
        mexErrMsgIdAndTxt("bchdec_mex:tNotScalar",
                "Input t must be a scalar.");
    }
    // rx is type double.
    if( !mxIsDouble(prhs[3]) ||
         mxIsComplex(prhs[3]))
    {
        mexErrMsgIdAndTxt("bchdec_mex:messageNotDouble",
                "Input rx must be type double.");
    }
    // rx is a row vector.
    if(mxGetM(prhs[3])!=1)
    {
        mexErrMsgIdAndTxt("bchdec_mex:messageNotRowVector",
                "Input rx must be a row vector.");
    }
    
    // Read Input Data.
    int m = (int)mxGetScalar(prhs[0]); // Order of the Galois Field, input
    int k = (int)mxGetScalar(prhs[1]); // Length of information bits, input
    int t = (int)mxGetScalar(prhs[2]); // Error correction capacity, input
    double *rx = mxGetPr(prhs[3]); // Received corrupted codeword, input
    int n = (int)mxGetN(prhs[3]); // Length of the codeword
    
    // Prepare decoding.
    // Allocate memory, type unsigned char for compatibility with bchdec.c
    // Remember to free.
    unsigned char *ucRx = (unsigned char*)malloc(n * sizeof(unsigned char));
    // Copy rx into ucRx.
    // Caution: Adjust the codeword s.t. it's compatible with the MATLAB built-in function bchdec.
    // In MATLAB, message: ->, redundancy: ===>, codeword: ->===>
    // In bchdec.c, the LSB and MSB are reversed w.r.t. that in MATLAB.
	// To obtain the same decoding result, the input should be codeword: <-<===, and the message is <-.
    int idxRx = 0;
    int idxUcRx = k - 1;
    for (idxRx = 0; idxRx < k; idxRx++) // Flip the message part.
    {
        ucRx[idxUcRx] = (unsigned char)rx[idxRx];
        idxUcRx--;
    }
    idxUcRx = n - 1;
    for (idxRx = k; idxRx < n; idxRx++) // Flip the redundancy part.
    {
        ucRx[idxUcRx] = (unsigned char)rx[idxRx];
        idxUcRx--;
    }
    
    // Decode
    int flag = bchdec(m, n, t, ucRx);
    
    // Return differently according to different decoding results
    // 0: no/undetectable error; 1: error (erroneously) corrected; -1/-2: error detected.
    if(flag == 0 || flag == 1) // 0: no/undetectable error; 1: error (erroneously) corrected.
    {
        // Prepare Output Data.
        plhs[0] = mxCreateDoubleMatrix(1,n+1,mxREAL); //ComplexFlag= mxREAL, only real data in mxArray; ComplexFlag= mxCOMPLEX, mxArray contains imaginary data
        double *message = mxGetPr(plhs[0]); // Decoded message, output
        
        // Post-processing after decoding
        // Copy unsigned char message into double message.
        // Caution: Adjust the codeword s.t. it's compatible with the MATLAB built-in function bchdec.
        // In MATLAB, message: ->, redundancy: ===>, codeword: ->===>
        // In bchdec.c, the LSB and MSB are reversed w.r.t. that in MATLAB.
        // To obtain the same decoding result, the input should be codeword: <-<===, and the message is <-.
        int idxMessage = k - 1;
        for (idxUcRx = 0; idxUcRx < k; idxUcRx++) // Flip the message part.
        {
            message[idxMessage] = ucRx[idxUcRx];
            idxMessage--;
        }

         idxMessage = n - 1;
         for (idxUcRx = k; idxUcRx < n; idxUcRx++) // Flip the redundancy part.
         {
            message[idxMessage] = ucRx[idxUcRx];
            idxMessage--;
         }
        message[n]=flag; // restore the error status.0= no errors, 1=exist error.
    }
    else if(flag == -1 || flag == -2) // -1/-2: error detected.
    {
        // Prepare Output Data.
      // This location has been modified by Yi Lei
        //plhs[0] = mxCreateDoubleScalar(flag); // Return flag
         plhs[0] = mxCreateDoubleMatrix(1,n+1,mxREAL); 
         double *message = mxGetPr(plhs[0]); // Decoded message, output
         for (idxUcRx = 0; idxUcRx < n; idxUcRx++) // Flip the message part.
         {
            message[idxUcRx] = (unsigned char)rx[idxUcRx];
         }
          message[n]=flag;// restore the error status.0= no errors, 1=exist error.
         
    }
    else // Decoding error happens.
    {
        // Prepare Output Data.
        //This location has been modified by Yi Lei
         plhs[0] = mxCreateDoubleMatrix(1,n+1,mxREAL); 
         double *message = mxGetPr(plhs[0]); // Decoded message, output
         for (idxUcRx = 0; idxUcRx < n; idxUcRx++) // Flip the message part.
         {
            message[idxUcRx] = (unsigned char)rx[idxUcRx];
         }
          message[n]=flag;// restore the error status.0= no errors, 1=exist error.
       // plhs[0] = mxCreateDoubleScalar(-3);
       // mexErrMsgIdAndTxt("bchdec_mex:error",
              //  "Decoding error happens.");
    } // if flag
    
    // Release memory.
    free(ucRx);
} // fxn mexFunction
