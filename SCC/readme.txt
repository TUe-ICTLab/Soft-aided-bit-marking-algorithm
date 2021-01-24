readme.txt
readme file for bch codec.
Author: Alan Baojian ZHOU
Created: 9:05 PM 10/22/2014
Last-Modified: 9:12 PM 10/22/2014


Usage

1. Follow instructions in bchenc_mex.cpp and bchdec_mex.cpp to
1) set up MATLAB MEX, and
2) compile the MEX functions bchenc_mex and bchenc_mex.

2) Run the script TestBchCodec.m to get a quick idea of how to use the generated MEX functions. You are expected to get an output as the SampleOutput.txt.


Note

If you do not have the communication toolbox in your MATLAB, the script TestBchCodec.m cannot run since the functions bchenc and bchdec are missing. However, you can still use the bchenc_mex and bchdec_mex to do BCH encoding and decoding, which is exactly what this project aims at.


Version Log

v0.1
1. First version.