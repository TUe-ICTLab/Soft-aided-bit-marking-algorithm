%% Project folder containing simulations of enhanced "soft-aided" algebraic decoding
%% of (degree-2) product codes using channel reliabilities (LLRs). 
%% It merges iBDD scaled reliability approach from Alireza Sheikh (see [1]) and Yi Lei's bit marking and flipping (see [2])

% [1] Sheikh, A., Amat, A. G. i, & Liva, G. (2019). Binary Message Passing Decoding of Product Codes Based on Generalized Minimum Distance Decoding. ArXiv:1901.02914. Retrieved from http://arxiv.org/abs/1901.02914
% [2] Lei, Y., Chen, B., Liga, G., Deng, X., Cao, Z., Li, J., … Alvarado, A. (2019). Improved Decoding of Staircase Codes: The Soft-aided Bit-marking (SABM) Algorithm, 1–10. Retrieved from http://arxiv.org/abs/1902.01178


%% Algorithm description
% The algorithm is based on iBDD of product codes (BCH as a compontent code) 
% with the help the channel reliabilities to improve the algebraic decoding perfomance. 
% BDD suffers from the following limitations: it only succeds when received
% codeword is within t errors from the transmitted one. Fails or incurs
% miscorrections otherwise. Failures on the component code can be improved (other than by the iterative decoding) also 
% using updated reliabilites and generating updated (at each iteration)hard bits based on them.  
% BDD output can indicate whether the decoder thinks a certain bits is correct (either 0 or 1) or if it fails just declares
% all the bits being equal to zero. The updated reliability is then found
% as ri=ui+wi*Li where ui:{+,1,0,-1} and wi is an optimised weigth and Li
% is the channel reliability. On top of this miscorrections are detected
% and prevented using syndrome information and bit marking based on the
% updated reliabilities. 

