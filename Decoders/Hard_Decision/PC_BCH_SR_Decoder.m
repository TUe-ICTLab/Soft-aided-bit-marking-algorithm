
function [hat_coded_bits,r]=PC_BCH_SR_Decoder(coded_bits,BCHnu,BCHt,BCHe,Ni,L,w)
%% PC decoder using iBDD and scaled reliabilities (see [1])
% Decoding function for a product code.
% It uses the function eBCHdecoder (from Yi Lei) which assumes the BCH code
% was extended by one bit. BCH codes with no extension or with more than
% one extension bits (C. Haeger uses 2 for example) has not been
% implemented yet.

% Reference 
% [1] Sheikh, A., Amat, A. G. i, & Liva, G. (2019). Binary Message Passing Decoding of Product Codes Based on Generalized Minimum Distance Decoding. ArXiv:1901.02914. Retrieved from http://arxiv.org/abs/1901.02914

% Inputs
% coded_bits         Input hard bits to the decoder
% BCHnu,BCHt,BCHe    Parameters of the eBCH component code
% Ni                 Number of iterations for iBDD decoding
% L                  Block of BCHn x BCHn LLRs passed to the decoder
% w                  Vector of weigths for scaled reliability update

% Outputs
% r                  BCHn x BCHn x Ni matrix containing scaled reliabilities for each code array nxn at each iteration

   
    
BCHn=(2^BCHnu)-1;                % Coded bits, without the extended bits
BCHk=(2^BCHnu)-1-BCHnu*BCHt;     % Information bits
r=zeros(BCHn+BCHe,BCHn+BCHe,Ni); % matrix of scaled reliabilities (for each iteration)

%ItThresh=2;                      % Iterations threshold for zero-syndrome decoding rejection

hat_coded_bits=coded_bits;
%CheckConflict=zeros(1,BCHn+BCHe);
%ZeroSyndromeIdx=zeros(1,BCHn+BCHe); % Initialise as none of the codewords having 0 syndrome
status=zeros(1,BCHn+BCHe);

%% Iterative decoding loop
for ll=1:Ni 
    % Decode Rows
    for rr=1:BCHn+BCHe
        [tmp_err_poly, tmp_decoded, status(rr)]=eBCHdecoder(hat_coded_bits(rr,:),BCHnu,BCHk,BCHt);
        
        if ~(status(rr)==0 || status(rr)==1)    % Decoding failure
         u=zeros(1,BCHn+BCHe);
        else
          u=2*tmp_decoded-1;
        end
         
        r(rr,:,ll)=w(ll)*u+L(rr,:);                                                 % scaled reliabilities 
        hat_coded_bits(rr,:)=(r(rr,:,ll)>0);
 
       
    end
       
    % Decode Columns
    for cc=1:BCHn+BCHe 
        [tmp_err_poly, tmp_decoded, status(cc)]=eBCHdecoder(hat_coded_bits(:,cc).',BCHnu,BCHk,BCHt);
        
         if ~(status(cc)==0 || status(cc)==1)    % Decoding failure
         u=zeros(BCHn+BCHe,1);
        else
          u=2*tmp_decoded.'-1;
        end
         
        r(:,cc,ll)=w(ll)*u+L(:,cc);                                                 % scaled reliabilities 
        hat_coded_bits(:,cc)=(r(:,cc,ll)>0);
                  
        
    end
         %ZeroSyndromeIdx=(status==0);

  %  hat_coded_bits(Mark)=coded_bits(Mark);    % keep reliable bits as they were before the decoding at the end of each iteration
end