%%                             BCH decoder                              %%
% This is only for extended BCH with 1 extra parity bit                 %
%  codeword is the received codeword                                    %
%  m is the field extension degree for Galois field                     %
%  k is the length of information bits                                  %
%  t is the error correct ability                                       % 
%  Error_location is the detected error location                        %
%  Word is the decoded codeword                                         %
%  Decoding_status=0,no errors, Decoding_status=1,decoding success,
%  Decoding_status = others, there is error, but without correcting      %
%% 
function [Error_location, Word, Decoding_status]=eBCHdecoder(codeword,m,k,t)

BCH_n=2^m-1;
Error_location=zeros(1,BCH_n+1);
Word=zeros(1,BCH_n+1);
Return=bchdec_mex(m,k,t,codeword(1:BCH_n)); 
Error_poly(1:BCH_n)=xor(Return(1,1:BCH_n),codeword(1:BCH_n));% get the error polynomial
Success=Return(1,length(Return));  % return the status of the decoding results, 0=no errors, 1= decode successfully, others= detecting errors, but decoding failure
if Success~=1  % indicate decoding failure
     Word=codeword;
     Error_location=zeros(1,BCH_n+1);
     Decoding_status=Success; % Sucess=0, -1, -2
else
     d=sum(Error_poly(1:BCH_n)); %indicate the distance between the received codeword and the transmitted codeword
     de=mod(d+sum(codeword(1:BCH_n+1)),2);
     if d+de<=t
          Word(1:BCH_n)=Return(1,1:BCH_n);
          Word(BCH_n+1)=xor(codeword(BCH_n+1),de); 
          Error_location(1,1:BCH_n)=Error_poly(1:BCH_n);
          Error_location(1,BCH_n+1)=xor(Word(BCH_n+1),codeword(BCH_n+1));
          Decoding_status=1;
      else
          Word=codeword;
          Error_location=zeros(1,BCH_n+1);
          Decoding_status=-1;
     end   
end
   