%Miscorrection-free performance of an iBDD decoder for product codes

iter=10;     % Number of iterations performed for iBDD
SNR=5:.1:8   ;     % 2*Eb/No
modord=1;    % Constellation cardinality 2^modord
v=7;         % Order of the extension field for BCH component code
t=2;         % Correction power  
s=0;         % shortening factor 
BER=zeros(1,length(SNR));
Save=1;

parfor ii=1:length(SNR)
 disp(['Calculating Miscorrection Free BER for SNR=' num2str(SNR(ii)) 'dB']);
[~,BER(ii)]=Product_IMP_code_BCH_AWGN_normnoise_eBCH_misfree(iter,SNR(ii)-3,v,t,s,modord);
disp(['SNR= ' num2str(SNR(ii)) 'dB BER=' num2str(BER(ii))]);
end

figure;
semilogy(SNR,BER,'sr--');
axis([5,7.5,1e-3,1e-1]);
grid on;

if Save
  if isunix
s_path='/home/gabriele/Dropbox/Work/TUe/MATLAB/MATLAB Sims Data/EnhancedHDdecoding/MarkedSR_PCdecoding/';
    elseif ispc
s_path='E:\Dropbox\Work\TUe\MATLAB\MATLAB Sims Data\EnhancedHDdecoding\MarkedSR_PCdecoding\';
    else  
     error('Unexpected operative system');   
   end
end 
    
filename=['MiscFree_PC_BCH_v=' num2str(v) 't=' ]
save()