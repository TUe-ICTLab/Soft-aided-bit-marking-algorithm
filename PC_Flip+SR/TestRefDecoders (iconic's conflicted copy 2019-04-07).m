%% Script to check the performance of reference decoders for PC, such as Chase-Pyndiah and Miscorrection-Free
It=10;
BCHnu=8;
BCHt=2;
BCHe=1;
BCHs=0;
BCHn=2^BCHnu-1+BCHe-BCHs;
BCHk=2^BCHnu-1-BCHnu*BCHt;

modord=1;
nblock=1e5;
EbNo=3.7:.1:4.3;
r=(BCHk/BCHn)^2;

Save=1;

BER=zeros(1,length(EbNo));

for ss=1:length(EbNo)
disp(['Estimating BER for EbNo=' num2str(EbNo(ss)) 'dB....']);
[~,BER(ss)]=Produ_IMP_cod_BCH_AWGN_chase_pyndiah_matlab_v1(It,EbNo(ss),BCHnu,BCHt,BCHs,modord,nblock);
disp(['EbNo=' num2str(EbNo(ss))   'dB_BER=' num2str(BER(ss))]);
end

semilogy(EbNo,BER,'sr-');
grid on;hold on;
axis([2.5,4,1e-9,1e-1]);

if Save
%spath='/home/gliga/Dropbox/Work/TUe/MATLAB/MATLAB Sims Data/EnhancedHDdecoding/MarkedSR_PCdecoding/NewResults/';    
spath='E:\Users\Gabriele\Dropbox\Work\TUe\MATLAB\MATLAB Sims Data\EnhancedHDdecoding\MarkedSR_PCdecoding\NewResults\';    
filename=['TPC_Chase-Pyndiah_PAM2_PC_BCH_n=' num2str(BCHn) '_t=' num2str(BCHt) '_e= ' num2str(BCHe) '_NoIt=10.mat' ];
fpath=[spath filename];
save(fpath, 'EbNo','BER','r');
end


