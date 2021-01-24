%% Plots results on different bit marking algorithms for product codes and PAM2 in AWGN
clear all;
close all;
tikz=1;

BCHn=127;
M=[2,4,8];

if ispc
l_path='E:\Users\Gabriele\Dropbox\Work\TUe\MATLAB\MATLAB Sims Data\EnhancedHDdecoding\MarkedSR_PCdecoding\OECC_Results\';
elseif isunix 
%l_path='/home/gliga/Dropbox/Work/TUe/MATLAB/MATLAB Sims Data/EnhancedHDdecoding/MarkedSR_PCdecoding/OECC_Results/';    
l_path='/home/gabriele/Dropbox/Work/TUe/MATLAB/MATLAB Sims Data/EnhancedHDdecoding/MarkedSR_PCdecoding/OECC_Results/';    

end

%Color=rand(,3);

for mm=1:length(M)

if isequal(M(mm),2)
filename1=['iBDD_PAM' num2str(M(mm)) '_PC_BCH_n=' num2str(BCHn) '_t=2_e=1_NoIt=10.mat'];
filename2=['iBDD-SABM_PAM'  num2str(M(mm)) '_PC_BCH_n=' num2str(BCHn) '_t=2_e=1_NoIt=10_R=5_SABM_NoIt=5.mat'];
else
filename1=['iBDD_PAM' num2str(M(mm)) '_PC_BCH_n=' num2str(BCHn) '_t=2_e=1_NoIt=10_InterLeaved.mat'];
filename2=['iBDD-SABM_PAM'  num2str(M(mm)) '_PC_BCH_n=' num2str(BCHn) '_t=2_e=1_NoIt=10_R=5_SABM_NoIt=5_InterLeaved.mat'];
end
    
filepath1=[l_path filename1];
filepath2=[l_path filename2];

load(filepath1);

EbNo=snrdB-3-10*(log10(r)+log10(log2(M(mm))));
%p1=semilogy(EbNo,BER,'-rs','LineWidth',1.5,'MarkerFaceColor','w');
p1=semilogy(snrdB,BER,'-rs','LineWidth',1.5,'MarkerFaceColor','w');
hold on;
grid on;

load(filepath2);
EbNo=snrdB-3-10*(log10(r)+log10(log2(M(mm))));
%p2=semilogy(EbNo,BER,'-bs','LineWidth',1.5,'MarkerFaceColor','w');
p2=semilogy(snrdB,BER,'-bs','LineWidth',1.5,'MarkerFaceColor','w');

if ~isequal(M(mm),2)
filename3=['iBDD-SABM_PAM'  num2str(M(mm)) '_PC_BCH_n=' num2str(BCHn) '_t=2_e=1_NoIt=10_R=5_SABM_NoIt=5_InterLeaved_MaxLog.mat'];    
filepath3=[l_path filename3];
load(filepath3);
EbNo=snrdB-3-10*(log10(r)+log10(log2(M(mm))));
%p3=semilogy(EbNo,BER,'--k*','LineWidth',1.5,'MarkerFaceColor','w');
p3=semilogy(snrdB,BER,'--k*','LineWidth',1.5,'MarkerFaceColor','w');
end
end


axis([3,21,1e-8,1e-1]);
xlabel('SNR [dB]');
ylabel('BER');

if tikz
    %addpath(genpath('/home/gliga/Dropbox/Work/Matlab code/matlab2tikz-matlab2tikz-816f875'));
    addpath(genpath('/home/gabriele/Dropbox/Work/Matlab code/matlab2tikz-matlab2tikz-816f875'));

    tikzfilename='OECCresults.tikz';
    %tpath='/home/gliga/Dropbox/Work/TUe/MATLAB/MATLAB Sims Data/EnhancedHDdecoding/MarkedSR_PCdecoding/NewResults/tikz/';
    tpath='/home/gabriele/Dropbox/Work/TUe/MATLAB/MATLAB Sims Data/EnhancedHDdecoding/MarkedSR_PCdecoding/NewResults/tikz/';

    tikzfilepath=[tpath,tikzfilename];
    matlab2tikz(tikzfilepath);
end


%legend(p,'1 It','2 It','3 It','4 It','5 It');
%legend(p,'4 It','5 It');