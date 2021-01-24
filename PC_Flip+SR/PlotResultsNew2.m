%% Plots results on iBDD/SABM algorithms for product codes and PAM2 in AWGN with different BCH component codes (after TCOMM revision)

clear all;
close all;

BCHv=9;
BCHn=2.^BCHv-1;
BCHt=2;

tikz=0;

if ispc
l_path1='E:\Users\Gabriele\Dropbox\Work\TUe\MATLAB\MATLAB Sims Data\EnhancedHDdecoding\MarkedSR_PCdecoding\NewResults\iBDD\';
l_path2='E:\Users\Gabriele\Dropbox\Work\TUe\MATLAB\MATLAB Sims Data\EnhancedHDdecoding\MarkedSR_PCdecoding\NewResults\SABM\';


elseif isunix 
%l_path1='/home/gabriele/Dropbox/Work/TUe/MATLAB/MATLAB Sims Data/EnhancedHDdecoding/MarkedSR_PCdecoding/NewResults/SABM/';
%l_path1='/home/gliga/Dropbox/Work/TUe/MATLAB/MATLAB Sims Data/EnhancedHDdecoding/MarkedSR_PCdecoding/NewResults/iBDD/';
%l_path2='/home/gliga/Dropbox/Work/TUe/MATLAB/MATLAB Sims Data/EnhancedHDdecoding/MarkedSR_PCdecoding/NewResults/SABM/';
l_path1='/home/gabriele/Dropbox/Work/TUe/MATLAB/MATLAB Sims Data/EnhancedHDdecoding/MarkedSR_PCdecoding/NewResults/iBDD/';
l_path2='/home/gabriele/Dropbox/Work/TUe/MATLAB/MATLAB Sims Data/EnhancedHDdecoding/MarkedSR_PCdecoding/NewResults/SABM/';
      

end

Color=[0 0 1; 1 0 1 ; 1 0.6 0.7 ; 0 1 1];

figure;

for ll=1:length(BCHt)

filename1=['iBDD_PAM2_PC_BCH_n=' num2str(BCHn) '_t=' num2str(BCHt(ll)) '_e=1_NoIt=10.mat'];
filename2=['iBDD-SABM_PAM2_PC_BCH_n=' num2str(BCHn) '_t=' num2str(BCHt(ll)) '_e=1_NoIt=10_R=5_SABM_NoIt=5.mat'];


filepath1=[l_path1 filename1];
load(filepath1);
p(1)=semilogy(snrdB,BER,'-s','Color', Color(1,:),'LineWidth',1.5);
hold on;
grid on;

filepath2=[l_path2 filename2];
load(filepath2);
p(2)=semilogy(snrdB,BER,'-s','Color', Color(2,:),'LineWidth',1.5);


%legend(p,'iBDD','SABM','SABM-SR','SABM-SR-ST','TPC');


if tikz
    addpath(genpath('/home/gabriele/Dropbox/Work/Matlab code/matlab2tikz'));
    tikzfilename=['PC_BERvsEbNo_PAM2_BCH_n=255_Variable_t_e=1_NoIt=10_R=5_SABM_NoIt=5.mat.tikz'];
    tpath='/home/gabriele/Dropbox/Work/TUe/MATLAB/MATLAB Sims Data/EnhancedHDdecoding/MarkedSR_PCdecoding/NewResults/tikz/';
    tikzfilepath=[tpath,tikzfilename];
    matlab2tikz(tikzfilepath);
      
end

end
%legend(p,'4 It','5 It');



