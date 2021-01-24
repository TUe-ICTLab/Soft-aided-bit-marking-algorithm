%% Plots results on different bit marking algorithms for product codes and PAM2 in AWGN
clear all;
%close all;
BCHn=127;
ItThresh=5;


if ispc
%l_path='E:\Dropbox\Work\TUe\MATLAB\MATLAB Sims Data\EnhancedHDdecoding\MarkedSR_PCdecoding\NewResults\SABM\';
l_path='E:\Users\Gabriele\Dropbox\Work\TUe\MATLAB\MATLAB Sims Data\EnhancedHDdecoding\MarkedSR_PCdecoding\NewResults\SABM\';
elseif isunix 
l_path='/home/gliga/Dropbox/Work/TUe/MATLAB/MATLAB Sims Data/EnhancedHDdecoding/MarkedSR_PCdecoding/NewResults/SABM/';    
%l_path='/home/gabriele/Dropbox/Work/TUe/MATLAB/MATLAB Sims Data/EnhancedHDdecoding/MarkedSR_PCdecoding/NewResults/SABM/';    

end

Color=rand(length(ItThresh),3);

for tt=1:length(ItThresh)
%filename=['iBDD-SABM_PAM2_PC_BCH_n=' num2str(BCHn) '_t=2_e=1_NoIt=10_R=6_SABM_NoIt=' num2str(ItThresh(tt)) '.mat'];
filename=['iBDD-SABM_PAM2_PC_BCH_n=' num2str(BCHn) '_t=2_e=1_NoIt=10_R=5_SABM_NoIt=' num2str(ItThresh(tt)) '_HighAccuracy.mat'];

filepath=[l_path filename];
load(filepath);

p(tt)=semilogy(snrdB,BER,'-s','Color', Color(tt,:),'LineWidth',1.5);
hold on;
grid on;

end

axis([5,8.5,1e-10,1e-1]);
%legend(p,'1 It','2 It','3 It','4 It','5 It');
%legend(p,'4 It','5 It');