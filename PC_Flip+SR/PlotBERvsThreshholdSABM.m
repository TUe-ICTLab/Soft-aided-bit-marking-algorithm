%% Plot results for BERvsSNR for different marking thresholds
%clear all
close all

if ispc
l_path='E:\Dropbox\Work\TUe\MATLAB\MATLAB Sims Data\EnhancedHDdecoding\MarkedSR_PCdecoding\';
elseif isunix
l_path='/home/gabriele/Dropbox/Work/TUe/MATLAB/MATLAB Sims Data/EnhancedHDdecoding/MarkedSR_PCdecoding/';
end

delta=[4.5:.5:6,10];
Color=[0 1 1; 1 1 0; 1 0 1; 1 0 0; 0 1 0];
figure
for dd=1:length(delta)
   filename=['iBDD-SABM_PAM2_PC_BCH_n=127_t=2_e=1_NoIt=10_R=' num2str(delta(dd)) '.mat' ];
   filepath=[l_path filename];
   load(filepath,'BER','snrdB');
   p(dd)=semilogy(snrdB,BER,'*-','Color',Color(dd,:));
   grid on; hold on; 
end

legend(p,'\delta=4.5','\delta=5','\delta=5.5','\delta=6','\delta=10');
axis([5.9,6.5,1e-7,1e-2]);


