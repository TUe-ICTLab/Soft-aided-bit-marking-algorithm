%% Plots results on different bit marking algorithms for product codes and PAM2 in AWGN

if ispc
l_path='E:\Dropbox\Work\TUe\MATLAB\MATLAB Sims Data\EnhancedHDdecoding\MarkedSR_PCdecoding\';
elseif isunix 
l_path='/home/gabriele/Dropbox/Work/TUe/MATLAB/MATLAB Sims Data/EnhancedHDdecoding/MarkedSR_PCdecoding/';    
end

%% Component BCH (7,2,1)

filename1='iBDD_PAM2_PC_BCH_n=127_t=2_e=1_NoIt=10.mat';
filepath1=[l_path filename1];
filename2='iBDD-SABM_PAM2_PC_BCH_n=127_t=2_e=1_NoIt=10_R=5_merger.mat';
filepath2=[l_path filename2];
filename3='iBDD-SABM-SR_PAM2_PC_BCH_n=127_t=2_e=1_NoIt=10_MarkingIter= 4.5_R=5.mat';
filepath3=[l_path filename3];
filename4='iBDD-SABM-SR-ST_PAM2_PC_BCH_n=127_t=2_e=1_NoIt=10_MarkingIter= 4.5_R=5.mat';
filepath4=[l_path filename4];

figure;
load(filepath1);
p1=semilogy(snrdB,BER,'sr-');
hold on;
grid on;

load(filepath2);
p2=semilogy(snrdB,BER,'sb-');

load(filepath3);
p3=semilogy(snrdB,BER,'sk-');

load(filepath4);
p4=semilogy(snrdB,BER,'sg-');


%% Component BCH (8,2,1)
filename1='iBDD_PAM2_PC_BCH_n=255_t=2_e=1_NoIt=10.mat';
filepath1=[l_path filename1];
filename2='iBDD-SABM_PAM2_PC_BCH_n=255_t=2_e=1_NoIt=10_R=5.mat';
filepath2=[l_path filename2];
filename3='iBDD-SABM-SR_PAM2_PC_BCH_n=255_t=2_e=1_NoIt=10_MarkingIter= 4.5_R=5.mat';
filepath3=[l_path filename3];
filename4='iBDD-SABM-SR-ST_PAM2_PC_BCH_n=255_t=2_e=1_NoIt=10_MarkingIter= 4.5_R=5.mat';
filepath4=[l_path filename4];

load(filepath1);
semilogy(snrdB,BER,'sr-');
hold on;
grid on;

load(filepath2);
semilogy(snrdB,BER,'sb-');

load(filepath3);
semilogy(snrdB,BER,'sk-');

load(filepath4);
p4=semilogy(snrdB,BER,'sg-');

axis([5,8.5,1e-10,1e-1]);
xlabel('SNR [dB]');
ylabel('BER');
legend([p1,p2,p3,p4],'iBDD','iBDD+SABM','iBDD+SABM-SR','iBDD+SABM-SR-ST');
txtfilename='Results.tikz';
matlab2tikz(txtfilename);
