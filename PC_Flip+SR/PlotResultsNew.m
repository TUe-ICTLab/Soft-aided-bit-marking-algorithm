%% Plots results on different bit marking algorithms for product codes and PAM2 in AWGN
clear all;
close all;

BCHv=[7,8];
BCHn=2.^BCHv-1;
ItThresh=5;
tikz=0;

if ispc
%l_path='E:\Dropbox\Work\TUe\MATLAB\MATLAB Sims Data\EnhancedHDdecoding\MarkedSR_PCdecoding\NewResults\SABM\';
l_path1='E:\Users\Gabriele\Dropbox\Work\TUe\MATLAB\MATLAB Sims Data\EnhancedHDdecoding\MarkedSR_PCdecoding\NewResults\SABM\';

l_path2='E:\Users\Gabriele\Dropbox\Work\TUe\MATLAB\MATLAB Sims Data\EnhancedHDdecoding\MarkedSR_PCdecoding\NewResults\SABM-SR\';


elseif isunix 
%l_path1='/home/gabriele/Dropbox/Work/TUe/MATLAB/MATLAB Sims Data/EnhancedHDdecoding/MarkedSR_PCdecoding/NewResults/SABM/';
l_path1='/home/gliga/Dropbox/Work/TUe/MATLAB/MATLAB Sims Data/EnhancedHDdecoding/MarkedSR_PCdecoding/NewResults/SABM/';

%l_path2='/home/gabriele/Dropbox/Work/TUe/MATLAB/MATLAB Sims Data/EnhancedHDdecoding/MarkedSR_PCdecoding/NewResults/SABM-SR/';
l_path2='/home/gliga/Dropbox/Work/TUe/MATLAB/MATLAB Sims Data/EnhancedHDdecoding/MarkedSR_PCdecoding/NewResults/SABM-SR/';    
%l_path3='/home/gabriele/Dropbox/Work/TUe/MATLAB/MATLAB Sims Data/EnhancedHDdecoding/MarkedSR_PCdecoding/NewResults/';    
l_path3='/home/gliga/Dropbox/Work/TUe/MATLAB/MATLAB Sims Data/EnhancedHDdecoding/MarkedSR_PCdecoding/NewResults/';    


%l_path1='/home/gabriele/Dropbox/Work/TUe/MATLAB/MATLAB Sims Data/EnhancedHDdecoding/MarkedSR_PCdecoding/NewResults/SABM/';
%l_path2='/home/gabriele/Dropbox/Work/TUe/MATLAB/MATLAB Sims Data/EnhancedHDdecoding/MarkedSR_PCdecoding/NewResults/SABM-SR/';  
%l_path='/home/gabriele/Dropbox/Work/TUe/MATLAB/MATLAB Sims Data/EnhancedHDdecoding/MarkedSR_PCdecoding/NewResults/SABM/';    

end

Color=[0 0 1; 1 0 1 ; 1 0.6 0.7 ; 0 1 1];


for ll=1:length(BCHn)
figure;

filename1=['iBDD_PAM2_PC_BCH_n=' num2str(BCHn(ll)) '_t=2_e=1_NoIt=10.mat'];
filename2=['iBDD-SABM_PAM2_PC_BCH_n=' num2str(BCHn(ll)) '_t=2_e=1_NoIt=10_R=5_SABM_NoIt=5.mat'];
filename3=['iBDD-SABM-SR_PAM2_PC_BCH_n=' num2str(BCHn(ll)) '_t=2_e=1_NoIt=10_R=5_SABM_NoIt=5.mat'];
filename4=['iBDD-SABM-SR-ST_PAM2_PC_BCH_n=' num2str(BCHn(ll)) '_t=2_e=1_NoIt=10_R=5_SABM_NoIt=5.mat'];
filename5=['TPC_Chase-Pyndiah_PAM2_PC_BCH_n=' num2str(BCHn(ll)+1) '_t=2_e= 1_NoIt=10.mat'];


filepath1=[l_path1 filename1];
load(filepath1);
EbNo=snrdB-3-10*(log10(r));
%ber1=BER;
%Dist=ClosedFormEGNModel(EbNo,r);
p(1)=semilogy(EbNo,BER,'-s','Color', Color(1,:),'LineWidth',1.5);
hold on;
grid on;

filepath2=[l_path1 filename2];
load(filepath2);
EbNo=snrdB-3-10*(log10(r));
p(2)=semilogy(EbNo,BER,'-s','Color', Color(2,:),'LineWidth',1.5);


filepath3=[l_path2 filename3];
load(filepath3);
EbNo=snrdB-3-10*(log10(r));
p(3)=semilogy(EbNo,BER,'-s','Color', Color(3,:),'LineWidth',1.5);


filepath4=[l_path2 filename4];
load(filepath4);
EbNo=snrdB-3-10*(log10(r));
p(4)=semilogy(EbNo,BER,'-s','Color', Color(4,:),'LineWidth',1.5);

if BCHn(ll)==127
filepath5=[l_path3 filename5];
load(filepath5);
p(5)=semilogy(EbNo(1:end-1),BER(1:end-1),'--k','LineWidth',1.5);
end


if BCHn(ll)==127
axis([2.5,5,1e-10,1e-1]);
else
axis([3.5,5.5,1e-10,1e-1]);
end

legend(p,'iBDD','SABM','SABM-SR','SABM-SR-ST','TPC');


if tikz
    addpath(genpath('/home/gliga/Dropbox/Work/Matlab code/matlab2tikz-matlab2tikz-816f875'));
    tikzfilename=['PC_BERvsEbNo_PAM2_BCH_n=' num2str(BCHn(ll)) '_t=2_e=1_NoIt=10_R=5_SABM_NoIt=5.mat.tikz'];
    tpath='/home/gliga/Dropbox/Work/TUe/MATLAB/MATLAB Sims Data/EnhancedHDdecoding/MarkedSR_PCdecoding/NewResults/tikz/';
    tikzfilepath=[tpath,tikzfilename];
    matlab2tikz(tikzfilepath);
      
end

  %close all;
  
    %b1=semilogy(Dist,ber1);
  if tikz
    tikzfilename=['DistMapEbNo_BCHn=' num2str(BCHn(ll)) 'test.tikz'];
    tpath='/home/gliga/Dropbox/Apps/Overleaf/ECOC2019_A novel soft-aided marking decoder for product codes/Figures/';
    tikzfilepath=[tpath,tikzfilename];
    matlab2tikz(tikzfilepath);
  end


end
%legend(p,'4 It','5 It');



