function [PosFECBERC]=merge_incomplete_results(R,R_str,m,BCH,EsN0dB)
% merge_incomplete_results(R_str,m,BCH,EsN0dB)
% This function merges the results saved in
% /results/R/MPAM/incomplete/ and save them in 
% /results/R/MPAM/
% The function should be used in cases where not enough bit errors were
% counted in the simulation and temporary files were saved. 
% 
% Alex Alvarado
% March 2018

pref_save_dir_in   = strcat('results/',R_str,'/',num2str(2^m),'PAM/incomplete/');
res_name=strcat(pref_save_dir_in,num2str(2^m),'PAM_',num2str(EsN0dB,4),'_dB_m_',num2str(BCH.m),'_t_',num2str(BCH.t),'_w_',num2str(BCH.w),'_iter_',num2str(BCH.iter),'*.mat');
list = dir(res_name);
Nfiles=size(list,1);
if Nfiles==0,fprintf('No files to merge!\n');return
else fprintf('Merging %i files...\n',Nfiles);end

% Do the merging
PreFECBERC=0;
PosFECBERC=0;

FERC=0;
MIav=0;
GMIav=0;
HDMIav=0;
total_info_bits=0;
total_coded_bits=0;
sentblocks=0;
cd(pref_save_dir_in)
for i=1:Nfiles
    %list(i).name
    file=load(list(i).name);
    
    PreFECBERC          = file.PreFECBERC+PreFECBERC;
    PosFECBERC          = file.PosFECBERC+PosFECBERC;
    FERC                = file.FERC+FERC;
    
    MIav                = file.MIav+MIav;
    GMIav               = file.GMIav+GMIav;
    HDMIav               = file.HDMIav+HDMIav;
    
    total_info_bits = file.total_info_bits+total_info_bits;
    total_coded_bits= file.total_coded_bits+total_coded_bits;
    sentblocks          = file.sentblocks+sentblocks;
    
end

%total_coded_bits= sentblocks*n;                     % Number of sent coded bits
%PreFECBER       = PreFECBERC/n;                     % Pre FEC BER
%total_info_bits = sentblocks*k;                     % Number of sent info bits
%PosFECBER       = PosFECBERC/k;                     % Post FEC BER

PreFECBERav     = sum(PreFECBERC)/total_coded_bits; % Average pre FEC BER
PosFECBERav     = sum(PosFECBERC)/total_info_bits;  % Average Post FEC BER
FERav           = FERC/sentblocks;                  % FER

MIav            = MIav/Nfiles;                         % Average MI
GMIav           = GMIav/Nfiles;                      % Average GMI
HDMIav          = HDMIav/Nfiles;                      % Average HDMI

% Display results
% % fprintf('\n-------------------------------------------------\n');
% % fprintf('Inf. Bit errors = %d \n',sum(PosFECBERC));
% % fprintf('Coded Bit errors = %d \n',sum(PreFECBERC));
% % fprintf('Frame errors = %d \n',FERC);
% % fprintf('MI = %d \n',MIav);
% % fprintf('GMI = %d \n',GMIav);
% % fprintf('HDMI = %d \n',HDMIav);
% % fprintf('Post FEC BER  = %d\n',PosFECBERav);
% % fprintf('Pre FEC BER = %d \n',PreFECBERav);
% % fprintf('FER  = %1.6f\n',FERav);
% % fprintf('-------------------------------------------------\n');

cd ../../../../

pref_save_dir   = strcat('results/',R_str,'/',num2str(2^m),'PAM/');
res_name=[pref_save_dir,num2str(2^m),'PAM_',num2str(EsN0dB,4),'_dB_m_',num2str(BCH.m),'_t_',num2str(BCH.t),'_w_',num2str(BCH.w),'_iter_',num2str(BCH.iter),'.mat'];

% Set a flag for the merged 
merged_flag=1;
% Save
X=file.X;L=file.L;M=file.M;m=file.m;EbN0dB=file.EbN0dB;EbN0=file.EbN0;EsN0=file.EsN0;EsN0dB=file.EsN0dB;k=file.k;kwoLast=file.kwoLast;
   

save(res_name,'X','L','M','m','k','kwoLast','BCH',...
        'PreFECBERav','PosFECBERav',...
        'FERav','FERC','EbN0','total_coded_bits','total_info_bits','sentblocks',...
        'EbN0dB','EsN0','EsN0dB','R','MIav','GMIav','HDMIav','merged_flag');

return
