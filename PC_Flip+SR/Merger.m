%% Merge results

%% Script to merge results in /home/uceelig/MATLAB/My Versioned code/EnhancedHDdecoding/Results
M=2;
%l_path='/home/uceelig/MATLAB/My Versioned code/EnhancedHDdecoding/Results/';
l_path='/home/gliga/Dropbox/Work/TUe/MATLAB/MATLAB Sims Data/EnhancedHDdecoding/MarkedSR_PCdecoding/NewResults/SABM-SR/';


    fname='iBDD-SABM-SR_PAM2_PC_BCH_n=127_t=2_e=1_NoIt=10_R=5_SABM_NoIt=5';
    l_filename1=[l_path fname '#1'];
    l_filename2=[l_path fname '#2'];
    %l_filename3=[l_path fname '#3'];
    
    L1=load(l_filename1);
%     L1.snrdB=L1.snrdB(1:end-1);
%     L1.BER=L1.BER(1:end-1);
%     L1.BERout=L1.BERout{1:end-1};
%     L1.BERunc=L1.BERunc(1:end-1);
%     L1.BERout_unc=L1.BERout_unc{1:end-1};
%     L1.SERunc=L1.SERunc(1:end-1);
%     L1.SERout_unc=L1.SERout_unc{1:end-1};
    
    
    L2=load(l_filename2);
   % L3=load(l_filename3);
    
    L4=[];
    %L4.snrdB=[L1.snrdB L2.snrdB L3.snrdB];
    L4.snrdB=[L1.snrdB L2.snrdB];

    %L4.BER=[L1.BER L2.BER L3.BER];
    L4.BER=[L1.BER L2.BER];

    %L3.snrdB=[L1.snrdB L2.snrdB];
    %[L4.snrdB, idx]=sort(L3.snrdB);
    %[L4.snrdB,idx2]=unique(L4.snrdB);
    %idx=idx(idx2);
    
    %L4.BER=[L1.BER L2.BER,L3.BER];
    
    %L4.BER=[L1.BER L2.BER];

    %L4.BER=L4.BER(idx);
    
    %L4.BERout=[L1.BERout L2.BERout L3.BERout];
    L4.BERout=[L1.BERout L2.BERout];

    %L4.BERout=L4.BERout{idx};
    
    %L4.BERout_unc=[L1.BERout_unc L2.BERout_unc L3.BERout_unc];
    L4.BERout_unc=[L1.BERout_unc L2.BERout_unc];

   % L4.BERout_unc=L4.BERout_unc{idx};
    
    
    %L4.BERunc=[L1.BERunc L2.BERunc L3.BERunc];
    L4.BERunc=[L1.BERunc L2.BERunc];

    %L4.BERunc=L4.BERunc(idx);
    
    %L4.SERout_unc=[L1.SERout_unc L2.SERout_unc L3.SERout_unc];
    L4.SERout_unc=[L1.SERout_unc L2.SERout_unc];

    %L4.SERout_unc=L4.SERout_unc{idx};
    
    
    %L4.SERunc=[L1.SERunc L2.SERunc L3.SERunc];
    L4.SERunc=[L1.SERunc L2.SERunc];

    %L4.SERunc=L4.SERunc(idx);
    
    L4.r=L1.r;
    s_filename=[l_filename1 '_merger'];
    save(s_filename, '-struct','L4');
    