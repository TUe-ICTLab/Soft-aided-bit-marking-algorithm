function [error_vec]=Produ_IMP_cod_BCH_AWGN_chase_pyndiah_matlab_v2(iter,SNR,v,t,s,modord,FrameErrStop)

% Chase Pyndiah decoder based on [1]
% BCH component code
% Alireza Sheikh,  asheikh@chalmers.se, May 2018
% [1] Near-Optimum Decoding of Product Codes: Block Turbo Codes
%%
% Note1: "beta_vec" and "alpha_vec" lengthes are 2*iter, correspondong 2*iter half iterations
% Note2: beta_vec=[0.2 0.4 0.6 0.8 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1], see [1]
% Note3: alpha_vec=[0.2 0.3 0.5 0.7 0.9 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1], see [1] 
%%
%clc
%tic
%% 
%%%%%%% ======================== Parity check (Reed Solomon) ========================= %%%%%%%%

EbN0=SNR;
tt=t;

rng('shuffle')
%rng(600)
%rng(1500)

%%%%%%%%%%========== Initiaization ==========%%%%%%%%%%%%%
%%
block_len_original=2^(v)-1;       % original block length
inf_len_original=2^(v)-v*tt-1;    % information length
block_len=2^(v)-1-s;              % shortened original block length
inf_len=2^(v)-v*tt-1-s;           % shortened information length
kk=inf_len;
nn=block_len;
%%
%%%%%%%%% ================ Sigma generation ==========================%%%%%

Rate_c=(inf_len/(block_len+1))^2;
Es_N0=((10^(EbN0/10))*modord)*Rate_c; 
sigma=sqrt(1/(2*Es_N0));

%% 
%%%%%%%%%%% ======================== Transmitt and receive loop ========================= %%%%%%%%

error_vec=zeros(1,length(SNR));  % Vector containing number of total errors per SNR 
error_tot=zeros(1,length(SNR));  % Vector containing number of total errors per SNR 
FrameErrTot=zeros(1,length(SNR));
%FrameErr=zeros(1,length(SNR));
Nblocks_min=10*FrameErrStop;
%Nblocks_min=1;

% %% PARFOR LOOP
MatWork=32;
 if isempty(gcp('nocreate'))  
 Par1=parpool(MatWork);   
 %Par1=parpool('local');   

 else
  Par1=gcp('nocreate');
 end

%error_vec1=zeros(1,length(SNR));

for kk=1:length(SNR)

Eb=SNR(kk); 

error=zeros(1,Nblocks_min);    % Vector containing error per decoded block
block_count=0;   % 
while FrameErrTot(kk)<FrameErrStop

disp(['Decoding Blocks ' num2str(block_count) ' to ' num2str(block_count+Nblocks_min) '...']);
FrameErr=0;    
parfor indt=1:Nblocks_min   
    
%%
%%%%%%%% ================= Encoder & Modulation Initialization ==================== %%%%%%%%
%enc = comm.BCHEncoder(block_len_original,inf_len_original);
%%%%%%%%%%% ============================= Encoder ================================== %%%%%%%%
%%
%inf_stair_block=zeros(inf_len,inf_len);            % Inf. staircase block
%encod_stair_block1=zeros(inf_len,block_len);       % Block staircase
%encod_stair_block=zeros(block_len,block_len);      % Block staircase
%%%%%%%%%% ============================= Inf block =========================== %%%%%%%%%%%%%%
%%
inf_stair_block = randi([0 1],inf_len*inf_len,1);     % Data symbols
encod_stair_block11=tpcenc(inf_stair_block,[block_len+1;block_len+1],[inf_len;inf_len]);
encod_stair_block=reshape(encod_stair_block11,block_len+1,block_len+1);
%%%%%%%%% ========================== Encode Inf block ============================ %%%%%%%%%%
%%

%%%% Interleaving
state=randi(100000000);
encod_stair_block1=randintrlv(encod_stair_block,state);


senddata=encod_stair_block1(:);
l_l=length(senddata);
temp=reshape(senddata,modord,l_l/modord);
senddata_t=symbol_binlabels(temp.',2^modord);

%% Channel
%%%%%%% ====================== AWGN channel channel ==================%%%%%%%%

yy_demod=senddata_t;
yy_demodqam=real(pammod(yy_demod,2^(modord)));

% %%% Creat labels for the information symbols
% yyindex=(yy_demodqam+2^(modord)-1)/2;
% yyindex_gray=bin2gray(yyindex,'pam',2^(modord));
% uncoded_bits=de2bi(yyindex_gray',modord);
% uncoded_bits1=uncoded_bits.';
% uncoded_bits2=uncoded_bits1(:);
% uncoded_bits3=reshape(uncoded_bits2,block_len,block_len);

data_modu_snr=yy_demodqam;
noise=sigma*(randn(length(yy_demodqam),1)+sqrt(-1)*randn(length(yy_demodqam),1));
receive_vec=real(data_modu_snr+noise);

%%% ============ Reliability ============ %%%
%%
reliability_vec=(receive_vec);

%%%%%%%%%%%%% ====================  Demodulation ================== %%%%%%%%%%%%%

%[region,label_out,opt_deta]=delta_out(ASK_ord,SNR);
%decode_out=zeros(length(receive_vec),modord);

decode_out=pamdemod(receive_vec,2^modord);
decode_out1=bin2gray(decode_out,'pam',2^modord);
decode_out2=de2bi(decode_out1);
decode_out3=decode_out2';
decode_out4=decode_out3(:);

%decodeblock=reshap(decode_out4,block_len,block_len);
%y_demodp=gray2bin(y_demod1,'qam',2^(modord));
%length(find(uncoded_bits~=decode_out))
%length(find(decode_out==2))

%%%%% ==================================================================== %%%%%%

% y_mod2vv=decode_out;
% y_mod2pvv=y_mod2vv';
% y_demod11=y_mod2pvv(:);

%% 
%%%%%%%%%%% =========================== Decoder ============================= %%%%%%%%

%%%% Deinterleaver should come here

y_demod=reshape(decode_out4,[block_len+1,block_len+1]);
reliability_demod=reshape(reliability_vec,[block_len+1,block_len+1]);

%%%% Deinterleaving

y_demod=randdeintrlv(y_demod,state);
reliability_demod=randdeintrlv(reliability_demod,state);

%y_decode_gf_inf=y_decode_func(y_demod,reliability_demod,block_len,inf_len,iter,v,s,beta_vec,alpha_vec);

y_decode_gf_inf1=tpcdec(reliability_demod(:),[block_len+1;block_len+1],[inf_len;inf_len],[],iter);
y_decode_gf_inf=reshape(y_decode_gf_inf1,inf_len,inf_len);

% length(find(y_demod~=encod_stair_block))           % error before decoding (block error)
% length(find(y_decode_gf_inf~=encod_stair_block))   % error after decoding (block error)

%%%%%%%%%% ================= Remove transient ============================ %%%%%%%%%%%%
%%
y_decode_gf_inft=y_decode_gf_inf(1:inf_len,1:inf_len);
inf_stair_blockt=encod_stair_block(1:inf_len,1:inf_len);

%%%%%%%% ============================ Inf. bit error rate ========================== %%%%%%%
 
vvabit=find(y_decode_gf_inft~=inf_stair_blockt);

%Eb
error(indt)=length(vvabit);   % Errors per block 
if ~isempty(vvabit)
FrameErr=FrameErr+1;  
end
block_count=block_count+1;     % Count decoded blocks
end
 
%%% Inf. Symbol error rate
%%

% vvabit=find(senddata_tt~=senddata_tt1);
error_tot(kk)=error_tot(kk)+sum(error);    % Total errors  
FrameErrTot(kk)=FrameErrTot(kk)+FrameErr;
disp([num2str(FrameErrTot(kk)) ' Block Errors']);
%disp(block_count);
end

% [az,bz]=size(senddata_tt);
% error_vec1(kk)=error1/(mo*az*bz);

%[az,bz]=size(y_decode_gf_inft);
error_vec(kk)=error_tot(kk)/(block_count*inf_len^2);

end
%% 
%%%%%%%%%%% ============================= Plot =============================== %%%%%%%%%%%%
% hold on
% semilogy(Eb,error_vec,'*k');
% xlabel('E_b/N_0(dB)')
% ylabel('BER')
%time=toc;
%time=0;
end

function x2=symbol_binlabels(x_in,pamord)
 
x1=bi2de(x_in);
x2=gray2bin(x1,'pam',pamord);

end