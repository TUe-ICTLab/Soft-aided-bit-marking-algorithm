%% Script simulating the performance of a product code over an AWGN channel
% Implements standard iBDD decoding and iBDD decoding using scaled
% reliabilities (from LLR) as in [1]. 

% Note that the BCH decoder is compiled for windows so it will only work in
% a windows machine...

%% References 
% [1] Yi L. et al. "Decoding Staircase Codes with Marked Bits", International Symposium on Turbo Codes & Iterative Information Processing 2018
startup;
clear all;
Save=0;

BCHnu=7;
BCHt=2;
BCHe=1;
BCHn=(2^BCHnu)-1;                   % Coded bits, without the extended bits
BCHk=(2^BCHnu)-1-BCHnu*BCHt;        % Information bits
I=10; % decoder iterations
Nb=(BCHn+BCHe)^2;                   % Number of bits per PC block
Mark=1;

if Mark
R=10;  % LLR reliability threshold for bit marking for decoding algorithm in [1]  
end

r=BCHk^2/Nb;                               % PC rate

% Modulation parameters 
M=2;                   % Constellation cardinality/per dimension 
m=log2(M);             % Bits/symbol

N=1;                   % Constellation dimensions
X=-(M-1):2:M-1;      % PAM constellation
Ns=ceil(Nb/m);
Es=mean(abs(X).^2);
X=X/sqrt(Es);   
Mode='soft';

flag_demapper=1;       % Max-log algorithm for LLR calculation
L=de2bi(0:M-1,m);

% Sweep parameters 
%snrdB=12.5:0.1:15;
%snrdB=18:.1:18.9;
snrdB=5.5:.1:7;
%snrdB=18:0.1:22.5;
ErrStop=1e3;

%Encoder/Decoder Parameters
Nblocks_min=ErrStop/10;
Nblocks_max=1e6;

% Bits errors per block after FEC
UncSymErr=cell(1,length(snrdB));
UncBitErr=cell(1,length(snrdB));
BitErr=cell(1,length(snrdB));

% Total errors per after FEC
UncBitErrTot=zeros(1,length(snrdB));
UncSymErrTot=zeros(1,length(snrdB));
BitErrTot=zeros(1,length(snrdB));

% Error rates per block
BERout_unc=cell(1,length(snrdB));
SERout_unc=cell(1,length(snrdB));
BERout=cell(1,length(snrdB));

% Average error rates
BERunc=zeros(1,length(snrdB));
SERunc=zeros(1,length(snrdB));
BER=zeros(1,length(snrdB));


% PARFOR LOOP
%MatWork=32;
%  if isempty(gcp('nocreate'))  
%  Par1=parpool('local');   
%  else
%      Par1=gcp('nocreate');
%  end

for ss=1:length(snrdB)
    %fprintf('x');
    k=0;
    
    % Initialize error cells
    UncBitErr{ss}=[];
    UncSymErr{ss}=[];
    BitErr{ss}=[];
    BERout_unc{ss}=[];
    BERout{ss}=[];
    snr=snrdB(ss);      % copy into broadcast variable for parfor loop
    
    while BitErrTot(ss)<ErrStop && (k+1)*Nblocks_min<Nblocks_max
      sym_err_block=zeros(1,Nblocks_min);
      err_block_unc=zeros(1,Nblocks_min);
      err_block=zeros(1,Nblocks_min);
      ber_block_unc=zeros(1,Nblocks_min);
      ser_block_unc=zeros(1,Nblocks_min);
      ber_block=zeros(1,Nblocks_min);
      disp(['Decoding Blocks ' num2str(k*Nblocks_min) ' to ' num2str((k+1)*Nblocks_min) '...']);
    for ii=1:Nblocks_min
        
        info_bits=randi([0,1],BCHk,BCHk);                               % Generate the info bit matrix
        
        %% PC encoding
        coded_bits=PC_BCH_Encoder(info_bits,BCHnu,BCHt,BCHe);           % Encode the data using a product code
        coded_bits_res=reshape(coded_bits,1,Nb);
        
       % Bit-to-symbol mapping 
       %info_bits_res=reshape(info_bits,1,BCHk^2);
       %sym=mapper(X.',info_bits_res,m,BCHk^2/m,N).';
       pad=mod(length(coded_bits_res),m);
      
       % Padding for number of bits in the block non-mulitple of m
        if ~isequal(pad,0)
         coded_bits_res=[coded_bits_res zeros(1,m-pad)];
        end
        
        
        sym=mapper(X.',coded_bits_res,m,Ns,N).';
        %sym=sym/sqrt(Es);                             % Normalise constellation to unitary energy
        y=awgn(sym,snr).';                             % Transmit through a AWGN channel
        
        %info_bits_awgn=zeros(BCHk);        
         
         %% HARD DEMAPPING
         hat_y_dec=pamdemod(y*sqrt(Es),M);
         info_bits_awgn=de2bi(hat_y_dec).';
         if ~isequal(pad,0)
         info_bits_awgn=info_bits_awgn(1:end-m+pad);
         end
         info_bits_awgn=reshape(info_bits_awgn,BCHn+BCHe,BCHn+BCHe); 
         hat_y=pammod(hat_y_dec,M);
                 
        %% SOFT symbol-to-bits demapping (LLR calculator)
        l=demapper(X.',L,snr,y,flag_demapper); 
        % histogram(l,100);
        
        % Bits marking (algorithm used in [1])
        %l=l(:,1:end-m+pad);
        if ~isequal(pad,0)
        l=reshape(l(1:end-m+pad),BCHn+BCHe,BCHn+BCHe);
        else
        l=reshape(l,BCHn+BCHe,BCHn+BCHe);
        end
        
        %l=reshape(l,1,BCHk^2);
        
        %HARD-DECISION on bits
        idx1=(l>0);
        
        % HARD-DECISION on symbols
        [sym_err_block(ii),ser_block_unc(ii)]=symerr(sym*sqrt(Es),hat_y.');
        
        coded_bits_awgn=zeros(BCHn+BCHe);
        coded_bits_awgn(idx1)=1;
                 
        %TP decoding 
        if Mark 
        marked_bits=(abs(l)>R);   % Marks bits based on LLRs
        hat_coded_bits=PC_BCH_Marked_Decoder(coded_bits_awgn,marked_bits,BCHnu,BCHt,BCHe,I);    % Decode
        else
        hat_coded_bits=PC_BCH_Decoder(coded_bits_awgn,BCHnu,BCHt,BCHe,I);% Decode
        end
        err_block_unc(ii)=sum(sum(coded_bits_awgn(1:BCHk,1:BCHk)~=info_bits));
        err_block(ii)=sum(sum(hat_coded_bits(1:BCHk,1:BCHk)~=info_bits));
        ber_block_unc(ii)=err_block_unc(ii)/(BCHk^2); 
        ber_block(ii)=err_block(ii)/(BCHk^2); 
        %Err{ss}(ii)=sum(sum(hat_coded_bits(1:BCHk,1:BCHk)~=info_bits));
        %BERout{ss}(ii)=Err{ss}(k*Nblocks_min+ii)/(BCHk^2); % Count errors
    
     end
    
    % Errors in each block
    UncSymErr{ss}=[UncSymErr{ss},sym_err_block];  
    UncBitErr{ss}=[UncBitErr{ss},err_block_unc]; 
    BitErr{ss}=[BitErr{ss},err_block];
    
    % Total errors up to block k*Nblock_min 
    UncBitErrTot(ss)=sum(UncBitErr{ss});
    UncSymErrTot(ss)=sum(UncSymErr{ss});
    BitErrTot(ss)=sum(BitErr{ss});
    
    % Error rates in each block
    BERout_unc{ss}=[BERout_unc{ss},ber_block_unc];
    SERout_unc{ss}=[SERout_unc{ss},ser_block_unc];
    BERout{ss}=[BERout{ss},ber_block];
    
    k=k+1;
    end
   
    BERunc(ss)=UncBitErrTot(ss)/(length(UncBitErr{ss})*Nb); 
    SERunc(ss)=UncBitErrTot(ss)/(length(UncBitErr{ss})*Ns);
    BER(ss)=BitErrTot(ss)/(length(BitErr{ss})*Nb); 
    fprintf('SNR=%4.2f dB, Uncoded BER=%e \n',snrdB(ss),BERunc(ss));
    fprintf('SNR=%4.2f dB, Uncoded SER=%e \n',snrdB(ss),SERunc(ss));
    fprintf('SNR=%4.2f dB, BER=%e \n',snrdB(ss),BER(ss));
   
end

% Plot results
figure(1);
semilogy(snrdB,SERunc,'go-'); hold on;
semilogy(snrdB,BERunc,'b*-');
semilogy(snrdB,BER,'r*-');
grid on;hold on;
axis([12,15,1e-10,1e-1]);

if Save 
s_path='/home/uceelig/MATLAB/My Versioned code/EnhancedHDdecoding/Results/';
if ~exist(s_path,'dir')
mkdir(s_path); 
end
disp('Saving...');
filename=['PAM' num2str(M) '_PC_BCH_n=' num2str(BCHn) '_t=' num2str(BCHt) '_e=' num2str(BCHe) '_NoIt=' num2str(I) '_marked_R=' num2str(R)  '.mat'] ;
save([s_path filename],'snrdB','BERout_unc','SERout_unc','BERout','BERunc','SERunc','BER','r','R');
end

