%% Script simulating the performance of product codes over an AWGN channel
% Implements iBDD decoding using scaled
% reliabilities (from LLR) and bit marking as in [1-2]. 



%% References 
% [1] Yi L. et al. "Decoding Staircase Codes with Marked Bits", International Symposium on Turbo Codes & Iterative Information Processing 2018
% [2] Sheikh, A., Amat, A. G. i, & Liva, G. (2019). Binary Message Passing Decoding of Product Codes Based on Generalized Minimum Distance Decoding. ArXiv:1901.02914. Retrieved from http://arxiv.org/abs/1901.02914
startup;

% Flags 
Save=0;
DecoderType='SABM';   % Options are: 'iBDD', 'iBDD-SR', 'SABM', 'SABM-SR', 'SABM-SR-ST'  
Inter=0;               % Enables interleaver 

BCHnu=7;
BCHt=2;
BCHe=1;
BCHn=(2^BCHnu)-1;                   % Coded bits, without the extended bits
BCHk=(2^BCHnu)-1-BCHnu*BCHt;        % Information bits
I=10;                               % iBDD decoder iterations
Nb=(BCHn+BCHe)^2;                   % Number of bits per PC block
%MarkSR=1;


if isequal(DecoderType, 'iBDD-SR') || isequal(DecoderType, 'SABM-SR') || isequal(DecoderType, 'SABM-SR-ST') 
w=[3.427276299321742   3.866157771295065   4.083911657040851 ...
    4.271760664386480   4.490522576596476   4.817956305103627 5.453885179993800 ...
    7.083463707657454  11.641714059713033  23.369586755750639];                        % Vector of weigths for SR update for iBDD-SR algorithm
end

if isequal(DecoderType, 'SABM')|| isequal(DecoderType, 'SABM-SR') || isequal(DecoderType, 'SABM-SR-ST') 
R=5;  % LLR reliability threshold for bit marking for SABM
ItThresh=5;
end


r=BCHk^2/Nb;                               % PC rate

% Modulation parameters 
M=2;                   % Constellation cardinality/per dimension 
m=log2(M);             % Bits/symbol

N=1;                   % Constellation dimensions
X=-(M-1):2:M-1;        % PAM constellation
Ns=ceil(Nb/m);
Es=mean(abs(X).^2);
X=X/sqrt(Es);   
Mode='soft';

flag_demapper=1;       % Max-log algorithm for LLR calculation
L=de2bi(0:M-1,m);

% Sweep parameters 
%snrdB=12.5:0.1:13.8;
%snrdB=18:.1:18.9;
snrdB=5:.1:6.2;
%snrdB=18.5:0.1:19.6;
FrameErrStop=5e1;          % Number of frame errors at which MC simulation stops

%Encoder/Decoder Parameters
Nblocks_min=2*FrameErrStop;
Nblocks_max=3e6;

% Bits errors per block after FEC
UncSymErr=cell(1,length(snrdB));
UncBitErr=cell(1,length(snrdB));
BitErr=cell(1,length(snrdB));

% Total errors per after FEC
UncBitErrTot=zeros(1,length(snrdB));
UncSymErrTot=zeros(1,length(snrdB));
BitErrTot=zeros(1,length(snrdB));
FrameErrTot=zeros(1,length(snrdB));

% Error rates per block
BERout_unc=cell(1,length(snrdB));
SERout_unc=cell(1,length(snrdB));
BERout=cell(1,length(snrdB));

% Average error rates
BERunc=zeros(1,length(snrdB));
SERunc=zeros(1,length(snrdB));
BER=zeros(1,length(snrdB));

% %% PARFOR LOOP
MatWork=16;
 if isempty(gcp('nocreate'))  
 %Par1=parpool(MatWork);   
 Par1=parpool('local');   

 else
     Par1=gcp('nocreate');
 end

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
    
    while FrameErrTot(ss)<FrameErrStop && (k+1)*Nblocks_min<Nblocks_max
      sym_err_block=zeros(1,Nblocks_min);
      err_block_unc=zeros(1,Nblocks_min);
      err_block=zeros(1,Nblocks_min);
      ber_block_unc=zeros(1,Nblocks_min);
      ser_block_unc=zeros(1,Nblocks_min);
      ber_block=zeros(1,Nblocks_min);
      disp(['Decoding Blocks ' num2str(k*Nblocks_min) ' to ' num2str((k+1)*Nblocks_min) '...']);
      FrameErr=0;         % Frame errors per chunk of code blocks decoded in parallel 

           
     parfor ii=1:Nblocks_min
        
        info_bits=randi([0,1],BCHk,BCHk);                               % Generate the info bit matrix
        
        %% PC encoding
        coded_bits=PC_BCH_Encoder(info_bits,BCHnu,BCHt,BCHe);           % Encode the data using a product code    
        coded_bits_res=reshape(coded_bits,1,Nb);
        
        % Bit interleaver
         if Inter
        [coded_bits_res,map]=RandInter(coded_bits_res);
         end
        
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
         hat_y=pammod(hat_y_dec,M);
         info_bits_awgn=de2bi(hat_y_dec).';
         if ~isequal(pad,0)
         info_bits_awgn=info_bits_awgn(1:end-m+pad);
         end
         info_bits_awgn=reshape(info_bits_awgn,BCHn+BCHe,BCHn+BCHe); 
                 
        %% SOFT symbol-to-bits demapping (LLR calculator)
        l=demapper(X.',L,snr-3,y,flag_demapper); % snr is converted to Es/No as per definition in demapper function
        % histogram(l,100);
        
        
        % Bit de-interleaver
         
         if ~isequal(pad,0) 
          l=l(1:end-m+pad);
         end
         
         if Inter
          l=RandInter(l,map);    
         end
                   
        l=reshape(l,BCHn+BCHe,BCHn+BCHe); % reshape into PC block
       
          
        %l=reshape(l,1,BCHk^2);
        
        %HARD-DECISION on bits
        idx1=(l>0);
        
        % HARD-DECISION SER
        [sym_err_block(ii),ser_block_unc(ii)]=symerr(sym*sqrt(Es),hat_y.');
        
        coded_bits_awgn=zeros(BCHn+BCHe);
        hat_coded_bits=zeros(BCHn+BCHe);
        coded_bits_awgn(idx1)=1;
        
        %% Decoding options
        switch DecoderType
            
        %marked_bits=(abs(l)>R);   % Marks bits based on LLRs
        %hat_coded_bits=PC_BCH_MarkedSR_Decoder(coded_bits_awgn,BCHnu,BCHt,BCHe,I,l,R,w);    % Decode
            case 'iBDD-SR'
        hat_coded_bits=PC_BCH_SR_Decoder(coded_bits_awgn,BCHnu,BCHt,BCHe,I,l,w);        % Decode
            case 'SABM'
        hat_coded_bits=PC_BCH_Marked_Decoder_New(coded_bits_awgn,BCHnu,BCHt,BCHe,I,l,R,ItThresh); % Decode
        %hat_coded_bits=PC_BCH_Marked_Decoder2(coded_bits_awgn,BCHnu,BCHt,BCHe,I,l,R); % Decode
            case 'iBDD'    
        hat_coded_bits=PC_BCH_Decoder(coded_bits_awgn,BCHnu,BCHt,BCHe,I);          % Decode        
            case 'SABM-SR'
        %hat_coded_bits=PC_BCH_MarkedSR_Decoder(coded_bits_awgn,BCHnu,BCHt,BCHe,I,l,R,w); % Decode
        hat_coded_bits=PC_BCH_MarkedSR_Decoder(coded_bits_awgn,BCHnu,BCHt,BCHe,I,l,R,w,ItThresh); % Decode
            %otherwise
             case 'SABM-SR-ST'
        hat_coded_bits=PC_BCH_MarkedSR_Decoder2(coded_bits_awgn,BCHnu,BCHt,BCHe,I,l,R,w,ItThresh); % Decode
        %hat_coded_bits=PC_BCH_Decoder(coded_bits_awgn,BCHnu,BCHt,BCHe,I,l,R,w);    % Decode
        end
        
        
        %% FER/BER/SER calculation
        %err_block_unc(ii)=sum(sum(coded_bits_awgn(1:BCHk,1:BCHk)~=info_bits));
        err_block_unc(ii)=sum(sum(coded_bits_awgn~=coded_bits));             % number of errors in 1 PC block (pre-FEC)
        
        err_block(ii)=sum(sum(hat_coded_bits(1:BCHk,1:BCHk)~=info_bits));    % number of errors in 1 PC block

        ber_block_unc(ii)=err_block_unc(ii)/(BCHk^2);                       % BER per PC block (pre-FEC)

        ber_block(ii)=err_block(ii)/(BCHk^2);                               % BER per PC block
        
        if err_block(ii)>0
        FrameErr=FrameErr+1;
%         if isequal(mod(FrameErr,FrameErrStop/5),0) && FrameErr~=0  % Visualise frame error accumulation
%            FrameErr)
%         end
        end
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
    FrameErrTot(ss)=FrameErr+FrameErrTot(ss);
    disp([num2str(FrameErrTot(ss)) ' Frame Errors']);
    k=k+1;
    end
   
    BERunc(ss)=UncBitErrTot(ss)/(length(UncBitErr{ss})*Nb); 
    SERunc(ss)=UncBitErrTot(ss)/(length(UncBitErr{ss})*Ns);
    BER(ss)=BitErrTot(ss)/(length(BitErr{ss})*(BCHk^2)); 
    fprintf('SNR=%4.2f dB, Uncoded BER=%e \n',snrdB(ss),BERunc(ss));
    fprintf('SNR=%4.2f dB, Uncoded SER=%e \n',snrdB(ss),SERunc(ss));
    fprintf('SNR=%4.2f dB, BER=%e \n',snrdB(ss),BER(ss));
   
end

% Plot results
%figure(1);
semilogy(snrdB,SERunc,'go-'); hold on;
semilogy(snrdB,BERunc,'b*-');
semilogy(snrdB,BER,'sr-');
grid on;hold on;
axis([19,21,1e-9,1e-1]);

if Save 
    if isunix
s_path='/home/gliga/Dropbox/Work/TUe/MATLAB/MATLAB Sims Data/EnhancedHDdecoding/MarkedSR_PCdecoding/NewResults/SABM-SR/';
%s_path='/home/gabriele/Dropbox/Work/TUe/MATLAB/MATLAB Sims Data/EnhancedHDdecoding/MarkedSR_PCdecoding/OECC_Results/';    
    elseif ispc
%s_path='E:\Users\Gabriele\Dropbox\Work\TUe\MATLAB\MATLAB Sims Data\EnhancedHDdecoding\MarkedSR_PCdecoding\NewResults\SABM\';
s_path='E:\Users\Gabriele\Dropbox\Work\TUe\MATLAB\MATLAB Sims Data\EnhancedHDdecoding\MarkedSR_PCdecoding\OECC_Results\';    
%s_path='E:\User Files\Gabriele\Dropbox\Work\TUe\MATLAB\MATLAB Sims Data\EnhancedHDdecoding\MarkedSR_PCdecoding\NewResults\';
    else  
     error('Unexpected operative system');   
    end
    
if ~exist(s_path,'dir')
mkdir(s_path); 
end
disp('Saving...');
switch DecoderType
    case 'iBDD'
    filename=['iBDD_PAM' num2str(M) '_PC_BCH_n=' num2str(BCHn) '_t=' num2str(BCHt) '_e=' num2str(BCHe) '_NoIt=' num2str(I) '.mat'] ;
    save([s_path filename],'snrdB','BERout_unc','SERout_unc','BERout','BERunc','SERunc','BER','r');
    case 'iBDD-SR'
    filename=['iBDD-SR_PAM' num2str(M) '_PC_BCH_n=' num2str(BCHn) '_t=' num2str(BCHt) '_e=' num2str(BCHe) '_NoIt=' num2str(I) '.mat'] ;
    save([s_path filename],'snrdB','BERout_unc','SERout_unc','BERout','BERunc','SERunc','BER','r','w');
    case 'SABM'
    filename=['iBDD-SABM_PAM' num2str(M) '_PC_BCH_n=' num2str(BCHn) '_t=' num2str(BCHt) '_e=' num2str(BCHe) '_NoIt=' num2str(I) '_R=' num2str(R) '_SABM_NoIt=' num2str(ItThresh) '_InterLeaved_MaxLog.mat'] ;
    save([s_path filename],'snrdB','BERout_unc','SERout_unc','BERout','BERunc','SERunc','BER','r','R');
    case 'SABM-SR'
    filename=['iBDD-SABM-SR_PAM' num2str(M) '_PC_BCH_n=' num2str(BCHn) '_t=' num2str(BCHt) '_e=' num2str(BCHe) '_NoIt=' num2str(I) '_MarkingIter= ' num2str(ItThresh) '_R='  num2str(R) '.mat'] ;
    save([s_path filename],'snrdB','BERout_unc','SERout_unc','BERout','BERunc','SERunc','BER','r','R','w');
    case 'SABM-SR-ST'
    filename=['iBDD-SABM-SR-ST_PAM' num2str(M) '_PC_BCH_n=' num2str(BCHn) '_t=' num2str(BCHt) '_e=' num2str(BCHe) '_NoIt=' num2str(I) '_MarkingIter= ' num2str(ItThresh) '_R='  num2str(R) '.mat'] ;
    save([s_path filename],'snrdB','BERout_unc','SERout_unc','BERout','BERunc','SERunc','BER','r','R','w');
end

end

