%% Script simulating the performance of a product code over an AWGN channel using different decoding algorithms
% The script implements the following FEC decoder options 
%  1: iterative bounded distance decoding (iBDD)
%  2: iBBD-scaled reliability (SR) [1]
%  3: Soft-aided Bit Marking (SABM) algorithm    [2]
%  4: SABM-SR [3]
%  5: SABM-SR-scaled threshold (ST) which is a variant of SABM-SR

%% References 
% [1] A. Sheikh, A. Graell i Amat and G. Liva, "Binary Message Passing Decoding of Product-Like Codes," in IEEE Transactions on Communications, vol. 67, no. 12, pp. 8167-8178, Dec. 2019, doi: 10.1109/TCOMM.2019.2940180.
% [2] Lei, Y., Chen, B., Liga, G., Deng, X., Cao, Z., Li, J., Xu, K., & Alvarado, A. (2019). Improved decoding of staircase codes: The soft-aided bit-marking (SABM) algorithm. In IEEE Transactions on Communications (Vol. 67, Issue 12, pp. 8220–8232). IEEE.   
% [3] G. Liga, A. Sheikh and A. Alvarado, "A novel soft-aided bit-marking decoder for product codes," 45th European Conference on Optical Communication (ECOC 2019), 2019, pp. 1-4, doi: 10.1049/cp.2019.0804.


startup;
clear all;
ParametersAWGN;


%% Constellation generation
X=-(M-1):2:M-1;      % PAM constellation
Es=mean(abs(X).^2);
X=X/sqrt(Es);   
%Mode='soft';

flag_demapper=1;       % Max-log algorithm for LLR calculation
L=de2bi(0:M-1,m);


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
    LLR=zeros(BCHn+BCHe,BCHn+BCHe,Nllr);
    
    while FrameErrTot(ss)<FrameErrStop && (k+1)*Nblocks_min<Nblocks_max
      sym_err_block=zeros(1,Nblocks_min);
      err_block_unc=zeros(1,Nblocks_min);
      err_block=zeros(1,Nblocks_min);
      ber_block_unc=zeros(1,Nblocks_min);
      ser_block_unc=zeros(1,Nblocks_min);
      ber_block=zeros(1,Nblocks_min);
      disp(['Decoding Blocks ' num2str(k*Nblocks_min) ' to ' num2str((k+1)*Nblocks_min) '...']);
      FrameErr=0;         % Frame errors per chunk of code blocks decoded in parallel 

           
     for ii=1:Nblocks_min
        
        info_bits=randi([0,1],BCHk,BCHk);                               % Generate the info bit matrix
        
        %% PC encoding
        coded_bits=PC_BCH_Encoder(info_bits,BCHnu,BCHt,BCHe);           % Encode the data using a product code    
        coded_bits_res=reshape(coded_bits,1,Nbb);
        
        % Bit interleaver
        map=0;
         if Inter
        [coded_bits_res,map]=RandInter(coded_bits_res);
         end
        
       % Bit-to-symbol mapping        
       pad=mod(length(coded_bits_res),m);
      
       % Padding for number of bits in the block non-mulitple of m
        if ~isequal(pad,0)
         coded_bits_res=[coded_bits_res zeros(1,m-pad)];
        end
        
        
        coded_bits_res=reshape(coded_bits_res,m,length(coded_bits_res)/m);
        sym_dec=bi2de(coded_bits_res.').';
        sym=X(:,sym_dec+1);
        %sym=mapper(X.',coded_bits_res,m,Ns,N).';
        %sym=sym/sqrt(Es);                             % Normalise constellation to unitary energy
        
        
        % N-dimensional AWGN channel
        y=awgn(sym,snr,'measured');                             % Transmit through a AWGN channel
        
         
       %% HARD DEMAPPING
       % Minimum ED symbol decision
       Xext=permute(repmat(X,1,1,length(y)),[1,3,2]);
       yExt=repmat(y,1,1,M);
       ED=sum(abs(yExt-Xext).^2,1);
       %clear RXsymExt Xext;                       % Free up a substantial amount of memory
       [~,y_hat_dec]=min(ED(1,:,:),[],3);
       %y_hat=X(:,y_hat_idx);              % Hard-symbol out 
       hat_y=X(:,y_hat_dec);
       y_hat_dec=y_hat_dec-1;    % Shifts indexes from 0 to M-1 
       
       
         info_bits_awgn=de2bi(y_hat_dec).';
         if ~isequal(pad,0)
         info_bits_awgn=info_bits_awgn(1:end-m+pad);
         end
         info_bits_awgn=reshape(info_bits_awgn,BCHn+BCHe,BCHn+BCHe); 
                 
        %% SOFT symbol-to-bits demapping (LLR calculator)
        l=demapper(X.',L,snr,y.',flag_demapper); % snr is converted to Es/No as per definition in demapper function
        % histogram(l,100);
        
        % Bit de-interleaver
         if ~isequal(pad,0) 
          l=l(1:end-m+pad);
         end
         
         if Inter
          l=RandInter(l,map);    
         end
                   
        l=reshape(l,BCHn+BCHe,BCHn+BCHe); % reshape into PC block
        LLR(:,:,ii)=l;
          
        %l=reshape(l,1,BCHk^2);
        
        %HARD-DECISION on bits
        idx1=(l>0);
        
        % HARD-DECISION SER
        [sym_err_block(ii),ser_block_unc(ii)]=symerr(sym*sqrt(Es),hat_y);
        
        coded_bits_awgn=zeros(BCHn+BCHe);
        hat_coded_bits=zeros(BCHn+BCHe);
        coded_bits_awgn(idx1)=1;
        
        %% Decoding options
        switch DecoderType
            case 'iBDD'    
        hat_coded_bits=PC_BCH_Decoder(coded_bits_awgn,BCHnu,BCHt,BCHe,I);          % Decode             
            case 'iBDD-SR'
        hat_coded_bits=PC_BCH_SR_Decoder(coded_bits_awgn,BCHnu,BCHt,BCHe,I,l,w);        % Decode
            case 'SABM'
        hat_coded_bits=PC_BCH_Marked_Decoder_New(coded_bits_awgn,BCHnu,BCHt,BCHe,I,l,R,ItThresh); % Decode      
            case 'SABM-SR'
        hat_coded_bits=PC_BCH_MarkedSR_Decoder(coded_bits_awgn,BCHnu,BCHt,BCHe,I,l,R,w,ItThresh); % Decode
             case 'SABM-SR-ST'
        hat_coded_bits=PC_BCH_MarkedSR_Decoder2(coded_bits_awgn,BCHnu,BCHt,BCHe,I,l,R,w,ItThresh); % Decode
        end
        
        
        %% FER/BER/SER calculation
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
   
    BERunc(ss)=UncBitErrTot(ss)/(length(UncBitErr{ss})*Nbb); 
    SERunc(ss)=UncBitErrTot(ss)/(length(UncBitErr{ss})*Ns);
    BER(ss)=BitErrTot(ss)/(length(BitErr{ss})*(BCHk^2)); 
    fprintf('SNR=%4.2f dB, Uncoded BER=%e \n',snrdB(ss),BERunc(ss));
    fprintf('SNR=%4.2f dB, Uncoded SER=%e \n',snrdB(ss),SERunc(ss));
    fprintf('SNR=%4.2f dB, BER=%e \n',snrdB(ss),BER(ss));
   
end

% Plot results
%figure(1);
p1=semilogy(snrdB,SERunc,'go-'); hold on;
p2=semilogy(snrdB,BERunc,'b*-');
p3=semilogy(snrdB,BER,'sr-');
grid on;hold on;
axis([min(snrdB),max(snrdB),1e-9,1e-1]);
xlabel('E_s/N_o');
ylabel('BER');
legend([p1,p2,p3],'Pre-FEC SER','Pre-FEC BER','Post-FEC BER');


