function [time,error_vec]=Product_IMP_code_BCH_AWGN_normnoise_eBCH_misfree(iter,SNR,v,t,s,modord)

% Intrinsic decoding
% Note: Iteration is bi-directional
% Extended BCH
% Alireza Sheikh,  asheikh@chalmers.se

clc
tic

%% 
%%%%%%% ======================== Parity check (Reed Solomon) ========================= %%%%%%%%
EbN0=SNR;
tt=t;
rng('shuffle')
rng(300)

block_len_original=2^(v)-1;       % original block length
inf_len_original=2^(v)-v*tt-1;    % information length
block_len=2^(v)-1-s+1;            % shortened original block length + extension
inf_len=2^(v)-v*tt-1-s;           % shortened information length 
kk=inf_len;
nn=block_len;

%%
%%%%% ================ scaling factor ==========================%%%%%

ASK_ord=2^modord;
kk=ASK_ord/2;
Popt=(1/ASK_ord)*(ones(1,ASK_ord));
vec=[-(2*kk-1):2:(2*kk-1)];
PP=10^(SNR/10);
delopt=sqrt(PP/(sum(Popt.*(vec.^2))));

Rate_c=(inf_len/block_len)^2;
Es_N0=((10^(EbN0/10))*modord)*Rate_c; 
sigma=sqrt(1/(2*Es_N0));

%% 
%%%%%%%%%%% ======================== Transmitt and receive loop ========================= %%%%%%%%

error_vec=zeros(1,length(SNR));
error_vec1=zeros(1,length(SNR));

for kk=1:length(SNR)

Eb=SNR(kk) 

%Rate_c=(inf_len/block_len)^2;  
%P=((10^(Eb/10))*modord)*Rate_c;   
        
%     if 2^modord==16
%         
%     norm=10;    
%            
%     elseif   2^modord==32 
%         
%     norm=20;    
%         
%     elseif   2^modord==4         
%         
%     norm=2;    
%     
%     elseif   2^modord==64 
%         
%     norm=42;    
%         
%     elseif   2^modord==128
%         
%     norm=82;    
%         
%     elseif   2^modord==256
%         
%     norm=170;    
%         
%     elseif   2^modord==512
%         
%     norm=330;    
%         
%         
%     end
% Symbol Length 

gamma_bar_dB=Eb;
    
  if (gamma_bar_dB<7)
     
    mo=10000;

 elseif (gamma_bar_dB<17.7 && gamma_bar_dB>=2)
     
    mo=10000;    

 elseif (gamma_bar_dB>17.7)
     
    mo=10000;       
        
 end          
     
error=0;
error1=0;  
count=1;

for indt=1:(mo)   
%rng('shuffle') 
%indt    
    
%%
%%%%%%%% ================= Encoder & Modulation Initialization ==================== %%%%%%%%

enc = comm.BCHEncoder(block_len_original,inf_len_original);
%dec = comm.BCHDecoder(block_len_original,inf_len_original,'NumCorrectedErrorsOutputPort',true, 'GeneratorPolynomialSource', 'Auto');

%%%%%%%%%%% ============================= Encoder ================================== %%%%%%%%

inf_stair_block=zeros(inf_len,inf_len);            % Inf. staircase block
encod_stair_block1=zeros(inf_len,block_len);       % Block staircase
encod_stair_block=zeros(block_len,block_len);      % Block staircase

%%%%%%%%%% ============================= Inf block =========================== %%%%%%%%%%%%%%

inf_stair_block = randi([0 1],inf_len,inf_len);     % Data symbols

%%%%%%%%% ========================== Encode Inf block ============================ %%%%%%%%%%
%timm=0;
%%%% Row encoding  
for jj=1:inf_len
en_temp1=inf_stair_block(jj,:)';
%tic
en_temp2=step(enc,[zeros(s,1);en_temp1])';
%ttime=toc;
%timm=timm+ttime;
encod_stair_block1(jj,:)=[en_temp2(s+1:end) mod(sum(en_temp2),2)];
end
%timm
%%%% Column encoding
for jj=1:block_len    
en_temp1=encod_stair_block1(:,jj);
en_temp2=step(enc,[zeros(s,1);en_temp1]);
encod_stair_block(:,jj)=[en_temp2(s+1:end);mod(sum(en_temp2),2)];
end

% Interleaving
% state=rand(1);
% encod_stair_block1=randintrlv(encod_stair_block,state);
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

y_demod=reshape(decode_out4,[block_len,block_len]);

% Deinterleaving
y_demod=randdeintrlv(y_demod,state);
y_decode_gf_inf=y_decode_func(encod_stair_block,y_demod,block_len,inf_len,iter,v,s);

% length(find(y_demod~=encod_stair_block))           % error before decoding (block error)
% length(find(y_decode_gf_inf~=encod_stair_block))   % error after decoding (block error)

%%%%%%%%%% ================= Remove transient ============================ %%%%%%%%%%%%

y_decode_gf_inft=y_decode_gf_inf(1:inf_len,1:inf_len);
inf_stair_blockt=encod_stair_block(1:inf_len,1:inf_len);

%%%%%%% ====================== Inf. symbols to Inf. bits ======================= %%%%%%%

% vec_sym=y_decode_gf_inft(:);
% vec_sym1=inf_stair_blockt(:);

% l_l=length(vec_sym);
% tempt=reshape(vec_sym,v,l_l/v);
% tempt1=reshape(vec_sym1,v,l_l/v);
% senddata_tt=bi2de(tempt');
% senddata_tt1=bi2de(tempt1');

%% 
%%%%%%%% =============================== Error calc =============================== %%%%%%%%
%%%%%%%% ============================ Inf. bit error rate ========================== %%%%%%%
 
vvabit=find(y_decode_gf_inft~=inf_stair_blockt);
error=error+length(vvabit);
count=count+1;
[az,bz]=size(y_decode_gf_inft);
%error/(count*az*bz)

if error>150000
    
    break;
    
end
 
 %% Inf. Symbol error rate
 
% vvabit=find(senddata_tt~=senddata_tt1);
% error1=error1+length(vvabit);


end

% [az,bz]=size(senddata_tt);
% error_vec1(kk)=error1/(mo*az*bz);

[az,bz]=size(y_decode_gf_inft);
error_vec(kk)=error/(count*az*bz);

end
%% 
%%%%%%%%%%% ============================= Plot =============================== %%%%%%%%%%%%
% hold on
% semilogy(Eb,error_vec,'*k');
% xlabel('E_b/N_0(dB)')
% ylabel('BER')
time=toc;
end

%%
%%%%%%%%%%% ============================== Decode staircase code ============================= %%%%%%%%%%%
function decodeout1=y_decode_func(vecginie,vec,block_comp,inf_comp,itermax,v,s)

old_vec=vec;
new_vec=vec;
vecginie1=vecginie.';

%block_comp=block_comp-1;

t_code=(block_comp-inf_comp-1)/v;
enc = comm.BCHEncoder(block_comp+s-1,inf_comp+s);
%dec = comm.BCHDecoder(block_comp+s,inf_comp+s,'NumCorrectedErrorsOutputPort',true);
%dec = comm.BCHDecoder(block_comp+s,inf_comp+s);
dec = comm.BCHDecoder(block_comp+s-1,inf_comp+s,'NumCorrectedErrorsOutputPort',true, 'GeneratorPolynomialSource', 'Auto');

for iter=1:itermax


for i=1:block_comp
    
out1=extendedBCHdecoder(vecginie(i,:),enc,dec,s,t_code,vec(i,:));    
vec(i,:)=out1;   

end

vec2=vec';

for i=1:block_comp
    
out2=extendedBCHdecoder(vecginie1(i,:),enc,dec,s,t_code,vec2(i,:));
vec2(i,:)=out2;

end

vec=vec2';

new_vec=vec;

if length(find(new_vec~=old_vec))==0
    
    break;
    
else
    
old_vec=vec;

end

end

decodeout1=vec;

end
%%
function output1=extendedBCHdecoder(vecginiee,enc,dec,s,t,input1)

output1=zeros(1,length(input1));

[cv,error_idicate]=step(dec,[zeros(s,1);input1(1,1:end-1)']);

if error_idicate==-1
    
 output1=input1;
 
else
     
 de1=sum(input1)+error_idicate;
 de=mod(de1,2);
 
 if (de+error_idicate)>t
  
 output1=input1;   
        
 else
 
 cvp=step(enc,cv);
 output1(1:end-1)=cvp(s+1:end)';   
 output1(end)=mod(input1(end)+de,2);    
 
 if (vecginiee~=output1)                       % miscorrection free (Ginei operation)
 output1=input1; 
 end    
 
 end
    
 
end

end

function x2=symbol_binlabels(x_in,pamord)

x1=bi2de(x_in);
x2=gray2bin(x1,'pam',pamord);

end