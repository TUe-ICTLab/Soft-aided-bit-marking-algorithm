function [PosFECBER,GMIav,HDMIav]=simulation_core(m,EsN0dB,BCH,R,flag_display)
% Assumes a channel interleaver generated randomly every block
% MI and GMI are calcualted via Monte Carlo for all constellations (improve GMI results for dense const.)
% The last argument (flag_display) is optional and when set to 1, only one
% block will be simulated (and no files will be saved). This is handy to
% find the SNR range to simulate.
%
% Alex Alvarado
% June 2016

if nargin==4, flag_display=0;end
rng('shuffle');
R_str=sprintf('%1.2g',R);R_str(find(R_str=='.'))='_';
% Modulation
M=2^m;  % Number of constellation symbols. Up to 4096 (for 3 bits/dim)
N=1;    % Number of dimensions
[X,N,L] = get_constellations(M);
[L,idx]=sort(L); % Synchronoulsy permute X and L so that L is the NBC (needed for the mapper)
X=X(idx,:);
Lbin = dec2bin(L)-48;
% Find Subconstellations defined by the labeling
Ik0 = zeros(M/2,m);Ik1 = zeros(M/2,m);
for kk=1:m
    pntr0=1;
    pntr1=1;
    for i=1:M
        if Lbin(i,kk)==0
            Ik0(pntr0,kk)=i;
            pntr0=pntr0+1;
        else
            Ik1(pntr1,kk)=i;
            pntr1=pntr1+1;
        end
    end
end
%SNR definitions
EsN0    = 10^(EsN0dB/10);
N0      = 1/EsN0;
EbN0    = EsN0/(R*m);
EbN0dB  = 10*log10(EbN0);
% Vector definitions
SCC.bl=(BCH.n+1)/2;                 % Length of the SCC block
NBlocks=50;                        % Number of SCC blocks to transmit
if (R>0.8 & flag_display);NBlocks=10;end
if (R>0.9 & flag_display);NBlocks=2*BCH.w;end
n=NBlocks*SCC.bl^2;                 % Number of coded bits to transmit
k=NBlocks*SCC.bl*(BCH.k-SCC.bl);    % Number of info bits to transmit
fprintf('We will transmit %i block and %i coded bits\n',NBlocks,n);
Ns          = n/m;          % n an integer multiple of bits/sym
NsToRemove  = 0; % By default we won't remove any symbols
if mod(n/m,1)~=0 % But we will, if not an integer number of bits/sym
    Ns = ceil(n/m);
    NsToRemove = Ns-floor(n/m);
end
infbits     = zeros(k,1);
kwoLast     = length(infbits(1:end-SCC.bl^2)); % Info bits without last block
codedbits   = zeros(n,1);
intcodedbits= zeros(n+(Ns*m-n),1); % Here we take into account the possibility that we need to add padding bits
int_seq     = zeros(1,n);
x           = zeros(Ns,N);
y           = zeros(Ns,N);
z           = zeros(Ns,N);
l           = zeros(m,Ns);
Lpunct      = zeros(n,1);
intLpunct   = zeros(n,1);
hat_codedbits= zeros(n+(Ns*m-n),1); % Here we take into account the possibility that we need to add padding bits
hat_codedbits_preSCC= zeros(n,1);
hat_infbits	 = zeros(1,k);
% Simulation parameters
Minberc     =5000;            % Count at least 5k bits in error
Maxblocks 	= 5000;           	% Sent maximum 5,000 blocks for each SNR point (should take a few hours)
Minblocks   = 200;              % Sent at least 200 blocks for each SNR point (should take a few hours)
if R>0.85
    Minberc=Minberc/10;
    Maxblocks=Maxblocks/10;
    Minblocks=Minblocks/10;
end
flag_end    = 1;
PreFECBERC	= zeros(1,Maxblocks);
PosFECBERC	= zeros(1,Maxblocks);
MI          = zeros(1,Maxblocks);
GMI         = zeros(1,Maxblocks);
BERC        = 0;
FERC        = 0;
sentblocks  = 0;
pref_save_dir   = strcat('results/',R_str,'/',num2str(2^m),'PAM/');
if ispc,pref_save_dir(pref_save_dir=='/')='\';end
dd=isdir(pref_save_dir);if ~dd,mkdir(pref_save_dir);end
res_name=[pref_save_dir,num2str(2^m),'PAM_',num2str(EsN0dB,4),'_dB_m_',num2str(BCH.m),'_t_',num2str(BCH.t),'_w_',num2str(BCH.w),'_iter_',num2str(BCH.iter)];

% If a normal simulation was saved, don't do it again.
fid=fopen([res_name,'.mat']);
if fid ~= -1
    f=load([res_name,'.mat']);
    if exist([res_name,'.mat'])==2 & ~isfield(f,'merged_flag') 
        fprintf('File already exist. Leaving simulation... \n');
        return
    end
end

% Main loop
ii=0;
while flag_end
    %% Display
    ii=ii+1;if mod(sentblocks,10)==0,fprintf('.');end
    %% Info bits and coding
    infbits     = (randi([0 1],k,1));
    codedbits=StaircaseEncoding(infbits,BCH.m,BCH.t).';
    %% Randomly generate interleaver
    int_seq     = randperm(n);
    intcodedbits(1:n,1)= codedbits(int_seq); % If there are padding bits, they are not touched.
    if m==1
        x   	= 2*intcodedbits-1;
    else
        x  = mapper(X,intcodedbits,m,Ns,N);
    end
    %% AWGN Channel
    z           = sqrt(N0/2)*randn(Ns,N);               % Es=1 for all constellations and noise power is N0/2.
    y           = x + z;
    %% Demodulation
    if m==1
        l	= 4*EsN0*y.';
    else
        l   =  demapper(X,Lbin,EsN0dB,y,Ik1,Ik0,1);  %%% MEX
    end  
    %% Hard-decision based on L-values
    hat_codedbits       = 1*(l(:)>0);
    PreFECBERC(1,ii)    = sum(intcodedbits(1:n)~=hat_codedbits(1:n)); % avoid padding bits
    %% MI calculation via Monte Carlo integration
    MI(1,ii)            = MI_uniform_MonteCarlo(x,y,N0); %%% MEX
    %% HD-MI calculation
    cop=PreFECBERC(1,ii)/n; % Crossover probability
    HDMI(1,ii)          = m*(1+(1-cop)*log2(1-cop)+cop*log2(cop)); % HD-MI
    %% GMI calculation via Monte Carlo integration
    for kk=1:m
        GMI(1,ii)= (1-1/Ns*sum(log2(1+exp((-1).^(intcodedbits(kk:m:end).').*l(kk,:))))) + GMI(1,ii);
    end
    % Deinterleaver
    hat_codedbits_preSCC(int_seq)= hat_codedbits(1:n);
    % Decoding
    %[hat_infbits,~]=StaircaseDecoding(hat_codedbits_preSCC,BCH.m,BCH.t,BCH.w,BCH.iter);
    [hat_infbits,~]=StaircaseDecodingGenie(hat_codedbits_preSCC,codedbits,BCH.m,BCH.t,R,BCH.w,BCH.iter);
    % BER without taking into account the last block
    errC=sum(hat_infbits(1:end-SCC.bl^2)~=infbits(1:end-SCC.bl^2).');
    % Error counting
    PosFECBERC(1,ii)= errC+PosFECBERC(1,ii);
    FERC        = ~isequal(infbits(1:end-SCC.bl^2),hat_infbits(1:end-SCC.bl^2))+FERC;
    sentblocks  = sentblocks+1;
    if flag_display
        fprintf('-------------------------------------------------\n');
        %fprintf('Eb/N0 = %f dB \n',EbN0dB);
        fprintf('Current pre FEC BER is %2.2e. ',PreFECBERC(ii)/n);
        fprintf('Current post FEC BER is %2.2e \n',PosFECBERC(ii)/kwoLast);
        %fprintf('Current (normalized) GMI is %2.5f \n',GMI(ii)/m);
        %fprintf('Current (normalized) HDMI is %2.5f \n',HDMI(ii)/m);
        %fprintf('Current (normalized) MI is %2.5f \n',MI(ii)/m);
        GMIav=GMI(ii);
        HDMIav=HDMI(ii);
        PosFECBER=PosFECBERC(ii)/kwoLast;
        return
    end
    
    if sum(PosFECBERC) >= Minberc && sentblocks >= Minblocks
        flag_end = 0;
    end
    if sentblocks >= Maxblocks   
        fprintf('We did not count enough errors... Saving temp file...');
        pref_save_dir_in   = strcat('results/',R_str,'/',num2str(2^m),'PAM/incomplete/');
        if ispc,pref_save_dir_in(pref_save_dir_in=='/')='\';end
        dd=isdir(pref_save_dir_in);
        if ~dd,mkdir(pref_save_dir_in);end
        res_name=[pref_save_dir_in,num2str(2^m),'PAM_',num2str(EsN0dB,4),'_dB_m_',num2str(BCH.m),'_t_',num2str(BCH.t),'_w_',num2str(BCH.w),'_iter_',num2str(BCH.iter),num2str(randi(10000,1))];
        flag_end = 0;
    end
end
res_name = strcat(res_name,'.mat');
% Variablas to save
PreFECBERC=PreFECBERC(1:ii);
PosFECBERC=PosFECBERC(1:ii);
MI=MI(1:ii);
GMI=GMI(1:ii);
% Final calculations
total_coded_bits= sentblocks*n;                     % Number of sent coded bits
PreFECBER       = PreFECBERC/n;                     % Pre FEC BER
PreFECBERav     = sum(PreFECBERC)/total_coded_bits; % Average pre FEC BER
total_info_bits = sentblocks*kwoLast;                     % Number of sent info bits
PosFECBER       = PosFECBERC/kwoLast;                     % Post FEC BER
PosFECBERav     = sum(PosFECBERC)/total_info_bits;  % Average Post FEC BER
FER             = FERC/sentblocks;                  % FER
MIav            = mean(MI);                 	% Average MI (only over nonzero values)
HDMIav          = mean(HDMI);                 	% Average HDMI (only over nonzero values)
GMIav           = mean(GMI);                	% Average GMI (only over nonzero values)
% Display results
fprintf('\n-------------------------------------------------\n');
fprintf('Inf. Bit errors = %d \n',sum(PosFECBERC));
fprintf('Coded Bit errors = %d \n',sum(PreFECBERC));
fprintf('Frame errors = %d \n',FERC);
fprintf('MI = %d \n',MIav);
fprintf('GMI = %d \n',GMIav);
fprintf('Pre FEC BER = %d \n',PreFECBERav);
fprintf('Post FEC BER  = %d\n',PosFECBERav);
fprintf('FER  = %1.6f\n',FER);
fprintf('-------------------------------------------------\n');
% Save
if sentblocks >= Maxblocks % to be merged, save less data
    save(res_name,'X','L','M','m','k','kwoLast','BCH',...
        'PreFECBERav','PosFECBERav','PreFECBERC','PosFECBERC',...
        'FER','FERC','EbN0','total_coded_bits','total_info_bits','sentblocks',...
        'EbN0dB','EsN0','EsN0dB','R','MIav','GMIav','HDMIav');
else
    save(res_name,'X','L','M','m','k','kwoLast','BCH',...
        'PreFECBER','PreFECBERav','PosFECBER','PosFECBERav','PreFECBERC','PosFECBERC',...
        'FER','FERC','EbN0','total_coded_bits','total_info_bits','sentblocks',...
        'EbN0dB','EsN0','EsN0dB','R','MI','MIav','GMI','GMIav','HDMIav');
end
return
