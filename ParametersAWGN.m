%% Parameters

PlotConst=0;

%% Encoder/Decoder parameters
DecoderType='iBDD';   % Options are: 'iBDD', 'iBDD-SR', 'SABM', 'SABM-SR', 'SABM-SR-ST'  
BCHnu=7;
BCHt=2;
BCHe=1;
BCHn=(2^BCHnu)-1;                   % Coded bits, without the extended bits
BCHk=(2^BCHnu)-1-BCHnu*BCHt;        % Information bits
I=10;                               % iBDD decoder iterations
Nbb=(BCHn+BCHe)^2;                   % Number of bits per PC block
r=BCHk^2/(BCHn+BCHe)^2;

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


Inter=1;
flag_demapper=0;
FrameErrStop=1e1;
Nblocks_min=10*FrameErrStop;
Nblocks_max=1e5;    % x Nblocks_exp 
Nllr=100;
%L=de2bi(0:M-1,m);


% Modulation parameters 
M=2;                   % Constellation cardinality/per dimension 
m=log2(M);             % Bits/symbol
N=1;                   % Constellation dimensions
Ns=ceil(Nbb/m);        % Number of symbols in a block


% Sweeping parameters 
snrdB=5.5:.1:6;

    