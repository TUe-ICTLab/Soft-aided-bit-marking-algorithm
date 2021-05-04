%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                 staircase code- by Yi Lei                   %%%
% Data is the input information data, it should be a row vector   %
%  m is the field extension degree for Galois field               %
%  t is the error correct ability                                 %
%  This is for extend BCH with 1 extra parity bit                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  StaircaseCodeWord=StaircaseEncoding(data,m,t)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                 Input parameters                   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Staircasem=m; %BCH code m
Staircaset=t; %BCH code t
InformationData = data; 
InformationDataLength=length(InformationData);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                   Initialization                   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
BCHn=(2^Staircasem)-1;                   % The codeword length for BCH encode
BCHk=(2^Staircasem)-1-Staircasem*Staircaset;  % The number of inforamtion bit in each BCH codeword 
%  s=round((2^Staircasem)-1-((2*Staircasem*Staircaset)/(1-StaircaseR))); % The number of information bits to shorten each BCH codeword
s=0;
BlockLength=(BCHn+1-s)/2;       % The length of each block
SuqareOfBlockLength=BlockLength*BlockLength;
DataLengthInBlock=BCHk-s-BlockLength; % The information data length in each block(except B0)
BlockLengthXDataLengthInBlock=BlockLength*DataLengthInBlock;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%    Divide these data into a sequence of blocks     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NumofBlocks = ceil(InformationDataLength/(BlockLengthXDataLengthInBlock))+1; % Compute the needed block numbers,including B0.
c{1,1}=zeros(BlockLength);   % Initilize B0 to all zeros. c is used to store the blocks
%divide the information data into a sequence of blocks
for i=1:NumofBlocks-1
        Block_temp(1:BlockLengthXDataLengthInBlock)=InformationData((1+(i-1)*BlockLengthXDataLengthInBlock):(BlockLengthXDataLengthInBlock+(i-1)*BlockLengthXDataLengthInBlock));   
        Block=reshape(Block_temp,[DataLengthInBlock,BlockLength]);
        Block=Block';
        c{1,i+1}=Block;  % c store these blocks, c{1,1} store B0 that is all zeros. 
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                Using BCH encoding                  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
 for i=1:NumofBlocks-1 
        CombineTwoBlocks=[c{1,i}',c{1,i+1}];
        CombineTwoBlocks=[zeros(BlockLength,s),CombineTwoBlocks];
        Paritybits=zeros(BlockLength,(BlockLength-DataLengthInBlock-1));  % Storing the caculated parity bits of each block
        for j=1:BlockLength
            BCHCodeword = bchenc_mex(Staircasem, BCHn,Staircaset,CombineTwoBlocks(j,1:BCHk)); %encoding for each row,return codeword of each row
            Paritybits(j,:)=BCHCodeword((BCHk+1):BCHn); % Storing the caculated paritybits 
        end    
        CombineTwoBlocks=[CombineTwoBlocks, Paritybits];
        c{1,i+1}=[c{1,i+1},Paritybits]; 
        c{1,i+1}=[c{1,i+1},mod(sum(CombineTwoBlocks,2),2)];
 end

%Using to store the encodedword row-by-row, then transmitted row-by-row
EncodedInformationData=[];
for i=2:NumofBlocks
     temp= (c{1,i})';
     EncodedInformationData=[EncodedInformationData (temp(:))'];
end 
    
StaircaseCodeWord=EncodedInformationData;