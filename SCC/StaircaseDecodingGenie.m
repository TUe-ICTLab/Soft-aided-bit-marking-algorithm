%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                Staircase decoding algorithm                        %%%
%%%                           --by Yi Lei                              %%%
%  Received_data is the received data, it should be a row vector         %
%  m is the field extension degree for Galois field                      %
%  t is the error correct ability                                        %
%  Window_Length is the decoding window size                             %
%  Iteration_Times is the iteration times for each decoding window       % 
%  Decode_Information is the decoded information bits and is a row vector%
%  Full_DecodeWord is the decoded codeword, including the information    %
%  bits and parity bits . It is a row vector.                            %
%  This is for extended BCH with 1 extra parity bit                      %                      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Decode_Information,Full_DecodeWord]=StaircaseDecodingGenie(Received_Codeword,Transmitted_data,Staircasemm,Staircasett,StaircaseRR,Window_Length,Iteration_Times)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                 Decoding parameters                %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
BCH_n=(2^Staircasemm)-1;                   % The length of  BCH codeword
BCH_k=(2^Staircasemm)-1-Staircasemm*Staircasett;      % The number of inforamtion bit in each codeword 
s=round(2^Staircasemm-((2*Staircasemm*Staircasett+2)/(1-StaircaseRR))); % the number of shorten bit
Block_Length=(BCH_n-s+1)/2;   % The length of each block
Square_of_Block_Length=Block_Length*Block_Length;    
DataLength_In_Block=BCH_k-Block_Length-s; % The information data length in each block(except B0)
Received_Data_Length=length(Received_Codeword);% the length of received data length 
Num_of_Blocks=ceil(Received_Data_Length/Square_of_Block_Length); % the received Block number

Temp_row=zeros(Block_Length,BCH_n+1);
Err_poly_row=zeros(1,BCH_n+1);

%reshape the received data into blocks
ReceivedData_in_Blocks=cell(1,Num_of_Blocks+1);
TransmitData_in_Blocks=cell(1,Num_of_Blocks+1);
for i=1:Num_of_Blocks  % read Slidding_Window_Length blocks
     Temp=(i-1)*Square_of_Block_Length;
     Block_Temp=Received_Codeword(Temp+1:Temp+Square_of_Block_Length);
     ReceivedData_in_Blocks{1,i+1}= (reshape(Block_Temp,[Block_Length,Block_Length]))';
     Block_Temp=Transmitted_data(Temp+1:Temp+Square_of_Block_Length);
     TransmitData_in_Blocks{1,i+1}= (reshape(Block_Temp,[Block_Length,Block_Length]))';
end
ReceivedData_in_Blocks{1,1}=zeros(Block_Length);   % Initilize B0 to all zeros.
TransmitData_in_Blocks{1,1}=zeros(Block_Length);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%            Starting staircase decoding             %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Blocks_in_Window=cell(1,Window_Length);
Data_in_Window=cell(1,Window_Length);

Num_of_decoding_window=Num_of_Blocks+1-(Window_Length-1);


for i=1:Num_of_decoding_window
    % select Window_Length Blocks as a decoding group
    for j=1:Window_Length 
        Blocks_in_Window{1,j}=ReceivedData_in_Blocks{1,(i-1)+j};
        Data_in_Window{1,j}=TransmitData_in_Blocks{1,(i-1)+j};
    end
    for ii=1:Iteration_Times    
        for jj=Window_Length:-1:2
            Temp_array_row=[zeros(Block_Length,s),Blocks_in_Window{1,jj-1}',Blocks_in_Window{1,jj}];
            Temp_array_data=[zeros(Block_Length,s),Data_in_Window{1,jj-1}',Data_in_Window{1,jj}]; 
            Decoding_Result=Temp_array_row;
            for kk=1:Block_Length
                 temp=xor(Temp_array_row(kk,:),Temp_array_data(kk,1:BCH_n+1));
                 if sum(temp)<= Staircasett % if the error bits <= t, correct it
                        Decoding_Result(kk,1:BCH_n+1)=Temp_array_data(kk,1:BCH_n+1);
                 else
                     if sum(temp) > (Staircasett+1)  % if the error bits > (t+1), then decode it, accept the decoding result of the BCH decoder
                        [Err_poly_row(kk,:), Temp_row(kk,:), status]=eBCHdecoder(Temp_array_row(kk,1:BCH_n+1),Staircasemm,BCH_k,Staircasett);
                        if (status==1) && (sum(Err_poly_row(kk,1:s))==0) % decoding successfully
                            Decoding_Result(kk,1:BCH_n+1)=Temp_row(kk,1:BCH_n+1);
                        end
                     end
                 end
                 % if the error bits == t , keep the received codeword
                 % unchange
             end   
            Blocks_in_Window{1,jj-1}=Decoding_Result(:,s+1:s+Block_Length);
            Blocks_in_Window{1,jj-1}=Blocks_in_Window{1,jj-1}';
            Blocks_in_Window{1,jj}= Decoding_Result(:,s+Block_Length+1:BCH_n+1);
        end
    end
    
    for iii=1:Window_Length
        ReceivedData_in_Blocks{1,(i-1)+iii}=Blocks_in_Window{1,iii};
    end  
end

%% recover the information data from the decoded Codeword

%information data in row
Recovered_Data=[];
for i=2:Num_of_Blocks+1
      temp= (ReceivedData_in_Blocks{1,i}(:,1:DataLength_In_Block))';
      Recovered_Data=[Recovered_Data (temp(:))'];
end
Decode_Information=Recovered_Data;  

%codeword in row
Recovered_Codeword=[];
for i=2:Num_of_Blocks+1
      temp= (ReceivedData_in_Blocks{1,i}(:,1:Block_Length))';
      Recovered_Codeword=[Recovered_Codeword (temp(:))'];
end
Full_DecodeWord=Recovered_Codeword;


        



