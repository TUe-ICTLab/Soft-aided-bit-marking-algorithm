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
%  This is for extended BCH with 1 extra parity bit                      %                                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Decode_Information,Full_DecodeWord]=StaircaseDecoding(Received_data,m,t,Window_Length,Iteration_Times)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                 Decoding parameters                %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
Received_Codeword= Received_data;
Staircasemm=m; %input('m =');   %BCH code m
Staircasett=t; %input('t =');   %BCH code t
Slidding_Window_Length=Window_Length; % decoding of staircase code is performing iteratively over a window of L blocks
Num_of_Iteration= Iteration_Times; % decoding of staircase code is performing iteratively over a window of L blocks

BCH_n=(2^Staircasemm)-1;                   % The length of  BCH codeword
BCH_k=(2^Staircasemm)-1-Staircasemm*Staircasett;      % The number of inforamtion bit in each codeword 
% s=round(2^Staircasemm-1-((2*Staircasemm*Staircasett)/(1-StaircaseRR))); % The number of information bits to shorten each component code
Block_Length=(BCH_n+1)/2;
Suqare_of_Block_Length=Block_Length*Block_Length;
% Block_Length=(BCH_n-s)/2;       % The length of each block
DataLength_In_Block=BCH_k-Block_Length; % The information data length in each block(except B0)
Received_Data_Length=length(Received_Codeword);% the length of received data length 
Num_of_Blocks=ceil(Received_Data_Length/Suqare_of_Block_Length); % caculate the received Block number

%reshape the received data into blocks
for i=1:Num_of_Blocks  % read Slidding_Window_Length blocks
     Temp=(i-1)*Suqare_of_Block_Length;
        Block_Temp=Received_Codeword(Temp+1:Temp+Suqare_of_Block_Length);
        ReceivedData_in_Blocks{1,i+1}= (reshape(Block_Temp,[Block_Length,Block_Length]))';
end
ReceivedData_in_Blocks{1,1}=zeros(Block_Length);   % Initilize B0 to all zeros.

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%            Starting staircase decoding             %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Blocks_in_Window=cell(1,Slidding_Window_Length);
Num_of_decoding_window=Num_of_Blocks+1-(Slidding_Window_Length-1);

for i=1:Num_of_decoding_window
    % select Slidding_Window_Length Blocks as a decoding group
    for j=1:Slidding_Window_Length  % read Slidding_Window_Length blocks
        Blocks_in_Window{1,j}=ReceivedData_in_Blocks{1,(i-1)+j};
    end
    for ii=1:Num_of_Iteration      
        for jj=Slidding_Window_Length:-1:2
            Error_poly=zeros(Block_Length,BCH_n);
            Temp1_array=zeros(Block_Length,BCH_n+1);
            Temp_array=[Blocks_in_Window{1,jj-1}',Blocks_in_Window{1,jj}];  
            for jjj=1:Block_Length
%                         Temp1=bchdec_mex(Staircasemm,BCH_k,Staircasett,Temp_array(jjj,1:BCH_n)); 
%                         Error_poly(jjj,:)=xor(Temp1(1,1:BCH_n),Temp_array(jjj,1:BCH_n));% get the error polynomial
%                         Flag1=Temp1(1,length(Temp1));  % return the status of the decoding results, 0=no errors, 1= decode successfully, others= detecting errors, but decoding failure
%                         if ~isequal(Flag1,1)  % indicate decoding failure
%                             Temp1_array(jjj,1:BCH_n+1)=Temp_array(jjj,1:BCH_n+1);
%                         else
%                             d=sum(Error_poly(jjj,:)); %indicate the distance between the received codeword and the transmitted codeword
%                             de=mod(d+sum(Temp_array(jjj,1:BCH_n+1)),2);
%                             if d+de<=Staircasett
%                                 Temp1_array(jjj,1:BCH_n)=Temp1(1,1:BCH_n);
%                                 Temp1_array(jjj,BCH_n+1)=xor(Temp_array(jjj,BCH_n+1),de);
%                             else
%                                 Temp1_array(jjj,1:BCH_n+1)=Temp_array(jjj,1:BCH_n+1);
%                             end   
%                         end
                [~, Temp1_array(jjj,1:BCH_n+1), ~]=eBCHdecoder(Temp_array(jjj,1:BCH_n+1),Staircasemm,BCH_k,Staircasett);
                
            end       
            Blocks_in_Window{1,jj-1}=Temp1_array(:,1:Block_Length);
            Blocks_in_Window{1,jj-1}=Blocks_in_Window{1,jj-1}';
            Blocks_in_Window{1,jj}= Temp1_array(:,Block_Length+1:BCH_n+1); 
        end
    end
    for iii=1:Slidding_Window_Length
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





        



