% This scriupt runs some simple simulations of a product code over a BSC
% Note that the BCH decoder is compiled for windows so it will only work in
% a windows machine...

startup;
clear all;
%close all
L=10; % iterations
pvec=[1.6:-0.05:1.6]*1e-2;
BCHnu=7;
BCHt=2;
BCHe=1;
BCHn=(2^BCHnu)-1;                   % Coded bits, without the extended bits
BCHk=(2^BCHnu)-1-BCHnu*BCHt;        % Information bits
for pp=1:length(pvec)
    fprintf('x');
   parfor ii=1:2e3 
        info_bits=randi([0,1],BCHk,BCHk);                               % Generate the info bit matrix
        coded_bits=PC_BCH_Encoder(info_bits,BCHnu,BCHt,BCHe);           % Encode the data using a product code
        coded_bits_bsc=bsc(coded_bits,pvec(pp));                        % Transmit through a BSC channel
        hat_coded_bits=PC_BCH_Decoder(coded_bits_bsc,BCHnu,BCHt,BCHe,L);% Decode
        BERout(pp,ii)=sum(sum(hat_coded_bits(1:BCHk,1:BCHk)~=coded_bits(1:BCHk,1:BCHk)))/(BCHk^2); % Count errors
    end
    fprintf('p=%f, BER=%f \n',pvec(pp),mean(BERout(pp,:)));
end
figure(1);
semilogy(pvec,mean(BERout,2),'r*-');
grid on;hold on;
axis([1e-2*0.9,2.7e-2,1e-10,1e-1]);


