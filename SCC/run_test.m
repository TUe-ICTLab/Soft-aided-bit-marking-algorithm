function run_test(m)
%close all
% Simulation for an SCC code based on an extended BCH 
%m=7;
t=2;

BCH.n=(2^m)-1;          % BCH Codeword length
BCH.k=(2^m)-1-m*t;      % Number of inforamtion bit in BCH codeword 
BCHR=(BCH.k-(BCH.n+1)/2)/((BCH.n+1)/2);
fprintf('Code rate of SCC is %1.3f\n',BCHR);
SCC.bl=(BCH.n+1)/2;     % Length of the SCC block
NBlocks=1e4;

figure(1);grid on;hold on;
set(gca,'Yscale','log');
axis([0.8e-2,2.8e-2,1e-11,1]);
pause(1);
p=[1.4:.1:2.8]*1e-2;
y=zeros(1,NBlocks*SCC.bl^2);
x=zeros(1,NBlocks*SCC.bl^2);
for pp=1:length(p),fprintf('o');end
fprintf('\n');

parfor pp=1:length(p)
    fprintf('x');
    d=randi([0 1],1,NBlocks*SCC.bl*(BCH.k-SCC.bl));
    x=StaircaseEncoding(d,m,t);
    y=bsc(x,p(pp));
    [dhat,yhat]=StaircaseDecoding(y,m,t,9,7);
    % BER without taking into account the last block
    BER(pp)=sum(dhat(1:end-SCC.bl^2)~=d(1:end-SCC.bl^2))/length(d(1:end-SCC.bl^2));
end
semilogy(p,BER,'^-');grid on;hold on;
drawnow


