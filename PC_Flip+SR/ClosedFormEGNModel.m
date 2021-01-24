%% Returns Optical TX Distance vs EbNo for a given code rate r and given system parameters  
function [DistOut, DistPP, EbNoOut,EbNoPP]= ClosedFormEGNModel(EbNoIn, EbNoPinPoint,r) 


h=6.6256e-34;
c=3e5;
RefWavelength=1550;
Fc=3e8/(RefWavelength*1e-9);            % Reference carrier frequency [Hz]
alpha=0.2;                         % Fibre attenuation [dB/km]
D=17;                               % Dispersion parameter [ps/nm/km]
beta2=-D*RefWavelength.^2/(2*pi*c);% Dispersion coefficient [ps2/km]
Gamma=1.2;                          % Nonlinear parameter [1/W/km]
NFdB=4.5;                           % EDFA noise figure [dB] 
SpanLength=80;                   % Span Length [km]
GdB=SpanLength*alpha;               % EDFA gain [dB] 
Rs=33e9;                            % Symbol rate [sym/s]
DF=33e9;                            % Channel spacing [Hz]
NF=10^(NFdB/10);                    % NF [dB]--> L.U.
G=10^(GdB/10);                      % Gain [dB]--> L.U.
Nsp=(NF*G)/(2*(G-1));               % Spont. emission factor
Pase=2*Nsp*(G-1)*h*Fc*Rs;           % ASE noise power (over 2 pols) [W] 
Leff=1/(alpha/4.3429);              % Fibre Effective length [km]
Nch=295;                              % Number of WDM channels
Ns=1:250;                           % Number of spans
Plch=(-10:.01:12);

%EbNo=3:.1:5;     % Eb/No 
%r=0.78;          % Code rate  

X=[1+1i,1-1i,-1+1i,-1-1i]; % QPSK

%Eta=8/27*Gamma^2*Leff^2/(pi*abs(beta2*(1e-12)^2)*Rs^2*Leff)*asinh(pi^2/2*abs(beta2*(1e-12)^2)*Leff*Rs^2*Nch^(2*Rs/DF));  % Nonlinear coeffiecient Eta [1/W^2]
%eps=3/10*log(1+(6/SpanLength)*(Leff/2/asinh(pi^2/2*abs(beta2)*(1e-12)^2*Rs^2*Leff/2*Nch^(2*Rs/DF))));
Eps1=3/10*log(1+(6/SpanLength)*(Leff/asinh(pi^2/2*abs(beta2)*(1e-12)^2*Rs^2*Leff*Nch^(2*Rs/DF))));
psi=Ns.^(2+Eps1)/(2+Eps1)+Ns.^(1+Eps1)/2;
k=0.7715*sqrt(2+Eps1);
EtaGN=8/27*Gamma^2*Leff^2/(pi*abs(beta2*(1e-12)^2)*Rs^2*Leff)*asinh(pi^2/2*abs(beta2*(1e-12)^2)*Leff*Rs^2*Nch^(2*Rs/DF))*Ns.^(1+Eps1);  % Nonlinear coeffiecient Eta [1/W^2]


Phi=2-mean(abs(X).^4)/mean(abs(X).^2)^2; % QPSK

EtaCorr=80/81*Phi*(Gamma*Leff)^2.*Ns/(Rs*DF*pi*abs(beta2)*(1e-12)^2*SpanLength)*sum(1./(1:(Nch-1)/2));
EtaEGN=EtaGN-EtaCorr;

%EtaEGN=EtaGN;
Save=0;
Case='SNRvsDist'; %'Gain'
%Case='Gain';


%% SNR vs Dist 
P=10.^((Plch-30)/10);
SNR=zeros(length(Ns),length(Plch));
SNRopt=zeros(1,length(Ns));

for ss=1:length(Ns)
SNR(ss,:)=10*log10(P)-10*log10(Ns(ss)*Pase+EtaEGN(ss)*P.^3+3*EtaEGN(1)*Pase*P.^2*psi(ss));  
%SNR(ss,:)=10*log10(P)-10*log10(Ns(ss)*Pase+EtaEGN(ss)*P.^3);  
SNRopt(ss)=max(SNR(ss,:));   
end

EbNo=SNRopt-3-10*log10(r);     % EbNo for each tx distance  
Dist=Ns*SpanLength;


[X,Y]=meshgrid(EbNo,EbNoIn);
[~,idx]=min(abs(X-Y),[],2);
EbNoOut=EbNo(idx);
DistOut=Dist(idx);


tempEbNo=repmat(EbNoPinPoint.',1,length(EbNo));
[~,idxpp]=min(abs(EbNo-tempEbNo),[],2);


EbNoPP=EbNo(idxpp);
DistPP=Dist(idxpp);

[X,Y]=meshgrid(EbNoPP,EbNoIn);
[~,idx]=min(abs(X-Y),[],2);
%EbNoOut=EbNo(idx);

%% Saving Data (both mat files and txt files for (pgf plot))
s_path='/home/gabriele/Dropbox/Apps/Overleaf/ECOC2019_A novel soft-aided marking decoder for product codes/Figures/';

  if ~exist(s_path,'dir')
  mkdir(s_path);
  end
  
   if Save==1
    A=DistOut;  
    %B=[Plch;SNR_DBP];  
    filename=['DistMapEbNo' num2str(Nch) 'ch_r=' num2str(r) '.txt'];
    %filename2=['SNRvsP_DBP_' num2str(Nch) 'ch_' num2str(Ns) 'Spans.txt'];
    %matlab2tikz([s_path filename]);
    fileID1=fopen([s_path filename],'w');   
    disp('Saving...');
    fprintf(fileID1,'%6.2f\n',A);
    %fprintf(fileID2,'%6.2f %12.8f\n',B);
   elseif Save==2
    A=[DistPP;EbNoPP]; 
    %B=[Dist;EbNo];
    filename=['Dist_vs_EbNo_PinPoint' num2str(Nch) 'ch_r=' num2str(r) '.txt'];
    %filename2=['SNRvsP_DBP_' num2str(Nch) 'ch_' num2str(Ns) 'Spans.txt'];
    %matlab2tikz([s_path filename]);
    fileID1=fopen([s_path filename],'w');   
    disp('Saving...');
    fprintf(fileID1,'%6.2f  %12.8f\n',A);
   
   end
end
%p1=plot(Dist,SNRopt,'r','LineWidth',2);
%figure

% plot(Dist,EbNo,'b','LineWidth',2);
% ylim([3 5]);



