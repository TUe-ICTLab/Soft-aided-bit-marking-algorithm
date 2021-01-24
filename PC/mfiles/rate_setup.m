function [n,p,k,rp,hEnc,hDec,hError,f_punct,n_pct]=rate_setup(R)

f_punct=0;      % Flag: punctured or not?
n_pct=0;        % Number of bits to be punctured
if R==0.42 % Punctured 2/5 code (n,k)=(61800,25920). OH=138.4%
    H       = dvbs2ldpc(2/5);
    f_punct = 1;
    n_pct   = 3000;
elseif R==0.45 % Punctured 2/5 code (n,k)=(57600,25920). OH=122.2%
    H       = dvbs2ldpc(2/5);
    f_punct = 1;
    n_pct   = 7200;
elseif R==0.47 % Punctured 2/5 code (n,k)=(55150,25920). OH=112.7%
    H       = dvbs2ldpc(2/5);
    f_punct = 1;
    n_pct   = 9650;
elseif R==0.53 % Punctured 1/2 code (n,k)=(60800,32400). OH=87.7%
    H       = dvbs2ldpc(1/2);
    f_punct = 1;
    n_pct   = 4000;
elseif R==0.57 % Punctured 1/2 code (n,k)=(57250,32400). OH=76.7%
    H       = dvbs2ldpc(1/2);
    f_punct = 1;
    n_pct   = 7550;
elseif R==0.62 % Punctured 3/5 code (n,k)=(62700,38880). OH=61.3%
    H       = dvbs2ldpc(3/5);
    f_punct = 1;
    n_pct   = 2100;
elseif R==0.64 % Punctured 3/5 code (n,k)=(60750,38880). OH=56.3%
    H       = dvbs2ldpc(3/5);
    f_punct = 1;
    n_pct   = 4050;
elseif R==0.69 % Punctured 2/3 code (n,k)=(63000,43200). OH=45.8%
    H       = dvbs2ldpc(2/3);
    f_punct = 1;
    n_pct   = 1800;
elseif R==0.71 % Punctured 2/3 code (n,k)=(61000,43200). OH=41.2%
    H       = dvbs2ldpc(2/3);
    f_punct = 1;
    n_pct   = 3800;
elseif R==0.72 % Punctured 2/3 code (n,k)=(60000,43200). OH=38.9%
    H       = dvbs2ldpc(2/3);
    f_punct = 1;
    n_pct   = 4800;
elseif R==0.77 % Punctured 3/4 code (n,k)=(63000,48600). OH=29.6%
    H       = dvbs2ldpc(3/4);
    f_punct = 1;
    n_pct   = 1800;
elseif R==0.81 % Punctured 4/5 code (n,k)=(64000,51840). OH=23.4%
    H       = dvbs2ldpc(4/5);
    f_punct = 1;
    n_pct   = 800;
elseif R==0.82 % Punctured 4/5 code (n,k)=(63000,51840). OH=21.5%
    H       = dvbs2ldpc(4/5);
    f_punct = 1;
    n_pct   = 1800;
elseif R==0.84          % Punctured 5/6 code (n,k)=(62000,54000). OH=14.8%
    H       = dvbs2ldpc(5/6);
    f_punct = 1;
    n_pct   = 800;
elseif R==0.86          % Punctured 5/6 code (n,k)=(63000,54000). OH=16.7%
    H       = dvbs2ldpc(5/6);
    f_punct = 1;
    n_pct   = 1800;
elseif R==0.87      % Punctured 5/6 code (n,k)=(64000,54000). OH=18.5%
    H       = dvbs2ldpc(5/6);
    f_punct = 1;
    n_pct   = 2800;
elseif R==0.91      % Punctured 9/10 code (n,k)=(64320,58320). OH=10.29%
    H       = dvbs2ldpc(9/10);
    f_punct = 1;
    n_pct   = 480;
elseif R==0.92      % Punctured 9/10 code (n,k)=(63320,58320). OH=8.57%
    H       = dvbs2ldpc(9/10);
    f_punct = 1;
    n_pct   = 1480;
elseif R==0.94      % Punctured 9/10 code (n,k)=(62320,58320). OH=6.86%
    H       = dvbs2ldpc(9/10);
    f_punct = 1;
    n_pct   = 2480;
else            % Regular (unpunctured) DVB-S2 code
    H       = dvbs2ldpc(R);
end
n       = size(H,2);    % n
p       = size(H,1);    % n-k
k       = n-p;          % k
rp      = randperm(p);  % random bits to be punctured. This in principle could be optimzied and saved.
hEnc    = comm.LDPCEncoder(H);
if gpuDeviceCount==0
    hDec    = comm.LDPCDecoder(H);
else
    hDec    = comm.gpu.LDPCDecoder(H);
end
hError  = comm.ErrorRate;