% Author: Alex Alvarado
% Modified: Gabriele Liga Oct 2018, comments added

% This function computes L-values for 
% a given constellation X, 
% labeling L, 
% snrdB based on the channel observations y.
% y is assumed the s

% Ik1,Ik0 are indexes for the subset of symbols in the constellation with
% bit k equal to either 1 or 0, respectively. 

%% Input dimensions
% X : M x N  matrix, where M is the constellation cardinality and N is the number of dimensions
% L : M x m  labeling matrix, where M=2^m is the number of constellation points
% y : Ns x N.input sequence  
% 
% The function returns a  m x Ns matrix of LLRs
% flag_demapper==0 Sum-Exp
% flag_demapper==1 Max-Log

function LLR =  demapper(X,L,snrdB,y,flag_demapper)

snr = 10^(snrdB/10);
N0  = 1/snr;
M   = size(X,1);
m   = size(L,2);
Ns  = size(y,1);

LLR = zeros(m,Ns);
for k=1:m
    Ik0=(find(L(:,k)==0));
    Ik1=(find(L(:,k)==1));
    for l=1:Ns
        if flag_demapper==0 % Sum-Exp
            LLR(k,l) =  log(sum(exp(-sum((repmat(y(l,:),M/2,1)-X(Ik1,:)).^2/N0,2))))-...
                        log(sum(exp(-sum((repmat(y(l,:),M/2,1)-X(Ik0,:)).^2/N0,2))));
        elseif flag_demapper==1 %Max-log
            LLR(k,l) =  log(max(exp(-sum((repmat(y(l,:),M/2,1)-X(Ik1,:)).^2/N0,2))))-...
                        log(max(exp(-sum((repmat(y(l,:),M/2,1)-X(Ik0,:)).^2/N0,2))));
        end
    end
end