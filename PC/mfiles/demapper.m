function LLR =  demapper(X,L,snrdB,y,Ik1,Ik0,flag_demapper)
% This function computes L-values for a given constellation X, labeling
% L, snr based on the channel observations Y.
% X is M times N, where N is the number of dimensions
% L is M times m, where M=2^m is the number of constellation points
% y is Ns times N.
%
% The function returns a matrix of m times Ns LLRs
%
% flag_demapper==0 Sum-Exp
% flag_demapper==1 Max-Log

snr = 10^(snrdB/10);
N0  = 1/snr;
M   = size(X,1);
m   = size(L,2);
Ns  = size(y,1);

LLR = zeros(m,Ns);
for k=1:m
    for l=1:Ns
        if flag_demapper==0 % Sum-Exp
            LLR(k,l) =  log(sum(exp(-sum((repmat(y(l,:),M/2,1)-X(Ik1(:,k),:)).^2/N0,2))))-...
                        log(sum(exp(-sum((repmat(y(l,:),M/2,1)-X(Ik0(:,k),:)).^2/N0,2))));
        elseif flag_demapper==1 %Max-log
            LLR(k,l) =  log(max(exp(-sum((repmat(y(l,:),M/2,1)-X(Ik1(:,k),:)).^2/N0,2))))-...
                        log(max(exp(-sum((repmat(y(l,:),M/2,1)-X(Ik0(:,k),:)).^2/N0,2))));
        end
    end
end