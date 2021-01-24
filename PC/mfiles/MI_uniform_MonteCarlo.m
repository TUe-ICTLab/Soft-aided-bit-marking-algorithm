function I=MI_uniform_MonteCarlo(X,Y,N0)
% I=MI_uniform_MonteCarlo(X,Y,N0)
% This function computes the mutual information (in bits) for arbitrary
% N-dimensional input and output vectors (x and y) via Monte-Carlo
% integration. The function works only for an N-dimensional Gaussian
% channel Y=X+Z where E[|X|^2]=1, Z=[Z1,...,ZN], where Zn are i.i.d.
% Gaussian random variables with zero mean and variance N0/2. It also 
% assumes a uniform input distribution on the symbols
%
% X is a Ns x N vector of transmitted symbols
% Y is a Ns x N vector of transmitted symbols
% N0 is twice the variance of the noise in each dimension
%
% Alex Alvarado
% Sep. 2015

C=unique(X,'rows');         % Constellation
M=size(C,1);                % Constellation size
m=log2(M);                  % Bits/symbol
Ns=size(X,1);            	% Number of samples
Z=Y-X;                      % This should be a matrix with i.i.d. N(0,N0/2)
I=0;
for i=1:M
    for l=1:Ns
        arg=0;
        % In principle the next if is not needed. But we want to calculate
        % the "true" MI passed to the decoder.
        if isequal(X(l,:),C(i,:))   
            for j=1:M
                dij = C(i,:)-C(j,:);
                arg = exp(-(2*sum(Z(l,:).*dij)+sum(dij.^2))/N0)+arg;
            end
            I       = -log2(arg) + I;
        end
    end
end
I=m+I/Ns*1/M*M;
% Note that the last *M comes from the if above, which makes the number of
% samples to go down by a factor of M.
return