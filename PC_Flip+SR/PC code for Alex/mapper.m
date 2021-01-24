function symbols = mapper(X,bits,m,Ns,N)
% This function maps bits to symbols.
% For 4D constellations, use mapper_mex (which assumes N=4).
% For other constellations, use mapper.m (it's fast anyway).
% This function assumes L is the NBC so that it can be used as pointers
%
% Alex Alvarado
% July 2013

symbols     = zeros(Ns,N);
bits_mat    = reshape(bits,m,Ns).';

for n=1:Ns
    idx = sum(bits_mat(n,:) .* 2.^(0:m-1));   % Gray labeling
    symbols(n,:) = X(idx+1,:);
end
return