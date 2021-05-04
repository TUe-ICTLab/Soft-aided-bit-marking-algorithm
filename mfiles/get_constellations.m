function [X,N,L] = get_constellations(M)
% This function returns the constellation points for M-PAM with the BRGC
%
% Alex Alvarado
% June 2017

m=log2(M);
L=[];
X       = MPAM(M);
Lb      = Get_Labeling(m,'BRGC');
L       = bin2dec(num2str(Lb));
X = normalize_Es(X,M);
N = size(X,2);

%plot(X(:,1),X(:,2),'*');axis square;for i=1:M,text(X(i,1),X(i,2),dec2bin(L(i),log2(M)),'FontSize',16),end;
%for i=1:M,text(X(i,1),X(i,2),num2str(i),'FontSize',16),end;


function X = normalize_Es(X,M)

Es=1/M*sum(sum(X.^2,2));
X=X./sqrt(Es);

function X = MPAM(M)

X=[-(M-1):2:M-1].';


