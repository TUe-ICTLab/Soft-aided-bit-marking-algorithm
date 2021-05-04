function L=Get_Labeling(m,type)

% L=Get_Labeling(m,type);
% Create the labeling L (M x m)
% type: 'BRGC', 'NBC', 'AGC', and 'FBC'

if m==1
    L=[0 1]';
else
    M=2^m;
    switch type
        case 'BRGC'
            L=zeros(M,m);
            L(1:M/2,2:m)=Get_Labeling(m-1,type);
            L(M/2+1:M,2:m)=flipud(L(1:M/2,2:m));
            L(M/2+1:M,1)=1;
        case 'NBC'
            L=fliplr(de2bi(0:M-1));
        case 'AGC'
            L=zeros(M,m);
            L(1:M/2,2:m)=Get_Labeling(m-1,type);
            L(M/2+1:M,2:m)=flipud(L(1:M/2,2:m));
            L(2:2:M,1)=1;
            L(M/2+1:M,:)=not(L(M/2+1:M,:));
        case 'FBC'
            L=zeros(M,m);
            L(1:M/2,2:m)=Get_Labeling(m-1,'NBC');
            L(M/2+1:M,2:m)=flipud(L(1:M/2,2:m));
            L(M/2+1:M,1)=1;
        otherwise
            error(sprintf('Only ''BRGC'', ''NBC'', ''AGC'' and ''FBC'' are supported and here type=''%s''',type))
    end
end