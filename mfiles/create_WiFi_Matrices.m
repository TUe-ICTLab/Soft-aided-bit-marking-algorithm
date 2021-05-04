function create_WiFi_Matrices

%[ H, G, Z ] = buildHG( 648, 1/2 );

for R=[1/2,2/3,3/4,5/6]
    for n=[648,1296,1944]
        [ H, Z ] = buildH( n, R );
    end
end


return

function [ H, Z ] = buildH( n, R )
%BUILDH Creates the parity check matrix H (reading the prototype from file)
% for codeword length n and rate R

Zsize = [27, 54, 81];   % Square submatrices available size

col = 24;               % Number of submatrices in a row
if R == 1/2
    row = 12;           % Number of submatrices in a column
elseif R == 2/3
    row = 8;
elseif R == 3/4
    row = 6;
else
    row = 4;
end
R_str   = get_Rstr(R);

% This block will be needed only if the support for different
% codeword size will be implemented
if n == 648
    Z = Zsize(1);
elseif n == 1296
    Z = Zsize(2);
else
    Z = Zsize(3);
end


if(exist(['matrix/H_',num2str(n),'_',R_str,'.mat'],'file'))
    load(['matrix/H_',num2str(n),'_',R_str],'H');
else
    H = zeros(row*Z,col*Z);

    % Prototype matrix
    % If the element (i,j) is >= 0 the corresponding submatrix will be a right column circular
    % shift of the identy matrix by the number of position indicated by the element
    % otherwise the submatrix will be a zero matrix

    protoH = load(['protoH/',num2str(n),'_',R_str]);

    for i = 1:row
        for j = 1:col
            if protoH(i,j) >= 0
                A = eye(Z);
                k = protoH(i,j);
                H((i-1)*Z+1:(i-1)*Z+Z,(j-1)*Z+1:(j-1)*Z+Z) = circshift(A,[0 k]);
            end
        end
    end

    if ~exist('matrix','dir')
        mkdir('matrix');
    end
    H=sparse(logical(H));   % To be used in matlab
    save(['matrix/H_',num2str(n),'_',R_str],'H');     % Store the matrix to save computation time
end

return


