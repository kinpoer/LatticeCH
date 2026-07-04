function B0 = GetB0(B0tran,type)
% Multiscale Homogenization Code for Periodic Lattice Metamaterial
% Jin BO, Email:kinpoer@nuaa.edu.cn
%% Get the B0 matrix through B0tran
[m,n] = size(B0tran);
switch type
    case '2D'
        k = 3;
        I = eye(3);
        O = zeros(3,3);
        B0 = zeros(m*k,n*k);
    case '3D'
        k = 6;
        I = eye(6);
        O = zeros(6,6);
        B0 = zeros(m*k,n*k);
end
for i=1:m
    for j=1:n
        if B0tran(i,j)==1
            B0((i-1)*k+1:(i-1)*k+k,(j-1)*k+1:(j-1)*k+k) = I;
        elseif B0tran(i,j)==0
            B0((i-1)*k+1:(i-1)*k+k,(j-1)*k+1:(j-1)*k+k) = O;
        end
    end
end
end
