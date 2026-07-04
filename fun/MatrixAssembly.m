function K = MatrixAssembly(elements)
% Multiscale Homogenization Code for Periodic Lattice Metamaterial
% Jin BO, Email:kinpoer@nuaa.edu.cn
%% Assembly element stiffness matrix
nElements = length(elements);
nDofs = 0;
for i=1:nElements
    nDofs2=max(elements{i}.DOFs);
    if nDofs<nDofs2; nDofs=nDofs2;end
end
K = zeros(nDofs,nDofs);
for i=1:nElements
    nn = elements{i}.DOFs;
    k = elements{i}.K;
    K(nn,nn) = K(nn,nn)+k;
end
end

