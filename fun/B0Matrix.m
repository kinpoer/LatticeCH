function B0 = B0Matrix(mat,mapping)
% Multiscale Homogenization Code for Periodic Lattice Metamaterial
% Jin BO, Email:kinpoer@nuaa.edu.cn
%% Calculate the B0 matrix
nodes = mat.NewNode;
dependentNodeNum = length(find(mapping(:,1) == mapping(:,2)));
B0tran = zeros(length(nodes),dependentNodeNum);
for i=1:length(mapping)
    B0tran(i,mapping(i,2))=1;
end
nonZeroCols = any(B0tran~= 0, 1);
B0tran = B0tran(:,nonZeroCols);
if norm(mat.a3)==0
B0 = GetB0(B0tran,'2D');
else
B0 = GetB0(B0tran,'3D');
end
end

