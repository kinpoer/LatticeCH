function kuc = MatrixAssembly2D(kuc,k,index)
% Multiscale Homogenization Code for Periodic Lattice Metamaterial
% Jin BO, Email:kinpoer@nuaa.edu.cn
%% Assembly of element matrices into the system matrix
eleDof = length(index);
for i=1:eleDof
   ii=index(i);
   for j=1:eleDof
      jj=index(j);
      kuc(ii,jj)=kuc(ii,jj)+k(i,j);
   end
end
end