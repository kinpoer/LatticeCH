function index = ElementDof(nd,noNel,noDof)
% Multiscale Homogenization Code for Periodic Lattice Metamaterial
% Jin BO, Email:kinpoer@nuaa.edu.cn
%% Compute system dofs associated with each element
k = 0;
for i=1:noNel
    start = (nd(i)-1)*noDof;
    for j=1:noDof
        k = k+1;
        index(k) = start+j;
    end
end
end
