function Be = CubicBeMatrix(a1,a2,a3)
% Multiscale Homogenization Code for Periodic Lattice Metamaterial
% Jin BO, Email:kinpoer@nuaa.edu.cn
%% compute Be matrix
if norm(a3)~=0
    Be1 = [a1(1) 0 0 a1(2)/2 0 a1(3)/2;0 a1(2) 0 a1(1)/2 a1(3)/2 0;0 0 a1(3) 0 a1(2)/2 a1(1)/2];
    Be2 = [a2(1) 0 0 a2(2)/2 0 a2(3)/2;0 a2(2) 0 a2(1)/2 a2(3)/2 0;0 0 a2(3) 0 a2(2)/2 a2(1)/2];
    Be3 = [a3(1) 0 0 a3(2)/2 0 a3(3)/2;0 a3(2) 0 a3(1)/2 a3(3)/2 0;0 0 a3(3) 0 a3(2)/2 a3(1)/2];
    Be = [Be1;Be2;Be3];
else
    Be1 = [a1(1) 0 a1(2)/2;0 a1(2) a1(1)/2];
    Be2 = [a2(1) 0 a2(2)/2;0 a2(2) a2(1)/2];
    Be = [Be1;Be2];
end
end

