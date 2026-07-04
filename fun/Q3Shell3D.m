function Kt = Q3Shell3D(nodes,elements,Es,vs,t)
% Multiscale Homogenization Code for Periodic Lattice Metamaterial
% Jin BO, Email:kinpoer@nuaa.edu.cn
%% Stiffness of 3D Shell structures
DOFm = size(nodes,1)*6;
Kt = zeros(DOFm,DOFm);
for iElem = 1:size(elements,1)
    N1 = elements(iElem,1); N2 = elements(iElem,2); N3 = elements(iElem,3);
    % local coordinate of element
    local = [N1*6-5 N1*6-4 N1*6-3 N1*6-2 N1*6-1 N1*6 ...
        N2*6-5 N2*6-4 N2*6-3 N2*6-2 N2*6-1 N2*6 ...
        N3*6-5 N3*6-4 N3*6-3 N3*6-2 N3*6-1 N3*6];
    xyz(:,1) = nodes(N1,:);xyz(:,2)=nodes(N2,:);xyz(:,3)=nodes(N3,:);
    xyz(:,4) = (xyz(:,1)+xyz(:,2)+xyz(:,3))./3;
    Vx = xyz(:,2)-xyz(:,4);
    Vy = xyz(:,3)-xyz(:,4);
    Vz = cross(Vx,Vy);
    L = zeros(3,3);
    L(:,1) = Vx./(norm(Vx));L(:,3)=Vz./(norm(Vz));L(:,2)=cross(L(:,3),L(:,1));
    xyz(:,1) = L'*(xyz(:,1)-xyz(:,4));xyz(:,2)=L'*(xyz(:,2)-xyz(:,4));xyz(:,3)=L'*(xyz(:,3)-xyz(:,4));
    % Transformation matrix between local and global coordinate system
    T = blkdiag(L,L,L,L,L,L);
    K0 = TriStiffness(Es,vs,t,xyz);
    Ke = T*K0*T';
    Kt(local,local)=Kt(local,local)+Ke;
end
end

function Km = TriStiffness(E,mu,t,xyz)
%% 3D Q3 Shell element stiffness
% plane stress stiffness
A1 = [1 xyz(1,1) xyz(2,1);1  xyz(1,2) xyz(2,2);1  xyz(1,3) xyz(2,3)];
detJ = det(A1)/2;
Bp(1,1) = xyz(2,2)-xyz(2,3);
Bp(2,2) = -xyz(1,2)+xyz(1,3);
Bp(3,1) = Bp(2,2);
Bp(3,2) = Bp(1,1);
Bp(1,3) = xyz(2,3)-xyz(2,1);
Bp(2,4) = -xyz(1,3)+xyz(1,1);
Bp(3,3) = Bp(2,4);
Bp(3,4) = Bp(1,3);
Bp(1,5) = xyz(2,1)-xyz(2,2);
Bp(2,6) = -xyz(1,1)+xyz(1,2);
Bp(3,5) = Bp(2,6);
Bp(3,6) = Bp(1,5);
Bp = Bp/detJ/2;
Dp = E*t/(1-mu^2)*[1 mu 0;mu 1 0;0 0 (1-mu)/2];
Kp = Bp'*Dp*Bp*detJ;
% plate bending stiffness
Bb(1,3) = xyz(2,2)-xyz(2,3);
Bb(2,2) = xyz(1,2)-xyz(1,3);
Bb(3,2) = -Bb(1,3);
Bb(3,3) = -Bb(2,2);
Bb(1,6) = xyz(2,3)-xyz(2,1);
Bb(2,5) = xyz(1,3)-xyz(1,1);
Bb(3,5) = -Bb(1,6);
Bb(3,6) = -Bb(2,5);
Bb(1,9) = xyz(2,1)-xyz(2,2);
Bb(2,8) = xyz(1,1)-xyz(1,2);
Bb(3,8) = -Bb(1,9);
Bb(3,9) = -Bb(2,8);
Bb = Bb/detJ/2;
Db = E*t^3/12/(1-mu^2)*[1 mu 0;mu 1 0;0 0 (1-mu)/2];
kb = Bb'*Db*Bb*detJ;
% plate shear stiffness
x1 = xyz(1,1);x2 = xyz(1,2);x3 = xyz(1,3);
y1 = xyz(2,1);y2 = xyz(2,2);y3 = xyz(2,3);
A1 = [1 xyz(1,1) xyz(2,1);1 xyz(1,2) xyz(2,2);1 xyz(1,3) xyz(2,3)];
A = det(A1)/2;
m1 = (y2-y3)/2/A;m2 = (y3-y1)/2/A;m3 = (y1-y2)/2/A;
n1 = -(x2-x3)/2/A;n2 = -(x3-x1)/2/A;n3 = -(x1-x2)/2/A;
Bs1 = 1/3*[2*m1-m2-m3 -2*y1*m1-y2*m2/2-y3*m3/2 2*x1*m1+x2*m2/2+x3*m3/2;
    2*n1-n2-n3 -2*y1*n1-y2*n2/2-y3*n3/2 2*x1*n1+x2*n2/2+x3*n3/2;];
Bs2 = 1/3*[-m1+2*m2-m3 -y1*m1/2-2*y2*m2-y3*m3/2 x1*m1/2+2*x2*m2+x3*m3/2;
    -n1+2*n2-n3 -y1*n1/2-2*y2*n2-y3*n3/2 x1*n1/2+2*x2*n2+x3*n3/2;];
Bs3 = 1/3*[-m1-m2+2*m3 -y1*m1/2-y2*m2/2-2*y3*m3 x1*m1/2+x2*m2/2+2*x3*m3;
    -n1-n2+2*n3 -y1*n1/2-y2*n2/2-2*y3*n3 x1*n1/2+x2*n2/2+2*x3*n3;];
Bs = [Bs1 Bs2 Bs3];
Ds = E*t/2/(1+mu)*eye(2)*5/6;
ks = Bs'*Ds*Bs*detJ/2;
Kb = kb+ks;
% shell stiffness
K11 = blkdiag(Kp(1:2,1:2),Kb(1:3,1:3),1e-3);
K12 = blkdiag(Kp(1:2,3:4),Kb(1:3,4:6),0);
K13 = blkdiag(Kp(1:2,5:6),Kb(1:3,7:9),0);
K21 = blkdiag(Kp(3:4,1:2),Kb(4:6,1:3),0);
K22 = blkdiag(Kp(3:4,3:4),Kb(4:6,4:6),1e-3);
K23 = blkdiag(Kp(3:4,5:6),Kb(4:6,7:9),0);
K31 = blkdiag(Kp(5:6,1:2),Kb(7:9,1:3),0);
K32 = blkdiag(Kp(5:6,3:4),Kb(7:9,4:6),0);
K33 = blkdiag(Kp(5:6,5:6),Kb(7:9,7:9),1e-3);
Km = [K11,K12,K13;
    K21,K22,K23;
    K31,K32,K33];
end
