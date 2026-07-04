function [k,k1,T] = EulerBeam(E,Iy,Iz,A,L,G,J,Cx,Cy,Cz)
% Multiscale Homogenization Code for Periodic Lattice Metamaterial
% Jin BO, Email:kinpoer@nuaa.edu.cn
%% 3D Euler beam elements stiffness
k1 = [(E*A)/L        0             0              0        0             0          (-E*A)/L            0                  0              0             0              0;... 
        0     (12*E*Iz)/L^3       0              0        0         (6*E*Iz)/L^2       0             (-12*E*Iz)/L^3       0              0             0         (6*E*Iz)/L^2;...
        0           0       (12*E*Iy)/L^3        0   (-6*E*Iy)/L^2       0             0                0           (-12*E*Iy)/L^3       0        (-6*E*Iy)/L^2       0;...
        0           0             0           (G*J)/L     0              0             0                0                 0           (-G*J)/L         0              0;...
        0           0        (-6*E*Iy)/L^2       0     (4*E*Iy)/L        0             0                0            (6*E*Iy)/L^2        0         (2*E*Iy)/L         0;...
        0      (6*E*Iz)/L^2       0              0        0         (4*E*Iz)/L         0           (-6*E*Iz)/L^2          0              0             0         (2*E*Iz)/L;...
    (-E*A)/L        0             0              0        0              0           (E*A)/L            0                 0              0             0              0;... 
	    0     (-12*E*Iz)/L^3      0              0        0        (-6*E*Iz)/L^2       0            (12*E*Iz)/L^3         0              0             0         (-6*E*Iz)/L^2;... 
        0           0       (-12*E*Iy)/L^3       0   (6*E*Iy)/L^2        0             0                0            (12*E*Iy)/L^3       0         (6*E*Iy)/L^2       0;... 
        0           0             0          (-G*J)/L     0              0             0                0                 0           (G*J)/L          0              0;... 
        0           0       (-6*E*Iy)/L^2        0     (2*E*Iy)/L        0             0                0             (6*E*Iy)/L^2       0         (4*E*Iy)/L         0;...
        0      (6*E*Iz)/L^2       0              0        0         (2*E*Iz)/L         0          (-6*E*Iz)/L^2            0             0             0         (4*E*Iz)/L];
% Arbitrary local coordinates
t = TransM(Cx,Cy,Cz);
z = zeros(3,3);
T = [t z z z; z t z z; z z t z; z z z t];
k = T'*k1*T ;
end

function t = TransM(Cx,Cy,Cz)
%% Establish a local coordinate system to calculate the transformation matrix
v = [Cx,Cy,Cz];
v = v/norm(v);
% Generate a random vector
r = rand(3,1);
% Make sure the random vector not parallel to the beam direction
while dot(v,r) == 1
    r = rand(3, 1);
end
v1 = cross(v,r);v1 = v1/norm(v1); 
v2 = cross(v,v1);v2 = v2/norm(v2);
crossProduct = cross(v,v1);
dotProduct = dot(crossProduct,v2);
if dotProduct > 0
    isRightHanded = true;
else
    isRightHanded = false;
end
if isRightHanded == 0
    v1 = -v1;
end
% Calculate the direction cosines of V1 and V2 and the global coordinate system
cosThetaX1 = v1(1)/norm(v1);
cosThetaY1 = v1(2)/norm(v1);
cosThetaZ1 = v1(3)/norm(v1);
cosThetaX2 = v2(1)/norm(v2);
cosThetaY2 = v2(2)/norm(v2);
cosThetaZ2 = v2(3)/norm(v2);
% Transformation Matrix
t = [Cx,Cy,Cz;cosThetaX1,cosThetaY1,cosThetaZ1;cosThetaX2,cosThetaY2,cosThetaZ2];
end