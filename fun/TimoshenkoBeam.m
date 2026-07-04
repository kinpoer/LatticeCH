function [k,k0,T]=TimoshenkoBeam(E,Iy,Iz,A,L,G,Jx,Cx,Cy,Cz)
% Multiscale Homogenization Code for Periodic Lattice Metamaterial
% Jin BO, Email:kinpoer@nuaa.edu.cn
%% 3D closed-form Timoshenko beam element stiffness
% ck = 1/kappa.
ck_y = 10/9;
ck_z = 10/9;
% Coordinate transformation
t = TransM(Cx,Cy,Cz);
z = zeros(3,3);
T = [t z z z;
     z t z z;
     z z t z;
     z z z t];
% Axial and torsional stiffness
EA_L = E*A/L;
GJ_L = G*Jx/L;
% Bending about local z-axis: transverse v and rotation theta_z
EIz = E*Iz;
phi_y = 12*EIz*ck_y/(G*A*L^2);
psi_y = 1/(1 + phi_y);
ky = EIz*psi_y/L^3;
% Bending about local y-axis: transverse w and rotation theta_y
EIy = E*Iy;
phi_z = 12*EIy*ck_z/(G*A*L^2);
psi_z = 1/(1 + phi_z);
kz = EIy*psi_z/L^3;
k0 = zeros(12,12);
% Axial stiffness
k0(1,1) =  EA_L;
k0(1,7) = -EA_L;
k0(7,1) = -EA_L;
k0(7,7) =  EA_L;
% Torsional stiffness
k0(4,4)   =  GJ_L;
k0(4,10)  = -GJ_L;
k0(10,4)  = -GJ_L;
k0(10,10) =  GJ_L;
% Bending about local z-axis: v - theta_z
id = [2 6 8 12];
kb_z = [ 12*ky,              6*L*ky,             -12*ky,              6*L*ky;
        6*L*ky,  (4+phi_y)*L^2*ky,        -6*L*ky,  (2-phi_y)*L^2*ky;
        -12*ky,             -6*L*ky,              12*ky,             -6*L*ky;
        6*L*ky,  (2-phi_y)*L^2*ky,        -6*L*ky,  (4+phi_y)*L^2*ky];
k0(id,id) = k0(id,id) + kb_z;
% Bending about local y-axis: w - theta_y
% Sign convention follows the original EulerBeam.m
id = [3 5 9 11];
kb_y = [ 12*kz,             -6*L*kz,             -12*kz,             -6*L*kz;
        -6*L*kz, (4+phi_z)*L^2*kz,         6*L*kz,  (2-phi_z)*L^2*kz;
        -12*kz,              6*L*kz,              12*kz,              6*L*kz;
        -6*L*kz, (2-phi_z)*L^2*kz,         6*L*kz,  (4+phi_z)*L^2*kz];
k0(id,id) = k0(id,id) + kb_y;
k0 = (k0 + k0')/2;
% Transform to global coordinates
k = T' * k0 * T;
k = (k + k')/2;
end

function t = TransM(Cx,Cy,Cz)
%% Establish a deterministic right-handed local coordinate system
ex = [Cx, Cy, Cz];
ex = ex/norm(ex);
% Choose a stable reference vector not parallel to the beam axis
ref = [0, 0, 1];
if abs(dot(ex,ref)) > 0.90
    ref = [0, 1, 0];
end
ey = cross(ref,ex);
ey = ey/norm(ey);
ez = cross(ex,ey);
ez = ez/norm(ez);
% Transformation matrix: d_local = t * d_global
t = [ex; ey; ez];
end