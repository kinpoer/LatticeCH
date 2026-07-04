function k = EulerBeam2D(prop,L,beta,optSection)
% Multiscale Homogenization Code for Periodic Lattice Metamaterial
% Jin BO, Email:kinpoer@nuaa.edu.cn
%% 2D Euler beam elements stiffness
% Element Properties
E = prop(1);% elastic modulus
u = prop(2);% Poisson's ratio
if optSection==1 % 'Circle'
    h = prop(4);% height of beam cross-section
    b = prop(5);% width of beam cross-section
elseif optSection==2 % 'Tube'
    D = prop(4);% outer diameter of beam cross-section
    d = prop(5);% inner diameter of beam cross-section
end
A = prop(6);% area of beam cross-section
Iy = prop(7);% 2nd moment of inertia of cross-section about axis y
Iz = prop(8);% 2nd moment of inertia of cross-section about axis z
% Rotation matrix for the coordinate transformation
c = cos(beta);s= sin(beta);
T=[ c    s    0    0    0    0;
    -s   c    0    0    0    0;
    0    0    1    0    0    0;
    0    0    0    c    s    0;
    0    0    0   -s    c    0;
    0    0    0    0    0    1];
% Stiffness matrix at the local axis
ka = E*A/L; kc=E*Iz/(L^3);
k0 = [ka   0            0       -ka     0                0;
     0   12*kc        6*L*kc     0    -12* kc       6*L*kc;
     0   6*L*kc    4*L^2*kc      0    -6*L*kc     2*L^2*kc;
    -ka   0            0         ka     0                0;
     0  -12*kc       -6*L*kc     0     12*kc        -6*L*kc;
     0   6*L*kc    2*L^2*kc      0    -6*L*kc    4*L^2*kc];
% Stiffness matrix at the global axis
k=T'*k0*T;
k = (k + k')/2;
end