function k = TimoshenkoBeam2D(prop,L,beta,optSection)
% Multiscale Homogenization Code for Periodic Lattice Metamaterial
% Jin BO, Email:kinpoer@nuaa.edu.cn
%% 2D closed-form Timoshenko beam element stiffness
E  = prop(1);   % Young's modulus
nu = prop(2);   % Poisson's ratio
A  = prop(6);   % cross-sectional area
Iz = prop(8);   % second moment of area about local z-axis
% Shear modulus
if numel(prop) < 10 || isempty(prop(10)) || prop(10) <= 0 || prop(10) == 1
    G = E/(2*(1+nu));
else
    G = prop(10);
end
% ck = 1/kappa
% optSection = 1: rectangular section, kappa = 5/6
% optSection = 2: circular/tube section, kappa ≈ 0.9
if optSection == 1
    ck = 6/5;
else
    ck = 10/9;
end
% Rotation matrix: d_local = T * d_global
c = cos(beta);
s = sin(beta);
T = [ c    s    0    0    0    0;
     -s    c    0    0    0    0;
      0    0    1    0    0    0;
      0    0    0    c    s    0;
      0    0    0   -s    c    0;
      0    0    0    0    0    1];
% Closed-form Timoshenko beam stiffness in local coordinates
% phi = 12EI/(kappa*G*A*L^2) = 12EI*ck/(G*A*L^2)
EI  = E*Iz;
EA_L = E*A/L;
phi = 12*EI*ck/(G*A*L^2);
psi = 1/(1 + phi);
kb = EI*psi/L^3;
k0 = [ EA_L       0                  0             -EA_L       0                  0;
          0    12*kb             6*L*kb              0   -12*kb             6*L*kb;
          0   6*L*kb     (4+phi)*L^2*kb              0  -6*L*kb     (2-phi)*L^2*kb;
      -EA_L       0                  0              EA_L       0                  0;
          0   -12*kb            -6*L*kb              0    12*kb            -6*L*kb;
          0   6*L*kb     (2-phi)*L^2*kb              0  -6*L*kb     (4+phi)*L^2*kb];
k0 = (k0 + k0')/2;
% Transform to global coordinates
k = T' * k0 * T;
k = (k + k')/2;
end
