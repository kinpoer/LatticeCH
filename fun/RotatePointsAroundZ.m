function rotatedPoints = RotatePointsAroundZ(points, angle)
% Multiscale Homogenization Code for Periodic Lattice Metamaterial
% Jin BO, Email:kinpoer@nuaa.edu.cn
%% Rotate Points Around Z axis
rad = deg2rad(angle);
Rz = [cos(rad), -sin(rad),   0;
     sin(rad),  cos(rad),    0;
        0,         0,        1;];
rotatedPoints = (Rz * points')';
end
