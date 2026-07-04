function cornerIndices = FindCubicCorner(points,th)
% Multiscale Homogenization Code for Periodic Lattice Metamaterial
% Jin BO, Email:kinpoer@nuaa.edu.cn
%% Find the corner points of cubic truss cell
if nargin < 2
    th = 1e-3; 
end
[xmin,xmax,ymin,ymax,zmin,zmax] = MaxMinCubic(points);
isCornerPoint = (abs(points(:,1)-xmin)<th | abs(points(:,1)-xmax)<th) & ...
                (abs(points(:,2)-ymin)<th | abs(points(:,2)-ymax)<th) & ...
                (abs(points(:,3)-zmin)<th | abs(points(:,3)-zmax)<th);
cornerIndices = find(isCornerPoint);
end

function [xmin,xmax,ymin,ymax,zmin,zmax] = MaxMinCubic(points)
%% Calculate the bounding box of the cube truss cell
xmin = min(points(:,1));xmax = max(points(:,1));
ymin = min(points(:,2));ymax = max(points(:,2));
zmin = min(points(:,3));zmax = max(points(:,3));
end