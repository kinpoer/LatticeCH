function surfaceIndices = FindCubicSurface(points,th)
% Multiscale Homogenization Code for Periodic Lattice Metamaterial
% Jin BO, Email:kinpoer@nuaa.edu.cn
%% Find the surface points of cubic truss cell
if nargin < 2
    th = 1e-3; 
end
[minx,maxx,miny,maxy,minz,maxz] = MaxMinCubic(points);
surfaceIndices = [];
for i = 1:size(points, 1)
    p = points(i, :);
    if abs(p(1)-minx)<th || abs(p(1)-maxx)<th || abs(p(2)-miny)<th || abs(p(2)-maxy)<th || abs(p(3)-minz)<th || abs(p(3)-maxz)<th
        surfaceIndices = [surfaceIndices; i]; 
    end
end
end

function [xmin,xmax,ymin,ymax,zmin,zmax] = MaxMinCubic(points)
%% Calculate the bounding box of the cube truss cell
xmin = min(points(:,1));xmax = max(points(:,1));
ymin = min(points(:,2));ymax = max(points(:,2));
zmin = min(points(:,3));zmax = max(points(:,3));
end