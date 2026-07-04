function mapping = MatchCubicBoundary(nodes,cornerIdx,edgeIdx,surfaceIdx,th)
% Multiscale Homogenization Code for Periodic Lattice Metamaterial
% Jin BO, Email:kinpoer@nuaa.edu.cn
%% Matching indices of node by PBC 
if nargin < 5
    th = 1e-3; 
end
[xmin,xmax,ymin,ymax,zmin,zmax] = MaxMinCubic(nodes);
% Different types of node indices
cornerindices = cornerIdx;
edgeindices = setdiff(edgeIdx,cornerIdx);
surfaceindices = setdiff(surfaceIdx,edgeIdx);
innerindices = setdiff([1:length(nodes)]',surfaceIdx);

if max(nodes(:,3))==min(nodes(:,3))
    innerindices = setdiff([1:length(nodes)]',edgeIdx);
end
mapping = zeros(length(nodes),2);
mapping(:,1) = 1:length(nodes);

% Determine whether different types of nodes exist
isexistcorner = ~isempty(cornerindices);
isexistedge = ~isempty(edgeindices);
isexistsurface = ~isempty(surfaceindices);
isexistinner = ~isempty(innerindices);

% First distinguish between independent points and mapping points
% ----------- corner ----------- %
% The 8 corner 
if isexistcorner
cornerindependent = find(abs(nodes(:,1)-xmax)<th & abs(nodes(:,2)-ymax)<th & abs(nodes(:,3)-zmin)<th);
RotatePointsAroundZ(nodes,pi/2);
cornerother = setdiff(cornerindices,cornerindependent);
mapping(cornerindependent,2)=cornerindependent;
mapping(cornerother,2)=cornerindependent;
end

% ----------- edge ----------- %
% The 12 edges are divided into three categories according to their directions
if isexistedge
    points = nodes(edgeindices,:);
    % find points with x=xmax and z=zmin 
    index = intersect(find(abs(points(:,1)-xmax)<th),find(abs(points(:,3)-zmin)<th));
    matcha1 =zeros(length(index),1);
    for i=1:length(index)
        matcha1(i,1) = edgeindices(index(i));
        mapping(matcha1(i,1),2) = matcha1(i,1);
        tempoints = nodes(edgeindices(index(i)),:);
        for j=1:length(nodes)
            if abs(nodes(j,2)-tempoints(1,2))<th && abs(nodes(j,1)-xmin)<th && abs(nodes(j,3)-zmin)<th
                mapping(j,2) = matcha1(i,1);
            elseif abs(nodes(j,2)-tempoints(1,2))<th && abs(nodes(j,3)-zmax)<th && abs(nodes(j,1)-tempoints(1,1))<th
                mapping(j,2) = matcha1(i,1);
            elseif abs(nodes(j,1)-xmin)<th && abs(nodes(j,2)-tempoints(1,2))<th && abs(nodes(j,3)-zmax)<th
                mapping(j,2) = matcha1(i,1);
            end
        end
    end

    % find points with y=yamx and z=zmin
    index = intersect(find(abs(points(:,2)-ymax)<th),find(abs(points(:,3)-zmin)<th));
    matcha2 =zeros(length(index),1);
    for i=1:length(index)
        matcha2(i,1) = edgeindices(index(i));
        mapping(matcha2(i,1),2) = matcha2(i,1);
        tempoints = nodes(edgeindices(index(i)),:);
        for j=1:length(nodes)
            if abs(nodes(j,1)-tempoints(1,1))<th && abs(nodes(j,2)-ymin)<th && abs(nodes(j,3)-zmin)<th
                mapping(j,2) = matcha2(i,1);
            elseif abs(nodes(j,1)-tempoints(1,1))<th && abs(nodes(j,3)-zmax)<th && abs(nodes(j,2)-tempoints(1,2))<th
                mapping(j,2) = matcha2(i,1);
            elseif abs(nodes(j,2)-ymin)<th && abs(nodes(j,1)-tempoints(1,1))<th && abs(nodes(j,3)-zmax)<th
                mapping(j,2) = matcha2(i,1);
            end
        end
    end

    % find points with x=xmax and y=yamx 
    index = intersect(find(abs(points(:,1)-xmax)<th),find(abs(points(:,2)-ymax)<th));
    matcha3 =zeros(length(index),1);
    for i=1:length(index)
        matcha3(i,1) = edgeindices(index(i));
        mapping(matcha3(i,1),2) = matcha3(i,1);
        tempoints = nodes(edgeindices(index(i)),:);
        for j=1:length(nodes)
            if abs(nodes(j,1)-tempoints(1,1))<th && abs(nodes(j,2)-ymin)<th && abs(nodes(j,3)-tempoints(1,3))<th
                mapping(j,2) = matcha3(i,1);
            elseif abs(nodes(j,2)-tempoints(1,2))<th && abs(nodes(j,1)-xmin)<th && abs(nodes(j,3)-tempoints(1,3))<th
                mapping(j,2) = matcha3(i,1);
            elseif abs(nodes(j,1)-xmin)<th && abs(nodes(j,3)-tempoints(1,3))<th && abs(nodes(j,2)-ymin)<th
                mapping(j,2) = matcha3(i,1);
            end
        end
    end
end
% ----------- surface ----------- %
% The 6 surfaces are divided into three categories according to their directions
if isexistsurface
    points = nodes(surfaceindices,:);
    % find points on plane of x=xmax
    index = find(abs(points(:,1)-xmax)<th);
    matchface1 =zeros(length(index),1);
    for i=1:length(index)
        matchface1(i,1) = surfaceindices(index(i));
        mapping(matchface1(i,1),2) = matchface1(i,1);
        tempoints = nodes(surfaceindices(index(i)),:);
        for j=1:length(nodes)
            if abs(nodes(j,2)-tempoints(1,2))<th && abs(nodes(j,1)-xmin)<th && abs(nodes(j,3)-tempoints(1,3))<th
                mapping(j,2) = matchface1(i,1);
            end
        end
    end
    % find points on plane of y=ymax
    index = find(abs(points(:,2)-ymax)<th);
    matchface2 =zeros(length(index),1);
    for i=1:length(index)
        matchface2(i,1) = surfaceindices(index(i));
        mapping(matchface2(i,1),2) = matchface2(i,1);
        tempoints = nodes(surfaceindices(index(i)),:);
        for j=1:length(nodes)
            if abs(nodes(j,1)-tempoints(1,1))<th && abs(nodes(j,2)-ymin)<th && abs(nodes(j,3)-tempoints(1,3))<th
                mapping(j,2) = matchface2(i,1);
            end
        end
    end
    % find points on plane of z=zmin
    index = find(abs(points(:,3)-zmin)<th);
    matchface3 =zeros(length(index),1);
    for i=1:length(index)
        matchface3(i,1) = surfaceindices(index(i));
        mapping(matchface3(i,1),2) = matchface3(i,1);
        tempoints = nodes(surfaceindices(index(i)),:);
        for j=1:length(nodes)
            if abs(nodes(j,1)-tempoints(1,1))<th && abs(nodes(j,3)-zmax)<th && abs(nodes(j,2)-tempoints(1,2))<th
                mapping(j,2) = matchface3(i,1);
            end
        end
    end
end
% ----------- inner ----------- %
if isexistinner
mapping(innerindices,2)=innerindices;
end
end

function [xmin,xmax,ymin,ymax,zmin,zmax] = MaxMinCubic(points)
%% Calculate the bounding box of the cube truss cell
xmin = min(points(:,1));xmax = max(points(:,1));
ymin = min(points(:,2));ymax = max(points(:,2));
zmin = min(points(:,3));zmax = max(points(:,3));
end

function rotatedPoints = RotatePointsAroundZ(points, angle)
%% Rotate the points around the z-axis by a certain angle
rad = deg2rad(angle);
Rz = [
    cos(rad), -sin(rad), 0;
    sin(rad),  cos(rad), 0;
    0,         0,        1;];
rotatedPoints = (Rz * points')';
end
