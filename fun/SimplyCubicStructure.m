function [nodes1, elements1] = SimplyCubicStructure(nodes,elements)
% Multiscale Homogenization Code for Periodic Lattice Metamaterial
% Jin BO, Email:kinpoer@nuaa.edu.cn
%% Preserve the primitive form of cubic truss cells
centerP = lineCenterPosition(nodes,elements);
[xmin,xmax,ymin,ymax,zmin,zmax] = MaxMinCubic(nodes);
[onZ, onX, onY] = checkPlanes(centerP,xmin,xmax,ymin,ymax,zmin,zmax);
onPlane = zeros(size(centerP,1),1);
if zmin==zmax
    for i=1:length(centerP)
        if centerP(i,1)==xmin || centerP(i,2)==ymin
            onPlane(i,1) = 1;
        end
    end
    onPlane = logical(onPlane);
else
    onPlane = (onZ+onX+onY)> 0;
end
% Remove related nodes and elements
tempElements = elements;tempElements(onPlane,:) = [];
% Update nodes and elements
[nodes1, elements1] =removeIsolatedNodes(nodes,tempElements);
centerP = lineCenterPosition(nodes1,elements1);
[uniquePoints, ~, idx] = unique(centerP, 'rows', 'stable');
% Remove elements with same midpoint
duplicateFlags = false(size(centerP, 1), 1);
for i = 1:size(uniquePoints, 1)
    pointIndices = find(idx == i);
    if numel(pointIndices) > 1
        duplicateFlags(pointIndices(2:end)) = true;
    end
end
elements1(duplicateFlags,:) = [];
end

function centerP = lineCenterPosition(nodes,elements)
%% Calculate the coordinates of the beam midpoint
centerP = zeros(length(elements),3);
for i=1:length(elements)
    startp = nodes(elements(i,1),:);
    endp = nodes(elements(i,2),:);
    centerP(i,:) = mean([startp;endp]);
end
end

function [onZ, onX, onY] = checkPlanes(points,xmin,xmax,ymin,ymax,zmin,zmax)
%% Determine whether a point is on a specific plane
th = 1e-4;
% Determine whether the point is on the plane with z = zmax
onZ = abs(points(:,3) - zmax) < th;
% Determine whether the point is on the plane with x = xmin
onX = abs(points(:,1) - xmin) < th;
% Determine whether the point is on the plane with y = ymin
onY = abs(points(:,2) - ymin) < th;
end

function [xmin,xmax,ymin,ymax,zmin,zmax] = MaxMinCubic(points)
%% Calculate the bounding box of the cube truss cell
xmin = min(points(:,1));xmax = max(points(:,1));
ymin = min(points(:,2));ymax = max(points(:,2));
zmin = min(points(:,3));zmax = max(points(:,3));
end

function [newNodes, newElements] = removeIsolatedNodes(nodes, elements)
%% Remove isolated nodes
% Find all connected nodes
connectedNodes = unique(elements(:));
newNodes = nodes(connectedNodes, :);
% Update node number
newNodeIdx = zeros(size(nodes, 1), 1);
newNodeIdx(connectedNodes) = 1:length(connectedNodes);
newElements = newNodeIdx(elements);
end