function edgeIndices = FindCubicEdge(points,th)
% Multiscale Homogenization Code for Periodic Lattice Metamaterial
% Jin BO, Email:kinpoer@nuaa.edu.cn
%% Find the edge points of cubic truss cell
if nargin < 2
    th = 1e-3; 
end
if max(points(:,3))~=min(points(:,3))
    maxCoord = max(points, [], 1);
    minCoord = min(points, [], 1);
    vertices = [minCoord(1), minCoord(2), minCoord(3);
        maxCoord(1), minCoord(2), minCoord(3);
        maxCoord(1), maxCoord(2), minCoord(3);
        minCoord(1), maxCoord(2), minCoord(3);
        minCoord(1), minCoord(2), maxCoord(3);
        maxCoord(1), minCoord(2), maxCoord(3);
        maxCoord(1), maxCoord(2), maxCoord(3);
        minCoord(1), maxCoord(2), maxCoord(3)];
    edges = [vertices(1,:), vertices(2,:);
        vertices(2,:), vertices(3,:);
        vertices(3,:), vertices(4,:);
        vertices(4,:), vertices(1,:);
        vertices(5,:), vertices(6,:);
        vertices(6,:), vertices(7,:);
        vertices(7,:), vertices(8,:);
        vertices(8,:), vertices(5,:);
        vertices(1,:), vertices(5,:);
        vertices(2,:), vertices(6,:);
        vertices(3,:), vertices(7,:);
        vertices(4,:), vertices(8,:)];
else
    maxCoord = max(points, [], 1);
    minCoord = min(points, [], 1);
    vertices = [minCoord(1), minCoord(2), 0;
        maxCoord(1), minCoord(2), 0;
        maxCoord(1), maxCoord(2), 0;
        minCoord(1), maxCoord(2), 0];
    edges = [vertices(1,:), vertices(2,:);
        vertices(2,:), vertices(3,:);
        vertices(3,:), vertices(4,:);
        vertices(4,:), vertices(1,:)];
end

id = zeros(size(points,1),1);
for i = 1:size(edges, 1)
    for j=1:size(points,1)
        distance = distancePointToLine(points(j,:),edges(i,1:3),edges(i,4:6));
        if distance<th
            id(j) = 1;
        end
    end
end
edgeIndices=find(id == 1);
end


function distance = distancePointToLine(point, linePoint1, linePoint2)
%% Calculate the distance from a point to a line defined by two points
lineVector = linePoint2 - linePoint1;
pointVector = point - linePoint1;
distance = norm(cross(lineVector, pointVector)) / norm(lineVector);
end

