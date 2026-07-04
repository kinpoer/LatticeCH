function d = PointToPolygonDistance(point,vertices)
% Multiscale Homogenization Code for Periodic Lattice Metamaterial
% Jin BO, Email:kinpoer@nuaa.edu.cn
%% Caculate ditance from point to line defined by two points
% Caculate normal vector of plane
v1 = vertices(2, :) - vertices(1, :);
v2 = vertices(3, :) - vertices(1, :);
normal = cross(v1, v2);
normal = normal / norm(normal); 
% Caculate distance from point to plane
pointOnPlane = vertices(1, :);
DPlane = abs(dot(normal, point - pointOnPlane));
% Check projection of point in polygon
projectedPoint = point - DPlane * normal;
if inpolygon3D(projectedPoint, vertices)
    d = DPlane;
else
    % Caculate distance from point to edges of polygon
    dEdges = inf;
    numVertices = size(vertices, 1);
    for i = 1:numVertices
        vStart = vertices(i, :);
        vEnd = vertices(mod(i, numVertices) + 1, :);
        dEdges = min(dEdges, PointToSegmentDistance(point, vStart, vEnd));
    end
    d = min(DPlane, dEdges);
end
end


function inside = inpolygon3D(point, vertices)
%% Determine if a point is inside a polygon (ignoring collinear vertices or non-planar polygons)
numVertices = size(vertices, 1);
angleSum = 0;
for i = 1:numVertices
    v1 = vertices(i, :) - point;
    v2 = vertices(mod(i, numVertices) + 1, :) - point;
    angleSum = angleSum + atan2(norm(cross(v1, v2)), dot(v1, v2));
end
inside = abs(angleSum - 2*pi) < 1e-5;
end
