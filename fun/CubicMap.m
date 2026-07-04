function mappedPoint = CubicMap(A, B, P)
% Multiscale Homogenization Code for Periodic Lattice Metamaterial
% Jin BO, Email:kinpoer@nuaa.edu.cn
%% Mapping points from one prism to another
tempA = A(1:4,1:2);tempB = B(1:4,1:2);
tempP = P(:,1:2);
mappedPoint=zeros(length(tempP),3);
for i=1:size(tempP,1)
mappedPoint(i,1:2) = MapCubicPointToParallelogramWithCentroid(tempP(i,:), tempA, tempB);
end
maxA = max(A(:,3));maxB = max(B(:,3));
ratio = maxB/maxA;
mappedPoint(:,3) = ratio.*P(:,3);
% Visualize the results
% figure; hold on;
% plot3(P(:,1),P(:,2),P(:,3), 'bo', 'MarkerSize', 10, 'LineWidth', 2);
% plot3(mappedPoint(:,1), mappedPoint(:,2),mappedPoint(:,3), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
% view(3)
% axis equal;grid on;hold off;
end

function mappedPoint = MapCubicPointToParallelogramWithCentroid(point, Vertices, targetVertices)
%% Mapping point from one prism to another
% Calculate the centroid of origin cubic
centroid1 = mean(Vertices);
centroid1(abs(centroid1)<1e-7) = 0;
centroid2 = mean(targetVertices);
centroid2(abs(centroid2)<1e-7) = 0;
% Find the triangle where points is located 
triIndex = findTriangleWithMappedP(point, Vertices, centroid1);
% Determine the transformed position of the points by transforming the tetrahedron
P = [centroid1; Vertices(triIndex, :)];
Q = [centroid2; targetVertices(triIndex, :)];
mappedPoint = mapPointToAnotherTriangle(point, P, Q);
end

function triIndex = findTriangleWithMappedP(point, Vertices, centroid)
%% Find the triangle index where the point lies
for i = 1:size(Vertices,1)
    triIndex = [i, mod(i, 4) + 1];
    if isPInTriangle(point, [centroid; Vertices(triIndex, :)])
        return;
    end
end
error('Point is not inside the hexagon.');
end

function isIn = isPInTriangle(point, triangle, th)
%% Determine whether a point is within the tetrahedron
if nargin < 3
    th = 1e-5; 
end
x1 = triangle(1, 1);y1 = triangle(1, 2);
x2 = triangle(2, 1);y2 = triangle(2, 2);
x3 = triangle(3, 1);y3 = triangle(3, 2);
xp = point(1);yp = point(2);
originalArea = abs(det([x1,y1,1;x2,y2,1;x3,y3,1])/2);
area1 = abs(det([xp,yp,1;x2,y2,1;x3,y3,1])/2);
area2 = abs(det([x1,y1,1;xp,yp,1;x3,y3,1])/2);
area3 = abs(det([x1,y1,1;x2,y2,1;xp,yp,1])/2);
totalArea = area1 + area2 + area3;
isIn = abs(totalArea-originalArea)<th;
end

function weights = computeBarycentricWeights(p, tri)
%% Calculate the center of gravity weight
xp = p(1);yp = p(2);
x1 = tri(1,1);y1 = tri(1,2);
x2 = tri(2,1);y2 = tri(2,2);
x3 = tri(3,1);y3 = tri(3,2);
denominator = ((y2-y3)*(x1-x3)+(x3-x2)*(y1-y3));
w1 = ((y2-y3)*(xp-x3)+(x3-x2)*(yp-y3))/denominator;
w2 = ((y3-y1)*(xp-x3)+(x1-x3)*(yp-y3))/denominator;
w3 = 1-w1-w2;
weights = max(0, min([w1, w2, w3], 1));
end

function newPoint = mapPointToAnotherTriangle(point, tri1, tri2)
%% Mapping points to another triangle
weights = computeBarycentricWeights(point, tri1);
newPoint = weights(1)*tri2(1,:)+weights(2)*tri2(2,:)+weights(3)*tri2(3,:);
end

