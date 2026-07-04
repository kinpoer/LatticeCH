function d = PointToSegmentDistance(point, vStart, vEnd)
% Multiscale Homogenization Code for Periodic Lattice Metamaterial
% Jin BO, Email:kinpoer@nuaa.edu.cn
%% Caculate ditance from point to line defined by two points
lineVec = vEnd - vStart;
pointVec = point - vStart;
lineLen = dot(lineVec, lineVec);
if lineLen == 0
    d = norm(pointVec); % vStart == vEnd
    return;
end
t = dot(pointVec, lineVec) / lineLen;
t = max(0, min(1, t)); 
nearestPoint = vStart + t * lineVec;
d = norm(point - nearestPoint);
end
