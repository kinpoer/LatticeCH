function [newNodes,newElements,info] = CheckCross(nodes,elements,threshold)
% Multiscale Homogenization Code for Periodic Lattice Metamaterial
% Jin BO, Email:kinpoer@nuaa.edu.cn
%% Add intersection nodes for crossed truss members and split elements.
% Compared with the original implementation, this version first collects all
% intersections and then splits every affected element once.  This avoids
% duplicated, missing, or overlapping segments when one member has multiple
% intersections.

if nargin < 3 || isempty(threshold)
    threshold = autoGeometryTolerance(nodes);
end

nodes = double(nodes);
elements = double(elements);
ndim = size(nodes,2);
if ndim ~= 2 && ndim ~= 3
    error('CheckCross:InvalidDimension','nodes must be an N-by-2 or N-by-3 array.');
end

if isempty(elements)
    newNodes = nodes;
    newElements = zeros(0,2);
    info.numIntersections = 0;
    return;
end

% Work in 3D internally.  2D inputs are padded with z=0.
if ndim == 2
    workNodes = [nodes,zeros(size(nodes,1),1)];
else
    workNodes = nodes;
end

elements = round(elements(:,1:2));
nNode = size(workNodes,1);
validRow = all(elements >= 1,2) & all(elements <= nNode,2) & elements(:,1) ~= elements(:,2);
elements = elements(validRow,:);
elements = unique(sort(elements,2),'rows');
nEle = size(elements,1);

newWorkNodes = workNodes;
splitT = cell(nEle,1);
splitId = cell(nEle,1);
for i = 1:nEle
    splitT{i} = [0;1];
    splitId{i} = elements(i,:)';
end

numIntersections = 0;
for i = 1:nEle
    p1 = workNodes(elements(i,1),:);
    p2 = workNodes(elements(i,2),:);
    L1 = norm(p2-p1);
    if L1 <= threshold
        continue;
    end

    for j = i+1:nEle
        % Already connected members do not require a new intersection node.
        if any(elements(i,1) == elements(j,:)) || any(elements(i,2) == elements(j,:))
            continue;
        end

        q1 = workNodes(elements(j,1),:);
        q2 = workNodes(elements(j,2),:);
        L2 = norm(q2-q1);
        if L2 <= threshold
            continue;
        end

        [s,t,c1,c2,dist] = closestPointsOnSegmentsParam(p1,p2,q1,q2);
        if dist > threshold
            continue;
        end

        % Ignore endpoint contacts; they are handled by node merging or by
        % BreakElements if an endpoint lies on another member.
        endpointTol1 = threshold/L1;
        endpointTol2 = threshold/L2;
        if s <= endpointTol1 || s >= 1-endpointTol1 || t <= endpointTol2 || t >= 1-endpointTol2
            continue;
        end

        newPoint = 0.5*(c1+c2);
        newNodeId = findOrAddNode(newPoint);

        splitT{i}(end+1,1) = s; %#ok<AGROW>
        splitId{i}(end+1,1) = newNodeId; %#ok<AGROW>
        splitT{j}(end+1,1) = t; %#ok<AGROW>
        splitId{j}(end+1,1) = newNodeId; %#ok<AGROW>
        numIntersections = numIntersections+1;
    end
end

newElements = zeros(0,2);
for i = 1:nEle
    ids = splitId{i};
    ts = splitT{i};
    [ts,order] = sort(ts);
    ids = ids(order);

    p1 = newWorkNodes(elements(i,1),:);
    p2 = newWorkNodes(elements(i,2),:);
    L = norm(p2-p1);
    if L <= threshold
        continue;
    end

    % Remove duplicated split positions.
    keep = true(size(ids));
    for k = 2:length(ids)
        if abs(ts(k)-ts(k-1))*L <= threshold || ids(k) == ids(k-1)
            keep(k) = false;
        end
    end
    ids = ids(keep);

    if length(ids) >= 2
        newElements = [newElements;[ids(1:end-1),ids(2:end)]]; %#ok<AGROW>
    end
end

if ndim == 2
    newNodes = newWorkNodes(:,1:2);
else
    newNodes = newWorkNodes;
end

[newNodes,newElements] = MergeCloseNodesAndElements(newNodes,newElements,threshold);
newElements = BreakElements(newNodes,newElements,threshold);
[newNodes,newElements] = MergeCloseNodesAndElements(newNodes,newElements,threshold);

info.numIntersections = numIntersections;

    function nodeId = findOrAddNode(point)
        distToNodes = sqrt(sum((newWorkNodes-point).^2,2));
        [minDist,minId] = min(distToNodes);
        if minDist <= threshold
            nodeId = minId;
        else
            newWorkNodes(end+1,:) = point; %#ok<AGROW>
            nodeId = size(newWorkNodes,1);
        end
    end
end

function [s,t,c1,c2,dist] = closestPointsOnSegmentsParam(p1,p2,q1,q2)
% Robust closest points between two 3D segments.
d1 = p2-p1;
d2 = q2-q1;
r = p1-q1;
a = dot(d1,d1);
e = dot(d2,d2);
f = dot(d2,r);
smallNumber = 1e-14;

if a <= smallNumber && e <= smallNumber
    s = 0;t = 0;c1 = p1;c2 = q1;dist = norm(c1-c2);return;
end
if a <= smallNumber
    s = 0;
    t = clamp(f/e,0,1);
else
    c = dot(d1,r);
    if e <= smallNumber
        t = 0;
        s = clamp(-c/a,0,1);
    else
        b = dot(d1,d2);
        denom = a*e-b*b;
        if abs(denom) > smallNumber
            s = clamp((b*f-c*e)/denom,0,1);
        else
            % Nearly parallel.  Choosing s=0 is sufficient for the closest
            % distance test; endpoint/overlap splitting is handled elsewhere.
            s = 0;
        end
        t = (b*s+f)/e;
        if t < 0
            t = 0;
            s = clamp(-c/a,0,1);
        elseif t > 1
            t = 1;
            s = clamp((b-c)/a,0,1);
        end
    end
end

c1 = p1+s*d1;
c2 = q1+t*d2;
dist = norm(c1-c2);
end

function x = clamp(x,a,b)
x = max(a,min(b,x));
end

function tol = autoGeometryTolerance(nodes)
boxSize = norm(max(nodes,[],1)-min(nodes,[],1));
tol = max(1e-10,1e-8*max(1,boxSize));
end
