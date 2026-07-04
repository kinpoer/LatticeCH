function newElements = BreakElements(nodes,elements,threshold)
% Multiscale Homogenization Code for Periodic Lattice Metamaterial
% Jin BO, Email:kinpoer@nuaa.edu.cn
%% Split elements when existing nodes lie on their interior.
% This function does not create new nodes.  It only uses the current node set
% and reconnects every element into smaller segments if intermediate nodes
% are found on the same straight line segment.

if nargin < 3 || isempty(threshold)
    threshold = autoGeometryTolerance(nodes);
end

nodes = double(nodes);
elements = double(elements);
if isempty(elements)
    newElements = zeros(0,2);
    return;
end

elements = round(elements(:,1:2));
nNode = size(nodes,1);
validRow = all(elements >= 1,2) & all(elements <= nNode,2) & elements(:,1) ~= elements(:,2);
elements = elements(validRow,:);
newElements = zeros(0,2);

for i = 1:size(elements,1)
    id1 = elements(i,1);
    id2 = elements(i,2);
    p1 = nodes(id1,:);
    p2 = nodes(id2,:);
    d = p2-p1;
    L = norm(d);
    if L <= threshold
        continue;
    end

    vec = nodes-p1;
    t = (vec*d')/(L^2);
    proj = p1+t*d;
    dist = sqrt(sum((nodes-proj).^2,2));

    % Nodes strictly inside the segment.  The tolerance is scaled by length
    % to avoid adding endpoints again.
    inside = find(t > threshold/L & t < 1-threshold/L & dist <= threshold);
    ids = [id1;inside(:);id2];
    ts = [0;t(inside);1];

    [ts,order] = sort(ts);
    ids = ids(order);

    % Collapse multiple nodes with almost identical projected positions.
    keep = true(size(ids));
    for j = 2:length(ids)
        if abs(ts(j)-ts(j-1))*L <= threshold
            keep(j) = false;
        end
    end
    ids = ids(keep);

    if length(ids) >= 2
        newElements = [newElements;[ids(1:end-1),ids(2:end)]]; %#ok<AGROW>
    end
end

if isempty(newElements)
    newElements = zeros(0,2);
else
    newElements = newElements(newElements(:,1) ~= newElements(:,2),:);
    newElements = unique(sort(newElements,2),'rows');
end
end

function tol = autoGeometryTolerance(nodes)
boxSize = norm(max(nodes,[],1)-min(nodes,[],1));
tol = max(1e-10,1e-8*max(1,boxSize));
end
