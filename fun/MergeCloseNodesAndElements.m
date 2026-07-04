function [newNodes,uniqueElements,map] = MergeCloseNodesAndElements(nodes,elements,th)
% Multiscale Homogenization Code for Periodic Lattice Metamaterial
% Jin BO, Email:kinpoer@nuaa.edu.cn
%% Merge close nodes and remove duplicated/zero-length elements.
% map(oldNodeId) gives the new node ID after merging.

if nargin < 3 || isempty(th)
    th = autoGeometryTolerance(nodes);
end
if isempty(nodes)
    newNodes = nodes;
    uniqueElements = elements;
    map = [];
    return;
end

nodes = double(nodes);
elements = double(elements);
nNode = size(nodes,1);
map = zeros(nNode,1);
newNodes = zeros(0,size(nodes,2));

% A simple tolerance-based clustering.  The cluster center is updated as the
% mean of all assigned nodes.  This is robust for exact duplicated nodes and
% small numerical noise at periodic/assembled interfaces.
for i = 1:nNode
    if map(i) ~= 0
        continue;
    end
    dist = sqrt(sum((nodes-nodes(i,:)).^2,2));
    closeIds = find(dist <= th & map == 0);
    newNodes(end+1,:) = mean(nodes(closeIds,:),1); %#ok<AGROW>
    map(closeIds) = size(newNodes,1);
end

if isempty(elements)
    uniqueElements = zeros(0,2);
    return;
end

% Remove invalid rows before remapping.
elements = elements(:,1:2);
validRow = all(isfinite(elements),2) & all(elements >= 1,2) & all(elements <= nNode,2);
elements = elements(validRow,:);
elements = round(elements);
if isempty(elements)
    uniqueElements = zeros(0,2);
    return;
end

uniqueElements = [map(elements(:,1)),map(elements(:,2))];
uniqueElements = uniqueElements(uniqueElements(:,1) ~= uniqueElements(:,2),:);
if isempty(uniqueElements)
    uniqueElements = zeros(0,2);
    return;
end

% For truss homogenization the element orientation is not used to distinguish
% two physical members.  Sorting removes duplicate connections robustly.
uniqueElements = unique(sort(uniqueElements,2),'rows');
end

function tol = autoGeometryTolerance(nodes)
boxSize = norm(max(nodes,[],1)-min(nodes,[],1));
tol = max(1e-10,1e-8*max(1,boxSize));
end
