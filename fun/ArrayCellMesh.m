function [nodes,elements] = ArrayCellMesh(nodes,elements,m,n,k)
% Multiscale Homogenization Code for Periodic Lattice Metamaterial
% Jin BO, Email:kinpoer@nuaa.edu.cn
%% Array mesh of UC and update mesh
if size(nodes,2)==3
    cmax = max(nodes);cmin = min(nodes);dc = cmax-cmin;
    dx = dc(1);dy = dc(2);dz = dc(3);
    % Left-Bottom corner
    center = cmin;
end
% Array center points first
centers = [];
if size(nodes,2)==3
for i=1:m
    for j=1:n
       for z=1:k
           tempCenter = center + [(i-1)*dx (j-1)*dy (z-1)*dz];
           centers = [centers; tempCenter];
       end
    end
end
end
% Array nodes and elements next
allcell = m*n*k;
allNodeNum = allcell*length(nodes);
allEleNum = allcell*length(elements);
newNodes = zeros(allNodeNum,3);
newEles = zeros(allEleNum,3);
newNodes(1:length(nodes),:) = nodes+center;
newEles(1:length(elements),:) = elements;

for i=2:length(centers)
    startId = (i-1)*length(nodes)+1;
    endId = i*length(nodes);
    tempNode = nodes+centers(i,:);
    newNodes(startId:endId,:) = tempNode;
end
lennode = length(nodes);
for i=2:length(centers)
    startId = (i-1)*length(elements)+1;
    endId = i*length(elements);
    tempEle = elements(:,:)+ones(size(elements,1),size(elements,2))*(i-1)*lennode;
    newEles(startId:endId,:) = tempEle;
end
nodes = newNodes;elements = newEles;
end