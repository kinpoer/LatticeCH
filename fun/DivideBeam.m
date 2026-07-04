function [newNodes, newElements] = DivideBeam(nodes, elements, numDivisions)
% Multiscale Homogenization Code for Periodic Lattice Metamaterial
% Jin BO, Email:kinpoer@nuaa.edu.cn
%% Refined Beam Elements
newNodes = nodes;
newElements = [];
currentNodeCount = size(nodes, 1);
for i = 1:size(elements, 1)
    % Get the start and end node index of the current beam
    startNodeIdx = elements(i, 1);
    endNodeIdx = elements(i, 2);
    % Get the start and end node coordinate of the current beam
    startNode = nodes(startNodeIdx, :);
    endNode = nodes(endNodeIdx, :);
    % Generates new nodes between the start and end points on beam
    newSegmentNodes = LinspaceNode(startNode, endNode, numDivisions + 1);
    newNodes = [newNodes; newSegmentNodes(2:end-1, :)];
    % Update the new elements
    newElementIdx = currentNodeCount + (1:numDivisions - 1)';
    newElementIdx = [startNodeIdx; newElementIdx; endNodeIdx];
    for j = 1:numDivisions
        newElements = [newElements; newElementIdx(j), newElementIdx(j+1)];
    end
    % Update current node count
    currentNodeCount = currentNodeCount + numDivisions - 1;
end
end

function newSegmentNodes = LinspaceNode(startNode, endNode, numPoints)
%% Generates equally spaced nodes between the start and end points
if size(startNode,2)==3
    newSegmentNodes = zeros(numPoints, 3);
    for i = 1:3
        newSegmentNodes(:, i) = linspace(startNode(i), endNode(i), numPoints);
    end
else
    newSegmentNodes = zeros(numPoints, 2);
    for i = 1:2
        newSegmentNodes(:, i) = linspace(startNode(i), endNode(i), numPoints);
    end
end
end

