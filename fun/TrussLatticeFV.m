function U = TrussLatticeFV(X,Y,Z,nodes,struts,rstrut,rnode)
% Multiscale Homogenization Code for Periodic Lattice Metamaterial
% Jin BO, Email:kinpoer@nuaa.edu.cn
%% caculate volumetric distance field of truss or tube lattices
U = ones(size(X))*10000;
% Get caculate distance from line to grid points
for i = 1:numel(X)
    center = [X(i) Y(i) Z(i)]; % voxel center position
    for j = 1:size(struts,1)
        startP = nodes(struts(j,1),:);  % start node coordinate
        endP = nodes(struts(j,2),:);    % end node coordinate
        distance2 = min(norm(center - startP),norm(center - endP));
        distance1 = PointToSegmentDistance(center, startP, endP);
        tempDist = min(distance1-rstrut(j),distance2-rnode);
        k = 4; % the smoothing parameter
        U(i) = -log(exp(-k*U(i)) + exp(-k*tempDist))/k;
    end
end
end
