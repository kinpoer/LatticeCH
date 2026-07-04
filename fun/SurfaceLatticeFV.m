function U = SurfaceLatticeFV(X,Y,Z,vn,fn)
% Multiscale Homogenization Code for Periodic Lattice Metamaterial
% Jin BO, Email:kinpoer@nuaa.edu.cn
%% Caculate volumetric distance field of surface lattices
U = ones(size(X))*10000;
for i = 1:numel(X)
    center = [X(i) Y(i) Z(i)];
    for j = 1:length(fn)
        distance = PointToPolygonDistance(center,vn(fn(j,:),:));
        % U(i) = min(U(i),distance);
        k = 4; % the smoothing parameter
        U(i) = -log(exp(-k*U(i)) + exp(-k*distance))/k;
    end
end
end

