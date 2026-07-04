function [Vn,Fn] = Symmetrize3DMesh(Vn,Fn,tolerance)
% Multiscale Homogenization Code for Periodic Lattice Metamaterial
% Jin BO, Email:kinpoer@nuaa.edu.cn
%% Symmetry operations of a octants in cube
pts = Vn;
tri = Fn;
% Generate combanation (±1, ±1, ±1)
[signX, signY, signZ] = ndgrid([1 -1], [1 -1], [1 -1]);
signCombinations = [signX(:), signY(:), signZ(:)]; % 8x3

% Extend vertices
numOriginalPts = size(pts, 1);
allPts = zeros(numOriginalPts * 8, 3);
for i = 1:8
    signs = signCombinations(i, :);
    idxRange = (i-1)*numOriginalPts + 1 : i*numOriginalPts;
    allPts(idxRange, :) = pts .* signs;
end

% remove dup vertices
[uniquePts, ~, ic] = uniquetol(allPts, tolerance, 'ByRows', true);

% Extend tri mesh
numOriginalTri = size(tri, 1);
allTri = zeros(numOriginalTri * 8, 3);
for i = 1:8
    offset = (i-1)*numOriginalPts;
    idxRange = (i-1)*numOriginalTri + 1 : i*numOriginalTri;
    allTri(idxRange, :) = tri + offset;
end

% Map new index
allTri = ic(allTri);

% Create triangulation
TRFull = triangulation(allTri, uniquePts);
Vn = TRFull.Points;
Fn = TRFull.ConnectivityList;
end