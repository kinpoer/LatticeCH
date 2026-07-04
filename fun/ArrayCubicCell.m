function [newNodes,newElements] = ArrayCubicCell(nodes,elements,m,n,k,threshold)
% Multiscale Homogenization Code for Periodic Lattice Metamaterial
% Jin BO, Email:kinpoer@nuaa.edu.cn
%% Array 2D/3D truss cell and merge coincident boundary nodes.
% The returned arrayed cell is centered at the origin.  For example, if the
% original cubic cell is in [-0.5,0.5]^3, a 2-by-2-by-2 array is returned in
% [-1,1]^3.  The caller may then divide the coordinates by 2 to recover the
% original bounding-box size.

if nargin < 6 || isempty(threshold)
    threshold = autoGeometryTolerance(nodes);
end
if nargin < 5 || isempty(k)
    k = 1;
end

m = round(m);n = round(n);k = round(k);
if m < 1 || n < 1 || k < 1
    error('ArrayCubicCell:InvalidArraySize','m, n, and k must be positive integers.');
end

nodes = double(nodes);
elements = double(elements);
ndim = size(nodes,2);
if ndim ~= 2 && ndim ~= 3
    error('ArrayCubicCell:InvalidDimension','nodes must be an N-by-2 or N-by-3 array.');
end

cmax = max(nodes,[],1);
cmin = min(nodes,[],1);
dc = cmax-cmin;
if any(dc < threshold)
    error('ArrayCubicCell:DegenerateCell','The input cell has a near-zero bounding-box dimension.');
end

baseNodes = nodes-cmin;
newNodes = [];
newElements = [];
nodeNum = size(nodes,1);
cellId = 0;

if ndim == 3
    for ii = 1:m
        for jj = 1:n
            for kk = 1:k
                cellId = cellId+1;
                shift = [(ii-1)*dc(1),(jj-1)*dc(2),(kk-1)*dc(3)];
                newNodes = [newNodes;baseNodes+shift]; %#ok<AGROW>
                newElements = [newElements;elements+(cellId-1)*nodeNum]; %#ok<AGROW>
            end
        end
    end
else
    for ii = 1:m
        for jj = 1:n
            cellId = cellId+1;
            shift = [(ii-1)*dc(1),(jj-1)*dc(2)];
            newNodes = [newNodes;baseNodes+shift]; %#ok<AGROW>
            newElements = [newElements;elements+(cellId-1)*nodeNum]; %#ok<AGROW>
        end
    end
end

[newNodes,newElements] = MergeCloseNodesAndElements(newNodes,newElements,threshold);
newNodes = newNodes-0.5*(max(newNodes,[],1)+min(newNodes,[],1));
end

function tol = autoGeometryTolerance(nodes)
boxSize = norm(max(nodes,[],1)-min(nodes,[],1));
tol = max(1e-10,1e-8*max(1,boxSize));
end
