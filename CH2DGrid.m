% Multiscale Homogenization Code for Periodic Lattice Metamaterial
% Jin BO, Email:kinpoer@nuaa.edu.cn
%% Initialization
clc;clear;close all;addpath('fun','src','data');

%% Importing truss cells
load geo2D.mat %for 2D grid UCs 

%% Calculate the effective performance of each cell
for j=1:length(tempFiles)
    %% Truss Topology
    VOC = tempFiles{j,1};
    nodes = VOC.N;nodes(:,3) = 0;% To avoid errors, it is recommended to scale the dimensions to 20
    Plot3DCell(VOC.N,VOC.L,'k-');hold on;

    %% Primitive Cell
    [nodes,elements] = SimplyCubicStructure(nodes,VOC.L);
    [nodes,elements] = DivideBeam(nodes,elements, 5);

    %% PBC Matching
    cornerIndices = FindCubicCorner(nodes);
    edgeIndices = FindCubicEdge(nodes);
    if max(nodes(:,3))~=min(nodes(:,3))
        surfaceIndices = FindCubicSurface(nodes);
    else
        surfaceIndices = [];
    end
    innerIndices = setdiff([1:length(nodes)]',surfaceIndices);
    mapping = MatchCubicBoundary(nodes,cornerIndices,edgeIndices,surfaceIndices);

    %% Assigning structural properties
    sideLength = 20;
    r = ConstantRo2D(nodes,elements,0.05,sideLength^2);
    r(2) = r;%for rect, r(1) = H/2;r(2) = b/2;
    Es = 1;vs = 0.3;multi = 1;
    crossSection = 'Rect';% 'Circle'
    boxType = 'Cubic';
    elementType = 'Euler' ;%'Euler' 'Timoshenko'

    property = BeamProperty(nodes,elements,r,Es,vs,multi,boxType,crossSection);
    property.NewNode = property.multi*nodes;
    property.NewLine = elements;
    property.BoundMap = mapping;
    property.ElementType = elementType;

    %% Independent Node Matrix
    B0 = B0Matrix(property,mapping);

    %% Vector Displacement Matrix
    Ba = BaCubicMatrix(property,mapping);

    %% Truss Stiffness Matrix
    [Kuc,property] = OutputKuc2D(property);

    %% Three-periodic vector matrix
    Be = CubicBeMatrix(property.a1,property.a2,property.a3);
    [D0,infoD0] = SolveD0(Kuc,B0,Ba);
    Da = B0*D0+Ba;

    %% Unit cell Effective stiffness
    Ch = (1/property.S)*Be'*Da'*Kuc*Da*Be/2/r(2);% 2*r(2) for height or say thickness
    Ch(abs(Ch)<1e-20) = 0;
    Sh = inv(Ch);
    property.K = Ch;
    property.E11 = 1/Sh(1,1);
    property.E22 = 1/Sh(2,2);
    property.G12 = 1/Sh(3,3);
    property.v12 = -Sh(1,2)*property.E22;
    property.v21 = -Sh(2,1)*property.E11;

    %% Visualizing Effective Parameters
    % This part of the code comes from:
    % Mingqing Liao, Yong Liu, Nan Qu e.t. al. ElasticPOST: A Matlab Toolbox for Post-processing of Elastic Anisotropy with
    % Graphic User Interface. submitted to computer physics communication (2019)
    n = 200;
    flagSave = 0; %save or not
    flag = {'E', 'G'};
    ElasticPlot_2DM(Sh, n, flag, flagSave, 'Property');
    close all
    fprintf("The %d-th unit cell is running.......\n",j);
end




