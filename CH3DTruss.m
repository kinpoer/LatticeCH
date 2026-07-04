% Multiscale Homogenization Code for Periodic Lattice Metamaterial
% Jin BO, Email:kinpoer@nuaa.edu.cn
%% Initialization
clc;clear;close all;addpath('fun','src','data');

%% Importing truss cells
load geo3D.mat %for truss UCs

%% Assigning structural properties
Es = 1;vs = 0.3;density = 0.005;
multi = 1;
crossSection = 'Circle';% 'Circle' or 'Tube'
boxType = 'Cubic';
elementType = 'Euler';%'Euler' or 'Timoshenko'
sideLength = 20;% to avoid errors, it is recommended to scale the dimensions to 20
% comcase = cell(length(tempFiles),1);

%% Calculate the effective performance of each cell
for j=1:length(ICell)% tempFiles

    tic
    VOC = tempFiles{j,1};
    nodes = VOC.N;elements = VOC.L;
    Plot3DCell(nodes, elements,'b-');title('Origin Cell');%ViewPer();

    %% Primitive Cell
    [nodes, elements] = SimplyCubicStructure(nodes, elements);
    Plot3DCell(nodes, elements,'b-');title('Primitive Cell');%ViewPer();
    r = ConstantRo(nodes, elements, density, sideLength^3);
    % [nodes, elements] = DivideBeam(nodes, elements, 5);
    property = BeamProperty(VOC.N, VOC.L, r, Es, vs, multi, boxType, crossSection);
    property.ElementType = elementType;
    property.NewNode = nodes;
    property.NewLine = elements;

    %% PBC Matching
    cornerIndices = FindCubicCorner(nodes);
    edgeIndices = FindCubicEdge(nodes);
    surfaceIndices = FindCubicSurface(nodes);
    innerIndices = setdiff([1:length(nodes)]',surfaceIndices);
    mapping = MatchCubicBoundary(nodes,cornerIndices,edgeIndices,surfaceIndices);
    property.BoundMap = mapping;

    %% Independent Node Matrix
    B0 = B0Matrix(property,mapping);% EQ.(8)

    %% Vector Displacement Matrix
    Ba = BaCubicMatrix(property,mapping);% EQ.(8)

    %% Truss Stiffness Matrix
    [Kuc,property] = OutputKuc(property);% EQ.(11)

    %% Three-periodic vector matrix
    Be = CubicBeMatrix(property.a1,property.a2,property.a3);% EQ.(21)
    [D0,infoD0] = SolveD0(Kuc,B0,Ba);% EQ.(14)
    Da = B0*D0+Ba;% EQ.(17)

    %% Unit cell Effective stiffness
    Ch = (1/property.V)*Be'*Da'*Kuc*Da*Be;% EQ.(23)
    % Ch = comcase2{1,1}*1000;
    % Ch(abs(Ch)<1e-10) = 0;
    Sh = pinv(Ch);
    % -cubic truss cells- %
    E = (Ch(1,1)^2+Ch(1,2)*Ch(1,1)-2*Ch(1,2)^2)/(Ch(1,1)+Ch(1,2));
    G = Ch(4,4);
    B = (Ch(1,1)+2*Ch(1,2))/3;
    mu = Ch(1,2)/(Ch(1,1)+Ch(1,2));
    E11 = 1/Sh(1,1);E22 = 1/Sh(2,2);E33 = 1/Sh(3,3);
    G12 = 1/Sh(4,4);G23 = 1/Sh(5,5);G13 = 1/Sh(6,6);
    v12 = -Sh(1,2)*E11;v32 = -Sh(2,3)*E33;v31 = -Sh(1,3)*E33;
    v21 = -Sh(1,2)*E22;v23 = -Sh(2,3)*E22;v13 = -Sh(1,3)*E11;
 
    %% Visualizing Effective Parameters
    % This part of the code comes from:
    % Mingqing Liao, Yong Liu, Nan Qu e.t. al. ElasticPOST: A Matlab Toolbox for Post-processing of Elastic Anisotropy with
    % Graphic User Interface. submitted to computer physics communication (2019)
    n = 300;
    flagSave = 0; %save or not
    % Bulk modulus(B),Young's modulus(E),Shear modulus(G),Poisson's ratio(v),Hardness(H)
    flag = {'E', 'G'};
    ElasticPlot_3D(pinv(Ch), n, flag, 'Max', flagSave, 'Property');
    % Ch = comcase2{1,1}*1000;

    toc
    close all
    fprintf("The %d-th unit cell is running.......\n",j);
end




