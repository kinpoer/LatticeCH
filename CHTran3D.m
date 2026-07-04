% Multiscale Homogenization Code for Periodic Lattice Metamaterial
% Jin BO, Email:kinpoer@nuaa.edu.cn
%% Initialization
clc;clear;close all;addpath('fun','src','data');

%% Importing truss cells
load geo3D.mat %for truss UCs 

%% Assigning structural properties
Es = 1;vs = 0.3;density = 0.01;
multi = 1;
crossSection = 'Circle';
boxType = 'Cubic';
elementType = 'Euler';%'Euler' 'Timoshenko' 
sideLength = 20;

%% Calculate the effective performance of each cell
id = 8;
VOC = tempFiles{id,1};
nodes = VOC.N;
elements = VOC.L;
Plot3DCell(nodes, elements,'b-');GenerateCubic(20,'on');title('Before transformation');axis on

%% Primitive Cell
[nodes, elements] = SimplyCubicStructure(nodes, elements);
Plot3DCell(nodes, elements,'b-');GenerateCubic(20,'on');
r = ConstantRo(nodes, elements, density, sideLength^3);
% [nodes, elements] = DivideBeam(nodes, elements, 5);
property = BeamProperty(VOC.N, VOC.L, r, Es, vs, multi, boxType, crossSection);
property.NewNode = nodes;
property.NewLine = elements;

%% PBC Matching (Geometric transformations do not change PBC)
cornerIndices = FindCubicCorner(nodes);
edgeIndices = FindCubicEdge(nodes);
surfaceIndices = FindCubicSurface(nodes);
innerIndices = setdiff([1:length(nodes)]',surfaceIndices);
mapping = MatchCubicBoundary(nodes,cornerIndices,edgeIndices,surfaceIndices);

%% Independent Node Matrix
B0 = B0Matrix(property,mapping);

%% Vector Displacement Matrix
Ba = BaCubicMatrix(property,mapping);

%% Geometric transformation
A = GenerateCubic(20,'off');
load TransCubic.mat
% load StretchCubic.mat

i = 1;
tranT = angleT;% Form1:angleT; Form2:lenT; Form3:rotateT
%----------------------------------------------------------%
% 1: Update nodes for Form3 rotateT
% NewNode = RotatePointsAroundZ(nodes,45);%nodes
% for j=1:length(NewNode)
%     NewNode(j,1) = NewNode(j,1)+tranT{i,2}.*NewNode(j,3)./2./tan(pi/2-tranT{i,1});
%     NewNode(j,3) = NewNode(j,3)*tranT{i,2};
% end
% TranNode = RotatePointsAroundZ(VOC.L,45);%nodes
% for j=1:length(TranNode)
%     TranNode(j,1) = TranNode(j,1)+tranT{i,2}.*TranNode(j,3)./2./tan(pi/2-tranT{i,1});
%     TranNode(j,3) = TranNode(j,3)*tranT{i,2};
% end
%----------------------------------------------------------%
% 2 and 3: Update nodes for Form1 angleT or Form2 lenT
B = tranT{i,3};
NewNode = CubicMap(A, B, nodes);
TranNode = CubicMap(A, B, VOC.N);
%----------------------------------------------------------%
r = ConstantRo(NewNode,elements,0.01,tranT{i,5});
property = BeamProperty(NewNode, VOC.L, r, Es, vs, multi, boxType, crossSection);
property.SideLength = sideLength;
property.ElementType = elementType;
property.BoundMap = mapping;
property.NewNode = NewNode;
property.NewLine = elements;
Plot3DCell(NewNode,elements,'b-');PlotCubicPrism(tranT{i,3});
Plot3DCell(TranNode, VOC.L,'b-');title('After transformation');axis on;PlotCubicPrism(tranT{i,3});

property.a1 = tranT{i,4}(1,:);
property.a2 = tranT{i,4}(2,:);
property.a3 = tranT{i,4}(3,:);
property.V = tranT{i,5};

%% Truss Stiffness Matrix
[Kuc,property] = OutputKuc(property);

%% Three-periodic vector matrix
Be = CubicBeMatrix(property.a1,property.a2,property.a3);
[D0,infoD0] = SolveD0(Kuc,B0,Ba);% EQ.(9)
Da = B0*D0+Ba;% d = Da * delta_e

%% Unit cell Effective stiffness
Ch = (1/property.V)*Be'*Da'*Kuc*Da*Be;
% Ch(abs(Kana)<1e-10) = 0;
Sh = inv(Ch);

%% Visualizing Effective Parameters
% This part of the code comes from:
% Mingqing Liao, Yong Liu, Nan Qu e.t. al. ElasticPOST: A Matlab Toolbox for Post-processing of Elastic Anisotropy with
% Graphic User Interface. submitted to computer physics communication (2019)
n = 150;
% flagSave = 0; %save or not
% flag = {'E', 'G', 'B', 'v'};
% ElasticPlot_3D(inv(Ch), n, flag, 'Max', flagSave, 'Property');

property.K = Ch;
property.E11 = 1/Sh(1,1);
property.E22 = 1/Sh(2,2);
property.E33 = 1/Sh(3,3);
property.G12 = 1/Sh(4,4);
property.G23 = 1/Sh(5,5);
property.G13 = 1/Sh(6,6);

