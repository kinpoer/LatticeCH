% Multiscale Homogenization Code for Periodic Lattice Metamaterial
% Jin BO, Email:kinpoer@nuaa.edu.cn
%% Initialization
clc;clear;close all;addpath('fun','src','data');

%% Assigning structural properties
Es = 1e6;vs = 0.3;
thick = 0.9;
multi = 0.8;
boxType = 'Cubic';
sideLength = 20;% to avoid errors, it is recommended to scale the dimensions to 20

%% Calculate the effective performance of each cell
tic
load Shell-N2.mat % for 3D surface UC
nodes(:,3) = nodes(:,3)-10;
Plot3DMesh(nodes,elements,1)

%% Primitive Cell
property = ShellProperty(nodes, elements, Es, vs, multi, boxType);
property.SideLength = sideLength;
property.NewNode = multi*nodes;
property.NewLine = elements;

%% PBC Matching (Users need to prepare periodic mesh)
cornerIndices = FindCubicCorner(nodes);
edgeIndices = FindCubicEdge(nodes);
hold on;scatter3(nodes(edgeIndices,1),nodes(edgeIndices,2),nodes(edgeIndices,3),'k','filled','SizeData',30);
surfaceIndices = FindCubicSurface(nodes);
hold on;scatter3(nodes(surfaceIndices,1),nodes(surfaceIndices,2),nodes(surfaceIndices,3),'b','filled','SizeData',30);
innerIndices = setdiff([1:length(nodes)]',surfaceIndices);
mapping = MatchCubicBoundary(nodes,cornerIndices,edgeIndices,surfaceIndices);
property.BoundMap = mapping;

%% Independent Node Matrix
B0 = B0Matrix(property,mapping);

%% Vector Displacement Matrix
Ba = BaCubicMatrix(property,mapping);

%% Truss Stiffness Matrix
Kuc = Q3Shell3D(nodes,elements,Es,vs,thick);

%% Three-periodic vector matrix
Be = CubicBeMatrix(property.a1,property.a2,property.a3);
[D0,infoD0] = SolveD0(Kuc,B0,Ba);
Da = B0*D0+Ba;% d = Da * delta_e

%% Unit cell Effective stiffness
Ch = (1/property.V)*Be'*Da'*Kuc*Da*Be;
Ch(abs(Ch)<1e-10) = 0;
Ch = (Ch+Ch')/2;
Sh = inv(Ch);
% Cubic truss cells
E = (Ch(1,1)^2+Ch(1,2)*Ch(1,1)-2*Ch(1,2)^2)/(Ch(1,1)+Ch(1,2));
G = Ch(4,4);
B = (Ch(1,1)+2*Ch(1,2))/3;
mu = Ch(1,2)/(Ch(1,1)+Ch(1,2));
% General truss cells
E11 = 1/Sh(1,1);E22 = 1/Sh(2,2);E33 = 1/Sh(3,3);
G12 = 1/Sh(4,4);G23 = 1/Sh(5,5);G13 = 1/Sh(6,6);
v12 = -Sh(1,2)*E11;v32 = -Sh(2,3)*E33;v31 = -Sh(1,3)*E33;
v21 = -Sh(1,2)*E22;v23 = -Sh(2,3)*E22;v13 = -Sh(1,3)*E11;

%% Visualizing Effective Parameters
% This part of the code comes from:
% Mingqing Liao, Yong Liu, Nan Qu e.t. al. ElasticPOST: A Matlab Toolbox for Post-processing of Elastic Anisotropy with
% Graphic User Interface. submitted to computer physics communication (2019)
n = 200;
flagSave = 0; %save or not
% Bulk modulus(B),Young's modulus(E),Shear modulus(G),Poisson's ratio(v),Hardness(H)
flag = {'E', 'G'};
ElasticPlot_3D(inv(Ch), n, flag, 'Max', flagSave, 'Property');
toc
