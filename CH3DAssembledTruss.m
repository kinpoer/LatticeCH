% Multiscale Homogenization Code for Periodic Lattice Metamaterial
% Jin BO, Email:kinpoer@nuaa.edu.cn
%% Initialization
% clc;clear;close all;addpath('fun','src','data');

%% Importing truss cells
load AssembledUC.mat %for truss UCs

%% Assigning structural properties
Es = 1e3;vs = 0.3;density = 0.01;
multi = 1;
crossSection = 'Circle';% 'Circle' or 'Tube'
boxType = 'Cubic';
elementType = 'Euler';%'Euler' or 'Timoshenko'
sideLength = 20;% to avoid errors, it is recommended to scale the dimensions to 20

%% Calculate the effective performance of each cell
for j=1:length(ICell)

    tic
    VOC = ICell{j,1};
    nodes = VOC.N;elements = VOC.L;
    Plot3DCell(nodes, elements,'b-');title('Origin Cell');%ViewPer();

    %% Primitive Cell
    [nodes, elements] = SimplyCubicStructure(nodes, elements);
    % Plot3DCell(nodes, elements,'b-');title('Primitive Cell');%ViewPer();
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
    % ElasticPlot_3D(pinv(Ch), n, flag, 'Max', flagSave, 'Property');

    toc
    close all
    fprintf("The %d-th unit cell is running.......\n",j);
end

ro = 0.01;Es = 1e3;vs = 0.3;Gs = Es/2/(1+vs);Ks = Es/3/(1-2*vs);
%% HSU
Ehsu = 2*ro*(5*vs-7)/(13*ro+12*vs-2*ro*vs-15*ro*vs^2+15*vs^2-27)*Es;
Khsu = (Ks+(1-ro)/(-1/Ks+ro*(Ks+(4/3)*Gs)^(-1)))*Es;
Ghsu = (Gs+(1-ro)/((-1/Gs+2*ro*(Ks+2*Gs)/(5*Gs*(Ks+4/3*Gs)))))*Es;
% Ehsu = 2*ro*Es*(7-5*vs)/(3*(1-vs)*(9+5*vs));
% Ghsu = ro*Es*(7-5*vs)/30/(1-vs)/(1+vs);
% Khsu = 2*ro*Es/9/(1-vs);
vhsu = (3*Khsu-2*Ghsu)/2/(3*Khsu+Ghsu);

%% Voigt
Evoigt = Es*ro;
Kvoigt = Ks*ro;
Gvoigt = Gs*ro;

comcase(:,1:3)=comcase(:,1:3)/Evoigt;comcase(:,4:6)=comcase(:,4:6)/Gvoigt;

figure;
for i=1:2744
scatter3(max(comcase(i,1:3)),max(comcase(i,4:6)),min(comcase(i,7:9)),'ro');hold on
end

figure;
for i=1:2744
scatter3((comcase(i,1)),(comcase(i,2)),(comcase(i,3)),'ro');hold on
end

figure;
for i=1:2744
scatter3((comcase(i,4)),max(comcase(i,5)),max(comcase(i,6)),'ro');hold on
end