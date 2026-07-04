function [Kuc,property] = OutputKuc2D(property)
% Multiscale Homogenization Code for Periodic Lattice Metamaterial
% Jin BO, Email:kinpoer@nuaa.edu.cn
%% Calculate the cell stiffness matrix
Ver = property.NewNode;
line = property.NewLine;
no = [[1:length(Ver)]',Ver];
ele = [[1:length(line)]',line];
% Assign parameters to each elements
prop(1) = property.Es;
prop(2) = property.v;
prop(3) = 1;% mass density (mass per unit volume)
switch property.SectionType
    case 'Circle'
        prop(4) = property.D;% height of beam cross-section
        prop(5) = 0;% width of beam cross-section
        optSection = 2;
    case 'Rect'
        prop(4) = property.H;% height of beam cross-section
        prop(5) = property.b;% width of beam cross-section
        optSection = 1;
        property.D = prop(5);
end
prop(6) = property.SectionArea;% area of beam cross-section
prop(7) = property.Iy;% 2nd moment of inertia of cross-section about axis y
prop(8) = property.Iz;% 2nd moment of inertia of cross-section about axis z
prop(9) = property.J;
prop(10) = property.Gs;
noNode = size(no,1);
noNel = 2;% number of nodes per element
noDof = 3;% number of dofs per node
sysDof = noNode*noDof;% total system dofs
k = zeros(noNel*noDof,noNel*noDof);% element stiffness matrix
Kuc = zeros(sysDof,sysDof);% initialization of system stiffness matrix

Sreal = 0;
for iel=1:length(ele) 
    nd(1) = ele(iel,2);
    nd(2) = ele(iel,3); 
    x(1) = no(nd(1),2); y(1) = no(nd(1),3); z(1) = no(nd(1),4);
    x(2) = no(nd(2),2); y(2) = no(nd(2),3); z(2) = no(nd(2),4);
    L = sqrt((x(2)-x(1))^2+(y(2)-y(1))^2+(z(2)-z(1))^2);
    Sreal = Sreal + property.D*L;
    beta = atan2(y(2)-y(1), x(2)-x(1));
    switch  property.ElementType
        case 'Euler'
            % compute element matrix for Euler-Bernoulli beam
            k = EulerBeam2D(prop,L,beta,optSection);
        case 'Timoshenko'
            % compute element matrix for Timoshenko beam
            k = TimoshenkoBeam2D(prop,L,beta,optSection);
    end
    index = ElementDof(nd,noNel,noDof);% extract system dofs associated with element
    Kuc = MatrixAssembly2D(Kuc,k,index);% assemble system stiffness matrix
end
property.Kuc = Kuc;
property.RelativeRoh = Sreal/property.S;
end
