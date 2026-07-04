function [Kuc,property] = OutputKuc(property)
% Multiscale Homogenization Code for Periodic Lattice Metamaterial
% Jin BO, Email:kinpoer@nuaa.edu.cn
%% Calculate the cell stiffness matrix
node = property.NewNode;
line = property.NewLine;
no = [[1:length(node)]',node];
ele = [[1:length(line)]',line];
% Assign parameters to each elements
for i=1:length(no)
    Node{i}.x = no(i,2);Node{i}.y = no(i,3);Node{i}.z = no(i,4);
    Node{i}.DOFx = (i-1)*6+1;
    Node{i}.DOFy = (i-1)*6+2;
    Node{i}.DOFz = (i-1)*6+3;
    Node{i}.DOFthetax = (i-1)*6+4;
    Node{i}.DOFthetay = (i-1)*6+5;
    Node{i}.DOFthetaz = (i-1)*6+6;
end
trussDirection = zeros(length(ele),3);
for i=1:length(ele)
    Element{i}.Nodes = ele(i,2:3);
    Element{i}.r=property.Radius;
    Element{i}.A=property.SectionArea;
    Element{i}.Iy=property.Iy;
    Element{i}.Iz=property.Iy;
    Element{i}.J=property.J;
    Element{i}.E=property.Es;
    Element{i}.G=property.Gs;
    Element{i}.v=property.v;
    startp = no(ele(i,2),2:4);
    endp = no(ele(i,3),2:4);
    trussDirection(i,:) = endp - startp;
end

% Assemble the stiffness matrix
NElements=length(Element);
Vreal = 0;
switch property.ElementType
    case 'Euler'
        for i=1:1:NElements
            x1 = Node{Element{i}.Nodes(1)}.x;
            y1 = Node{Element{i}.Nodes(1)}.y;
            z1 = Node{Element{i}.Nodes(1)}.z;
            x2 = Node{Element{i}.Nodes(2)}.x;
            y2 = Node{Element{i}.Nodes(2)}.y;
            z2 = Node{Element{i}.Nodes(2)}.z;
            DOFx1 = Node{Element{i}.Nodes(1)}.DOFx;
            DOFy1 = Node{Element{i}.Nodes(1)}.DOFy;
            DOFz1 = Node{Element{i}.Nodes(1)}.DOFz;
            DOFthetax1 = Node{Element{i}.Nodes(1)}.DOFthetax;
            DOFthetay1 = Node{Element{i}.Nodes(1)}.DOFthetay;
            DOFthetaz1 = Node{Element{i}.Nodes(1)}.DOFthetaz;
            DOFx2 = Node{Element{i}.Nodes(2)}.DOFx;
            DOFy2 = Node{Element{i}.Nodes(2)}.DOFy;
            DOFz2 = Node{Element{i}.Nodes(2)}.DOFz;
            DOFthetax2 = Node{Element{i}.Nodes(2)}.DOFthetax;
            DOFthetay2 = Node{Element{i}.Nodes(2)}.DOFthetay;
            DOFthetaz2 = Node{Element{i}.Nodes(2)}.DOFthetaz;
            Element{i}.L=sqrt((x2-x1)^2+(y2-y1)^2+(z2-z1)^2);
            Element{i}.ratio = Element{i}.L/2/Element{i}.r;
            Vreal = Vreal + property.SectionArea*Element{i}.L;
            Element{i}.Cx=(x2-x1)/Element{i}.L;
            Element{i}.Cy=(y2-y1)/Element{i}.L;
            Element{i}.Cz=(z2-z1)/Element{i}.L;

            [Element{i}.K,Element{i}.LocalK,Element{i}.T]=EulerBeam(Element{i}.E,Element{i}.Iy,Element{i}.Iz,Element{i}.A,Element{i}.L,Element{i}.G,Element{i}.J,Element{i}.Cx,Element{i}.Cy,Element{i}.Cz);
            Element{i}.DOFs=[DOFx1; DOFy1; DOFz1; DOFthetax1; DOFthetay1; DOFthetaz1; DOFx2; DOFy2; DOFz2; DOFthetax2; DOFthetay2; DOFthetaz2];
        end
    case 'Timoshenko'
        for i=1:1:NElements
            x1 = Node{Element{i}.Nodes(1)}.x;
            y1 = Node{Element{i}.Nodes(1)}.y;
            z1 = Node{Element{i}.Nodes(1)}.z;
            x2 = Node{Element{i}.Nodes(2)}.x;
            y2 = Node{Element{i}.Nodes(2)}.y;
            z2 = Node{Element{i}.Nodes(2)}.z;
            DOFx1 = Node{Element{i}.Nodes(1)}.DOFx;
            DOFy1 = Node{Element{i}.Nodes(1)}.DOFy;
            DOFz1 = Node{Element{i}.Nodes(1)}.DOFz;
            DOFthetax1 = Node{Element{i}.Nodes(1)}.DOFthetax;
            DOFthetay1 = Node{Element{i}.Nodes(1)}.DOFthetay;
            DOFthetaz1 = Node{Element{i}.Nodes(1)}.DOFthetaz;
            DOFx2 = Node{Element{i}.Nodes(2)}.DOFx;
            DOFy2 = Node{Element{i}.Nodes(2)}.DOFy;
            DOFz2 = Node{Element{i}.Nodes(2)}.DOFz;
            DOFthetax2 = Node{Element{i}.Nodes(2)}.DOFthetax;
            DOFthetay2 = Node{Element{i}.Nodes(2)}.DOFthetay;
            DOFthetaz2 = Node{Element{i}.Nodes(2)}.DOFthetaz;
            Element{i}.L = sqrt((x2-x1)^2+(y2-y1)^2+(z2-z1)^2);
            Vreal = Vreal + property.SectionArea*Element{i}.L;
            Element{i}.ratio = Element{i}.L/2/Element{i}.r;
            Element{i}.Cx = (x2-x1)/Element{i}.L;
            Element{i}.Cy = (y2-y1)/Element{i}.L;
            Element{i}.Cz = (z2-z1)/Element{i}.L;

            [Element{i}.K,Element{i}.LocalK,Element{i}.T] = TimoshenkoBeam(Element{i}.E,Element{i}.Iy,Element{i}.Iz,Element{i}.A,Element{i}.L,Element{i}.G,Element{i}.J,Element{i}.Cx,Element{i}.Cy,Element{i}.Cz);
            Element{i}.DOFs = [DOFx1; DOFy1; DOFz1; DOFthetax1; DOFthetay1; DOFthetaz1; DOFx2; DOFy2; DOFz2; DOFthetax2; DOFthetay2; DOFthetaz2];
        end
end
Kuc = MatrixAssembly(Element);
% property.Kuc = MatrixAssembly(Element);
property.RelativeRoh = Vreal/property.V;
% property.NewElement = Element;
end

