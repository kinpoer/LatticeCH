function property = BeamProperty(nodes,elements,r,E,v,multi,boxtype,crosssection)
% Multiscale Homogenization Code for Periodic Lattice Metamaterial
% Jin BO, Email:kinpoer@nuaa.edu.cn
%% Beam element information
dx = max(nodes(:,1))-min(nodes(:,1));
dy = max(nodes(:,2))-min(nodes(:,2));
dz = max(nodes(:,3))-min(nodes(:,3));

property.SectionType = crosssection;
property.OriginNodes = nodes;
property.OriginElements = elements;

switch crosssection
    case 'Circle'
        property.Radius = r(1);
        property.D =2*property.Radius;
        property.SectionArea = pi*r(1)^2;
        property.Iy = pi*(2*r)^4/64;
        property.Iz = property.Iy;
        property.J = property.Iz + property.Iy;
    case 'Rect'
        property.H = r(1)*2;
        property.b = r(2)*2;
        property.SectionArea = property.H*property.b;
        property.Iy = property.H*property.b^3/12;
        property.Iz = property.b*property.H^3/12;
        property.J = property.Iz + property.Iy;
    case 'Tube'
        property.MaxRadius = r(1);
        property.MinRadius = r(2);
        property.D =2*property.MaxRadius;
        property.d =2*property.MinRadius;
        property.SectionArea = pi*(r(1)^2-r(2)^2);
        property.Iy = pi*(property.D^4-property.d^4)/64;
        property.Iz = property.Iy;
        property.J = property.Iz + property.Iy;

end
property.Es = E;
property.v = v;
property.Gs = E/2/(1+v);

switch boxtype
case 'Cubic'
        if ~isempty(multi)
            property.multi = multi;
            property.a1 = property.multi*[-dx 0 0];
            property.a2 = property.multi*[0 -dy 0];
            property.a3 = property.multi*[0 0 dz];
        else
            property.multi = 1;
            property.a1 = property.multi*[-dx 0 0];
            property.a2 = property.multi*[0 -dy 0];
            property.a3 = property.multi*[0 0 dz];
        end
        property.BoxType = 'Cubic';
        if norm(property.a3)==0
            property.S = norm(property.a1)*norm(property.a2);
        else
            property.V = norm(property.a1)*norm(property.a2)*norm(property.a3);
        end
end
end
