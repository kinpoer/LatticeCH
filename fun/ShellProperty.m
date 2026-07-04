function property = ShellProperty(nodes,elements,E,v,multi,boxtype)
% Multiscale Homogenization Code for Periodic Lattice Metamaterial
% Jin BO, Email:kinpoer@nuaa.edu.cn
%% Shell element information
dx = max(nodes(:,1))-min(nodes(:,1));
dy = max(nodes(:,2))-min(nodes(:,2));
dz = max(nodes(:,3))-min(nodes(:,3));

property.OriginNodes = nodes;
property.OriginElements = elements;

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