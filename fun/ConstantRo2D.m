function r = ConstantRo2D(nodes,elements,density,area)
% Multiscale Homogenization Code for Periodic Lattice Metamaterial
% Jin BO, Email:kinpoer@nuaa.edu.cn
%% Calculate beam radius with fixed relative density
d = nodes(elements(:,1),:)- nodes(elements(:,2),:);
Ltotal = 0;
for i =1:size(d,1)
Ltotal = norm(d(i,:))+Ltotal;
end
r = area*density/2/Ltotal;
end
