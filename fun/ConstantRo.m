function r = ConstantRo(nodes,elements,density,volume)
% Multiscale Homogenization Code for Periodic Lattice Metamaterial
% Jin BO, Email:kinpoer@nuaa.edu.cn
%% Calculate beam radius with fixed relative density
d = nodes(elements(:,1),:)- nodes(elements(:,2),:);
Ltotal = 0;
for i =1:length(d)
Ltotal = norm(d(i,:))+Ltotal;
end
r = sqrt(volume*density/Ltotal/pi);
end