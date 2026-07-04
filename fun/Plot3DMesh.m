function Plot3DMesh(nodes,elements,facealpha)
% Multiscale Homogenization Code for Periodic Lattice Metamaterial
% Jin BO, Email:kinpoer@nuaa.edu.cn
%% Plot 3D mesh
figure;
patch('vertices',nodes,'Faces',elements,'FaceColor','r','EdgeColor','black','facealpha',facealpha)%,'EdgeColor','none'
view(3)
axis off;axis equal
end