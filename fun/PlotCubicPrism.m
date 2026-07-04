function PlotCubicPrism(vertices)
% Multiscale Homogenization Code for Periodic Lattice Metamaterial
% Jin BO, Email:kinpoer@nuaa.edu.cn
%% Draw a prism
line = [1 2;2 3;3 4;4 1;5 6;6 7;7 8;8 5;1 5;2 6;3 7;4 8];
% for i=1:length(vertices)
%     scatter3(vertices(i,1),vertices(i,2),vertices(i,3),'filled');hold on;axis equal;
%     ax = gca;set(ax, 'TickDir', 'in');box on;xlabel('X');ylabel('Y');zlabel('Z');grid on
%     text(vertices(i,1),vertices(i,2),vertices(i,3), num2str(i), 'HorizontalAlignment','center', 'VerticalAlignment','bottom');
% end
for i=1:length(line)
    sp = vertices(line(i,1),:);ep = vertices(line(i,2),:);
    plot3([sp(1);ep(1)],[sp(2);ep(2)],[sp(3);ep(3)],'LineWidth',2,'Color','k','LineStyle','--');hold on
end
axis equal;grid on;
xlabel('X');ylabel('Y');zlabel('Z');view(3);
end
