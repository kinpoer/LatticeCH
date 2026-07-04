function Plot3DCell(points,lineCon,lins)
% Multiscale Homogenization Code for Periodic Lattice Metamaterial
% Jin BO, Email:kinpoer@nuaa.edu.cn
%% Visualizing truss Cells
if size(points,2)==2
   points(:,3) = 0;
end
figure;
% Nodes Number
% for i=1:size(points,1)
%     scatter3(points(i,1),points(i,2),points(i,3),'k','filled');hold on;axis equal;
%     text(points(i,1),points(i,2),points(i,3),num2str(i),'HorizontalAlignment','center', 'VerticalAlignment','bottom','FontSize',12);
% end
xlabel('X');ylabel('Y');zlabel('Z');grid on;view(3)
% Elements Number
% for i=1:size(lineCon,1)
%     startP = points(lineCon(i,1),1:3);
%     endP = points(lineCon(i,2),1:3);
%     midP = [(startP(1) + endP(1))/2, (startP(2) + endP(2))/2, (startP(3) + endP(3))/2];
%     text(midP(1),midP(2),midP(3),num2str(i),'HorizontalAlignment','center','VerticalAlignment','bottom','Color','red','FontSize',10);
% end
x = points(:,1);y = points(:,2);z = points(:,3);
node1 = lineCon(:, 1);node2 = lineCon(:, 2);
XLines = [x(node1), x(node2)]';
YLines = [y(node1), y(node2)]';
ZLines = [z(node1), z(node2)]';
plot3(XLines, YLines, ZLines, lins,'LineWidth', 2);
hold on;axis equal;ax = gca;set(ax, 'TickDir', 'in');
xlabel('X');ylabel('Y');zlabel('Z');grid on;axis off;
end
