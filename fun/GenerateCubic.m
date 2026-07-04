function coords = GenerateCubic(sideLength,turn)
% Multiscale Homogenization Code for Periodic Lattice Metamaterial
% Jin BO, Email:kinpoer@nuaa.edu.cn
%% Calculate and visualizethe node coordinates of the cube
halfSide = sideLength/2;
coords = [-halfSide, -halfSide, -halfSide;
    halfSide, -halfSide, -halfSide;
    halfSide,  halfSide, -halfSide;
    -halfSide,  halfSide, -halfSide;
    -halfSide, -halfSide,  halfSide;
    halfSide, -halfSide,  halfSide;
    halfSide,  halfSide,  halfSide;
    -halfSide,  halfSide,  halfSide];

edges = [1, 2; 2, 3; 3, 4; 4, 1;
    5, 6; 6, 7; 7, 8; 8, 5;
    1, 5; 2, 6; 3, 7; 4, 8];

linw = 2;
if strcmp(turn,'on')
    for i = [11 2 3]
        plot3([coords(edges(i, 1), 1), coords(edges(i, 2), 1)], ...
            [coords(edges(i, 1), 2), coords(edges(i, 2), 2)], ...
            [coords(edges(i, 1), 3), coords(edges(i, 2), 3)], 'k--', 'LineWidth', linw);hold on
    end
    for i = [1,4:10,12]
        plot3([coords(edges(i, 1), 1), coords(edges(i, 2), 1)], ...
            [coords(edges(i, 1), 2), coords(edges(i, 2), 2)], ...
            [coords(edges(i, 1), 3), coords(edges(i, 2), 3)], 'k--', 'LineWidth', linw);hold on
    end
    view(3)
end
end
