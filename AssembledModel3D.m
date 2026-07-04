% Multiscale Homogenization Code for Periodic Lattice Metamaterial
% Jin BO, Email:kinpoer@nuaa.edu.cn
%% Initialization
clc;clear;close all;addpath('fun','src','data');

%% Combinatorial information
load AssembledUCs.mat
ICell = cell(length(combineInfo),1);
for i=1:length(combineInfo)
    % --First layer-- %
    T1 = combineInfo(i,1);
    N1 = UCs{T1}.N;
    L1 = UCs{T1}.L;
    A1 = combineInfo(i,2);
    if A1==2
        [N1,L1] = ArrayCubicCell(N1,L1,2,2,2,1e-8);
        N1 = N1/2;
    end
    % Plot3DCell(N1,L1,'b-');axis on;title('Layer1')

    % --Second layer-- %
    T2 = combineInfo(i,3);
    N2 = UCs{T2}.N;
    L2 = UCs{T2}.L;
    A2 = combineInfo(i,4);
    if A2==2
        [N2,L2] = ArrayCubicCell(N2,L2,2,2,2,1e-8);
        N2 = N2./2;
    end
    % Plot3DCell(N2,L2,'b-');axis on;title('Layer2')
    L2 = L2+length(N1)*ones(size(L2));

    % --Third layer-- %
    T3 = combineInfo(i,5);
    N3 = UCs{T3}.N;
    L3 = UCs{T3}.L;
    A3 = combineInfo(i,6);
    if A3==2
        [N3,L3] = ArrayCubicCell(N3,L3,2,2,2,1e-8);
        N3 = N3./2;
    end
    % Plot3DCell(N3,L3,'b-');axis on;title('Layer3')
    L3 = L3+(length(N1)+length(N2))*ones(size(L3));

    % --Merge all layers-- %
    N = [N1;N2;N3];
    L = [L1;L2;L3];
    [N,L] = MergeCloseNodesAndElements(N,L,1e-8);
    L = BreakElements(N,L);
    [N,L] = MergeCloseNodesAndElements(N,L,1e-8);
    [N,L] = CheckCross(N,L,0.1);
    ICell{i,1}.N = N;
    ICell{i,1}.L = L;
    ICell{i,1}.T = 'Assembled Modeling';
    Plot3DCell(N,L,'b-');axis on;title('Assembled Cell');
    close all
end

