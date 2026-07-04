function SavePlotAsFIG(fileName,subfolderName)
% Multiscale Homogenization Code for Periodic Lattice Metamaterial
% Jin BO, Email:kinpoer@nuaa.edu.cn
%% Save Figure
currentFolder = pwd;
folderPath = fullfile(currentFolder, subfolderName);
if ~exist(folderPath, 'dir')
    mkdir(folderPath);
end
% The complete path to save the file
fullSavePath = fullfile(folderPath, fileName);
% Save picture in 'fig' format
saveas(gcf, fullSavePath);
end