function SavePlotAsPNG(fileName,subfolderName)
% Multiscale Homogenization Code for Periodic Lattice Metamaterial
% Jin BO, Email:kinpoer@nuaa.edu.cn
%% Save Figure
if nargin < 3
    pixels = '-r200';
end
currentFolder = pwd;
folderPath = fullfile(currentFolder, subfolderName);
if ~exist(folderPath, 'dir')
    mkdir(folderPath);
end
% The complete path to save the file
fullSavePath = fullfile(folderPath, fileName);
% Save picture in 'png' format
print(fullSavePath, '-dpng', pixels);
end
