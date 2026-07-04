function [nodes,elements] = ReadInp(filepath,filename)
% Multiscale Homogenization Code for Periodic Lattice Metamaterial
% Jin BO, Email:kinpoer@nuaa.edu.cn
%% Read nodes and elements in .Inp file
[fin,message] = fopen(strcat(filepath,filename),'r');%打开inp文件
if fin==-1
    error([message,': ',fileAbaqus]);
end
% Node info
tline = fgetl(fin);
while ~strncmp(tline,'*Node',5)
    tline = fgetl(fin);
end
tline = fgetl(fin);
tNodeCoor = textscan(tline,'%f','delimiter',',');
Nodes = textscan(fin,repmat('%f',1,numel(tNodeCoor{1})),'delimiter', ',');%批量读取节点信息
Nodes = [tNodeCoor{1}';cat(2,Nodes{:})];
Nodes(:,1) = [];
disp('Node information has been imported');
% Elements info
while ~strncmp(tline,'*Element',8)
    tline = fgetl(fin);
end
typeAbaqusEle = textscan(tline,'*Element, type=%s');%读取ABAQUS单元类型
typeAbaqusEle = typeAbaqusEle{1}{1};
switch typeAbaqusEle
    case {'S3R','S3'}
        typeEle = 'S3';
    case {'S4R','S4','S4R5'}
        typeEle = 'S4R';
    case {'CPS4','CPS4R','CPS4I'}
        typeEle = 'PS4';
    case {'CPS8','CPS8R'}
        typeEle = 'PS8';
    case {'C3D4','C3D4H'}
        typeEle = '3D4';
    case {'C3D8','C3D8R','C3D8H','C3D8I','C3D8RH','C3D8IH'}
        typeEle = '3D8';
    otherwise
        error('Non-supported type of element');
end
tline = fgetl(fin);
tEleNode = textscan(tline,'%f','delimiter',',');
Elements = textscan(fin,repmat('%f',1,numel(tEleNode{1})),'delimiter', ',');
Elements = [tEleNode{1}';cat(2,Elements{:})];
Elements(:,1) = [];
disp('Element information has been imported');
nodes = Nodes;elements = Elements;
end

