% Multiscale Homogenization Code for Periodic Lattice Metamaterial
% Jin BO, Email:kinpoer@nuaa.edu.cn
%% Initialization
clc;clear;close all;addpath('fun','src','data');

%% case 1: for Truss or tube lattices
% Importing truss cells
load geo3D.mat %for truss UCs 
% Read geometry
opt = 1;
n = 30;rstrut = 2.2;rnode = 0.8;cellsize = 20;
numArray = 3;
nodes = tempFiles{opt,1}.N;struts = tempFiles{opt,1}.L;
Plot3DCell(nodes,struts,'b-');
% Caculate volumetric distance field
[x,y,z] = meshgrid(-cellsize/2:cellsize/2/n:cellsize/2);
rstrut = rstrut*ones(length(struts),1);
U = TrussLatticeFV(x,y,z,nodes,struts,rstrut,rnode);
% Construct a level set function
% For Tube lattices
t = 0.3;U1 = -(U+t).*(U-t);% thickness
[fs,v] = isosurface(x,y,z,U1,0);
[fc,v2,~] = isocaps(x,y,z,U1,0);
fn = [fs ; fc+length(v(:,1))];
vn = [v ; v2];
Plot3DMesh(vn,fn,1);
% For Truss lattices
U2 = U;
[fs,v] = isosurface(x,y,z,-U2,0);
[fc,v2,~] = isocaps(x,y,z,-U2,0);
fn = [fs ; fc+length(v(:,1))];
vn = [v ; v2];
Plot3DMesh(vn,fn,1);
% Generate .stl file for AM
% [vn,fn] = ArrayCellMesh(vn,fn,numArray,numArray,numArray);% array UC
[stlVol,~]=StlVolume(vn',fn');
relativeDensity = abs(stlVol)/20/20/20;vn = 0.001*vn;
StlWrite('TrussLattice.stl',fn,vn); % generate .STL file of the final lattice

%% case 2: for Plate or Shell lattices
% Note: implicit formulas of volumetric distance fields suit for TPMS lattices
% Importing Plate cells
load PlateGeo.mat
% Read geometry
n = 30;
numArray = 3;elements = faces;
Plot3DMesh(nodes,elements,1);
cellsize = max(nodes)-min(nodes);cellsize = cellsize(1);
% Caculate volumetric distance field
[x,y,z] = meshgrid(-cellsize/2:cellsize/2/n:cellsize/2);% for cubic bound
U = SurfaceLatticeFV(x,y,z,nodes,elements);
% Construct a level set function
% For surface lattices
t = 0.3;U3 = (U-t);
[fs,v] = isosurface(x,y,z,-U3,0);
[fc,v2,~] = isocaps(x,y,z,-U3,0);
fn = [fs ; fc+length(v(:,1))];
vn = [v ; v2];
Plot3DMesh(vn,fn,1);
StlWrite('SurfaceLattice.stl',fn,vn);


