function Ba = BaCubicMatrix(mat,mapping)
% Multiscale Homogenization Code for Periodic Lattice Metamaterial
% Jin BO, Email:kinpoer@nuaa.edu.cn
%% Calculate Ba matrix
otherNodes = find(mapping(:,1) ~= mapping(:,2));
a1 = mat.a1;a2 = mat.a2;a3 = mat.a3;
th = 1e-5;
if norm(a3) ~= 0
    nodes = mat.NewNode;
    Batran =zeros(length(nodes),3);
    for j=1:length(otherNodes)
        indeNodes = nodes(mapping(otherNodes(j),2),:);% Independent Node
        deNodes = nodes(mapping(otherNodes(j),1),:);% Dependent Node
        if norm(indeNodes + a1 - deNodes)<th
            Batran(otherNodes(j),:)=[1 0 0];
        elseif norm(indeNodes + a1 + a2- deNodes)<th
            Batran(otherNodes(j),:)=[1 1 0];
        elseif norm(indeNodes +a2 - deNodes)<th
            Batran(otherNodes(j),:)=[0 1 0];
        elseif norm(indeNodes + a1 + a2 + a3 - deNodes)<th
            Batran(otherNodes(j),:)=[1 1 1];
        elseif norm(indeNodes + a3 - deNodes)<th
            Batran(otherNodes(j),:)=[0 0 1];
        elseif norm(indeNodes + a1 + a3 - deNodes)<th
            Batran(otherNodes(j),:)=[1 0 1];
        elseif norm(indeNodes + a2 + a3 - deNodes)<th
            Batran(otherNodes(j),:)=[0 1 1];
        end
    end
else
    nodes = mat.NewNode;
    Batran =zeros(length(nodes),2);
    for j=1:length(otherNodes)
        indeNodes = nodes(mapping(otherNodes(j),2),:);% Independent Node
        deNodes = nodes(mapping(otherNodes(j),1),:);% Dependent Node
        if norm(indeNodes + a1 - deNodes)<th
            Batran(otherNodes(j),:)=[1 0];
        elseif norm(indeNodes + a1 + a2- deNodes)<th
            Batran(otherNodes(j),:)=[1 1];
        elseif norm(indeNodes +a2 - deNodes)<th
            Batran(otherNodes(j),:)=[0 1];
        end
    end
end
% Periodic vector: represents the relationship between the resultant 
% external force and deformation. Ba represents the specific relative 
% displacement between the corresponding nodes represented by the periodic 
% vector.
if max(nodes(:,3))~=min(nodes(:,3))
    It = [eye(3);zeros(3)];
    [m,n] = size(Batran);
    Ba = [];
    for i=1:m
        temp = [];
        for j=1:n
            if Batran(i,j)~=0
                temp = [temp,Batran(i,j).*It];
            else
                temp = [temp,zeros(6,3)];
            end
        end
        Ba = [Ba;temp];
    end
else
    It = [1 0 ;0 1 ;0 0 ];
    [m,n] = size(Batran);
    Ba = [];
    for i=1:m
        temp = [];
        for j=1:n
            if Batran(i,j)~=0
                temp = [temp,Batran(i,j).*It];
            else
                temp = [temp,zeros(3,2)];
            end
        end
        Ba = [Ba;temp];
    end
end
end

