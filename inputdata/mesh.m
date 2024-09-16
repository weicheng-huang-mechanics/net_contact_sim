clear all;
close all;
clc;

l = 8;
dl = 6/24;
interNodeNum = 2;
h = sqrt(3) * dl / 2;

nv = l / dl;

nodes = zeros(3,3);

nodes(1,1) = l * cos(pi/3);
nodes(1,2) = l * sin(pi/3);
nodes(1,3) = 0.0;

figure(1)
hold on;

temp = 2;
for i = 2:nv+1
    startPoint = zeros(1,3);
    
    startPoint(1) = nodes(1,1) - (i-1) * dl / 2;
    startPoint(2) = nodes(1,2) - (i-1) * h;
    startPoint(3) = nodes(1,3);
    
    for j = 1:i
        nodes(temp, 1) = startPoint(1) + (j-1) * dl;
        nodes(temp, 2) = startPoint(2);
        nodes(temp, 3) = 0.0;
        temp = temp + 1;
    end
end

conStep = 1;
consNodes = zeros(3,3);

consNodes(1,1) = l * cos(pi/3);
consNodes(1,2) = l * sin(pi/3);
consNodes(1,3) = 0.0;

temp = 1;
for i = 1:conStep:nv+1
    startPoint = zeros(1,3);
    
    startPoint(1) = nodes(1,1) - (i-1) * dl / 2;
    startPoint(2) = nodes(1,2) - (i-1) * h;
    startPoint(3) = nodes(1,3);
    
    for j = 1:conStep:i
        consNodes(temp, 1) = startPoint(1) + (j-1) * dl;
        consNodes(temp, 2) = startPoint(2);
        consNodes(temp, 3) = 0.0;
        temp = temp + 1;
    end
end

nodesTotal = zeros(3,3);

theta = 2 * pi / 6;

[numVertex,~] = size(nodes);

temp = 1;
for i = 1:6
    for j = 1:numVertex
        nodesTotal(temp,1) = nodes(j,1) * cos((i-1)*theta) + nodes(j,2) * sin((i-1)*theta);
        nodesTotal(temp,2) =-nodes(j,1) * sin((i-1)*theta) + nodes(j,2) * cos((i-1)*theta);
        nodesTotal(temp,3) = 0.0;
        temp = temp + 1;
    end
end

[numVertexCons,~] = size(consNodes);
consNodesTotal = zeros(3,3);
temp = 1;
for i = 1:6
    for j = 1:numVertexCons
        consNodesTotal(temp,1) = consNodes(j,1) * cos((i-1)*theta) + consNodes(j,2) * sin((i-1)*theta);
        consNodesTotal(temp,2) =-consNodes(j,1) * sin((i-1)*theta) + consNodes(j,2) * cos((i-1)*theta);
        consNodesTotal(temp,3) = 0.0;
        temp = temp + 1;
    end
end

view(0,90);
axis equal;
box on;

[totalNv, ~] = size(nodesTotal);

nodesNew = zeros(1,3);

temp = 2;
for i = 1:totalNv
    [currentNv, ~] = size(nodesNew);
    
    node1(1) = nodesTotal(i,1);
    node1(2) = nodesTotal(i,2);
    
    flag = 0;
    for j = 1:currentNv
        node2(1) = nodesNew(j,1);
        node2(2) = nodesNew(j,2);
        
        if abs(norm(node2-node1)) < 1e-4
            flag = 1;
        end
    end
    
    if flag == 0
        nodesNew(temp,1) = nodesTotal(i,1);
        nodesNew(temp,2) = nodesTotal(i,2);
        nodesNew(temp,3) = 0.0;
        temp = temp + 1;
    end
end


[totalNvCons, ~] = size(consNodesTotal);

nodesNewCons = zeros(1,3);

temp = 2;
for i = 1:totalNvCons
    [currentNv, ~] = size(nodesNewCons);
    
    node1(1) = consNodesTotal(i,1);
    node1(2) = consNodesTotal(i,2);
    
    flag = 0;
    for j = 1:currentNv
        node2(1) = nodesNewCons(j,1);
        node2(2) = nodesNewCons(j,2);
        
        if abs(norm(node2-node1)) < 1e-4
            flag = 1;
        end
    end
    
    if flag == 0
        nodesNewCons(temp,1) = consNodesTotal(i,1);
        nodesNewCons(temp,2) = consNodesTotal(i,2);
        nodesNewCons(temp,3) = 0.0;
        temp = temp + 1;
    end
end
plot3(nodesNewCons(:,1),nodesNewCons(:,2),nodesNewCons(:,3),'blacko');

[rowsConsNode,~] = size(nodesNewCons);
[totalNv, ~] = size(nodesNew);

consIndex = zeros(3,1);
temp = 1;
for i = 1:rowsConsNode
    
    node1(1) = nodesNewCons(i,1);
    node1(2) = nodesNewCons(i,2);
    for j = 1:totalNv
        node2(1) = nodesNew(j,1);
        node2(2) = nodesNew(j,2);
        if norm(node1-node2) < 1e-5
            consIndex(temp, 1) = j;
            %consIndex(temp, 2) = norm(node1) / l;
            temp = temp + 1;
        end
    end
end

view(0,90);
axis equal;
box on;

epsilon = 1e-2;

nodes = nodesNew;

[numVertex,~] = size(nodes);

temp = 1;
stretchIndex = zeros(3,2);
for i = 1:numVertex
    for j = i+1:numVertex
        node1 = nodes(i,:);
        node2 = nodes(j,:);
        
        distance = norm(node2 - node1);
        
        tangent = (node2 - node1) / norm(node2 - node1);
        
        if distance < 1.1 * dl && distance > 0.9 * dl
            
            consTangent = zeros(3,1);
        
            if ( node1(2) - sqrt(3)*node1(1) <= epsilon && node1(2) >= -epsilon && node1(2) + sqrt(3)*node1(1) - sqrt(3)*l <= -epsilon)
                if ( node2(2) - sqrt(3)*node2(1) <= epsilon && node2(2) >= -epsilon && node2(2) + sqrt(3)*node2(1) - sqrt(3)*l <= -epsilon)
%                     plot3(node1(1),node1(2),node1(3),'blacko');
%                     plot3(node2(1),node2(2),node2(3),'blacko');
                    consTangent(1) = 1.0;
                    consTangent(2) = -sqrt(3);
                    consTangent(3) = 0.0;
            
                    consTangent = consTangent / norm(consTangent);
                end
            end
        
            if ( node1(2) - sqrt(3)*node1(1) >= -epsilon && node1(2) + sqrt(3)*node1(1) >= -epsilon && node1(2) - sqrt(3)/2*l <= -epsilon)
                if ( node2(2) - sqrt(3)*node2(1) >= -epsilon && node2(2) + sqrt(3)*node2(1) >= -epsilon && node2(2) - sqrt(3)/2*l <= -epsilon)
%                     plot3(node1(1),node1(2),node1(3),'blacko');
%                     plot3(node2(1),node2(2),node2(3),'blacko');
                    consTangent(1) = 1.0;
                    consTangent(2) = 0.0;
                    consTangent(3) = 0.0;
            
                    consTangent = consTangent / norm(consTangent);
                end
            end
        
            if ( node1(2) + sqrt(3)*node1(1) < epsilon && node1(2) > -epsilon && node1(2) - sqrt(3)*node1(1) - sqrt(3)*l < -epsilon)
                if ( node2(2) + sqrt(3)*node2(1) < epsilon && node2(2) > -epsilon && node2(2) - sqrt(3)*node2(1) - sqrt(3)*l < -epsilon)
%                     plot3(node1(1),node1(2),node1(3),'blacko');
%                     plot3(node2(1),node2(2),node2(3),'blacko');
                    consTangent(1) = 1;
                    consTangent(2) = sqrt(3);
                    consTangent(3) = 0;
            
                    consTangent = consTangent / norm(consTangent);
                end
            end
        
            if ( node1(2) - sqrt(3)*node1(1) > -epsilon && node1(2) < epsilon && node1(2) + sqrt(3)*node1(1) + sqrt(3)*l > epsilon)
                if ( node2(2) - sqrt(3)*node2(1) > -epsilon && node2(2) < epsilon && node2(2) + sqrt(3)*node2(1) + sqrt(3)*l > epsilon)
%                     plot3(node1(1),node1(2),node1(3),'blacko');
%                     plot3(node2(1),node2(2),node2(3),'blacko');
                    consTangent(1) = 1.0;
                    consTangent(2) = -sqrt(3);
                    consTangent(3) = 0.0;
            
                    consTangent = consTangent / norm(consTangent);
                end
            end
        
            if ( node1(2) - sqrt(3)*node1(1) < epsilon && node1(2) + sqrt(3)*node1(1) < epsilon && node1(2) + sqrt(3)/2*l > epsilon)
                if ( node2(2) - sqrt(3)*node2(1) < epsilon && node2(2) + sqrt(3)*node2(1) < epsilon && node2(2) + sqrt(3)/2*l > epsilon)
%                     plot3(node1(1),node1(2),node1(3),'blacko');
%                     plot3(node2(1),node2(2),node2(3),'blacko');
                    
                    consTangent(1) = 1.0;
                    consTangent(2) = 0.0;
                    consTangent(3) = 0.0;
            
                    consTangent = consTangent / norm(consTangent);
                end
            end
        
            if ( node1(2) + sqrt(3)*node1(1) > -epsilon && node1(2) < epsilon && node1(2) - sqrt(3)*node1(1) + sqrt(3)*l > epsilon)
                if ( node2(2) + sqrt(3)*node2(1) > -epsilon && node2(2) < epsilon && node2(2) - sqrt(3)*node2(1) + sqrt(3)*l > epsilon)
%                     plot3(node1(1),node1(2),node1(3),'blacko');
%                     plot3(node2(1),node2(2),node2(3),'blacko');
                    consTangent(1) = 1;
                    consTangent(2) = sqrt(3);
                    consTangent(3) = 0;
            
                    consTangent = consTangent / norm(consTangent);
                end
            end
            
            if (abs(dot(consTangent,tangent)) < 0.999)
                stretchIndex(temp, 1) = i;
                stretchIndex(temp, 2) = j;
                temp = temp + 1;
            end           
        end
    end
end

[stretchNum, ~] = size(stretchIndex);

stretchIndex2 = zeros(3,2);
temp2 = 1;

temp = numVertex + 1;

for i = 1:stretchNum
    index1 = stretchIndex(i,1);
    index2 = stretchIndex(i,2);
    node1 = nodes(index1,:);
    node2 = nodes(index2,:);
    tangent = (node2 - node1) / norm(node2 - node1);
    
    len = norm(node2 - node1);
    deltal = len / (interNodeNum + 1);
    
    for j = 1:interNodeNum
        
        nodes(temp,1) = node1(1) + tangent(1) * j * deltal;
        nodes(temp,2) = node1(2) + tangent(2) * j * deltal;
        nodes(temp,3) = 0.0;
        
        if j == 1
            stretchIndex2(temp2, 1) = index1;
            stretchIndex2(temp2, 2) = temp;
            temp2 = temp2 + 1;
        end
        
        if j == interNodeNum
            stretchIndex2(temp2, 1) = temp;
            stretchIndex2(temp2, 2) = index2;
            temp2 = temp2 + 1;
        end
        
        if j < interNodeNum
            stretchIndex2(temp2, 1) = temp;
            stretchIndex2(temp2, 2) = temp + 1;
            temp2 = temp2 + 1;
        end
        
        temp = temp + 1;
    end
end

[numVertex,~] = size(nodes);

stretchIndex = stretchIndex2;

[stretchNum, ~] = size(stretchIndex);

for i = 1:stretchNum
    index1 = stretchIndex(i,1);
    index2 = stretchIndex(i,2);
    
    plot3([nodes(index1,1), nodes(index2,1)],[nodes(index1,2), nodes(index2,2)],[nodes(index1,3), nodes(index2,3)],'r-');
end

dx = 0.5;
dy = 0.5;

connerIndex = zeros(6,1);
tempConner = 1;
for i = 1:numVertex
    
    nodeLocal = nodes(i,:);
    
    if (sqrt(nodeLocal(1)*nodeLocal(1)+nodeLocal(2)*nodeLocal(2))<l*1.01 && sqrt(nodeLocal(1)*nodeLocal(1)+nodeLocal(2)*nodeLocal(2))>l*0.99)
        b = num2str(i);
        c = cellstr(b);
        text(nodeLocal(1)+dx, nodeLocal(2)+dy, c);
        connerIndex(tempConner) = i;
        tempConner = tempConner + 1;
    end
end


for i = 1:numVertex
    %plot3(nodes(i,1),nodes(i,2),nodes(i,3),'blacko');
    
    if ( abs(sqrt(nodes(i,1)*nodes(i,1)+nodes(i,2)*nodes(i,2)) - l) < 0.001 )
        b = num2str(i);
        c = cellstr(b);
        %text(nodes(i,1)+dx, nodes(i,2)+dy, c);
    end
end
%plot3(nodes(:,1),nodes(:,2),nodes(:,3),'blacko');


% edgeNode = zeros(3,2);
% edgeNode(1,1) = 1;
% 
% temp3 = 2;
% for i = 1:numVertex
%     if ( abs(nodes(i,1)*sqrt(3) - nodes(i,2)) < 1e-3 && nodes(i,1) > 0.001 )
%         %plot3(nodes(i,1),nodes(i,2),nodes(i,3),'blacko');
%         edgeNode(temp3,1) = i;
%         edgeNode(temp3,2) = sqrt(nodes(i,1)*nodes(i,1)+nodes(i,2)*nodes(i,2));
%         temp3 = temp3 + 1;
%     end
% end
% 
% edgeConsIndex = zeros(temp3-1,6);
% 
% aaa = sortrows(edgeNode,2);
% edgeConsIndex(:,1) = aaa(:,1);
% 
% temp3 = 2;
% for i = 1:numVertex
%     if ( abs(nodes(i,1)*sqrt(3) + nodes(i,2)) < 1e-3 && nodes(i,1) < -0.001 )
%         %plot3(nodes(i,1),nodes(i,2),nodes(i,3),'blacko');
%         edgeNode(temp3,1) = i;
%         edgeNode(temp3,2) = sqrt(nodes(i,1)*nodes(i,1)+nodes(i,2)*nodes(i,2));
%         temp3 = temp3 + 1;
%     end
% end
% 
% aaa = sortrows(edgeNode,2);
% edgeConsIndex(:,2) = aaa(:,1);
% 
% temp3 = 2;
% for i = 1:numVertex
%     if ( abs(nodes(i,2)) < 1e-3 && nodes(i,1) < -0.001 )
%         %plot3(nodes(i,1),nodes(i,2),nodes(i,3),'blacko');
%         edgeNode(temp3,1) = i;
%         edgeNode(temp3,2) = sqrt(nodes(i,1)*nodes(i,1)+nodes(i,2)*nodes(i,2));
%         temp3 = temp3 + 1;
%     end
% end
% 
% aaa = sortrows(edgeNode,2);
% edgeConsIndex(:,3) = aaa(:,1);
% 
% temp3 = 2;
% for i = 1:numVertex
%     if ( abs(nodes(i,1)*sqrt(3) - nodes(i,2)) < 1e-3 && nodes(i,1) < -0.001 )
%         %plot3(nodes(i,1),nodes(i,2),nodes(i,3),'blacko');
%         edgeNode(temp3,1) = i;
%         edgeNode(temp3,2) = sqrt(nodes(i,1)*nodes(i,1)+nodes(i,2)*nodes(i,2));
%         temp3 = temp3 + 1;
%     end
% end
% 
% aaa = sortrows(edgeNode,2);
% edgeConsIndex(:,4) = aaa(:,1);
% 
% temp3 = 2;
% for i = 1:numVertex
%     if ( abs(nodes(i,1)*sqrt(3) + nodes(i,2)) < 1e-3 && nodes(i,1) > 0.001 )
%         %plot3(nodes(i,1),nodes(i,2),nodes(i,3),'blacko');
%         edgeNode(temp3,1) = i;
%         edgeNode(temp3,2) = sqrt(nodes(i,1)*nodes(i,1)+nodes(i,2)*nodes(i,2));
%         temp3 = temp3 + 1;
%     end
% end
% 
% aaa = sortrows(edgeNode,2);
% edgeConsIndex(:,5) = aaa(:,1);
% 
% temp3 = 2;
% for i = 1:numVertex
%     if ( abs(nodes(i,2)) < 1e-3 && nodes(i,1) > 0.001 )
%         %plot3(nodes(i,1),nodes(i,2),nodes(i,3),'blacko');
%         edgeNode(temp3,1) = i;
%         edgeNode(temp3,2) = sqrt(nodes(i,1)*nodes(i,1)+nodes(i,2)*nodes(i,2));
%         temp3 = temp3 + 1;
%     end
% end
% 
% aaa = sortrows(edgeNode,2);
% edgeConsIndex(:,6) = aaa(:,1);


[rowsConsNode,~] = size(nodesNewCons);

edgeNode2 = zeros(3,2);

temp3 = 1;
for i = 1:rowsConsNode
    if ( abs( nodesNewCons(i,2) + sqrt(3)*nodesNewCons(i,1) - sqrt(3)*l ) < 1e-3 )
        %plot3(nodes(i,1),nodes(i,2),nodes(i,3),'blacko');
        edgeNode2(temp3,1) = i;
        edgeNode2(temp3,2) = nodesNewCons(i,1);
        temp3 = temp3 + 1;
    end
end

edgeConsIndex2 = zeros(temp3-1,6);

aaaa = sortrows(edgeNode2,2);
edgeConsIndex2(:,1) = aaaa(:,1);

temp3 = 1;
for i = 1:rowsConsNode
    if ( abs( nodesNewCons(i,2) - sqrt(3)/2*l ) < 1e-3 )
        %plot3(nodes(i,1),nodes(i,2),nodes(i,3),'blacko');
        edgeNode2(temp3,1) = i;
        edgeNode2(temp3,2) = nodesNewCons(i,1);
        temp3 = temp3 + 1;
    end
end
aaaa = sortrows(edgeNode2,2);
edgeConsIndex2(:,2) = aaaa(:,1);

temp3 = 1;
for i = 1:rowsConsNode
    if ( abs( nodesNewCons(i,2) - sqrt(3)*nodesNewCons(i,1) - sqrt(3)*l ) < 1e-3  )
        %plot3(nodes(i,1),nodes(i,2),nodes(i,3),'blacko');
        edgeNode2(temp3,1) = i;
        edgeNode2(temp3,2) = nodesNewCons(i,1);
        temp3 = temp3 + 1;
    end
end
aaaa = sortrows(edgeNode2,2);
edgeConsIndex2(:,3) = aaaa(:,1);

temp3 = 1;
for i = 1:rowsConsNode
    if ( abs( nodesNewCons(i,2) + sqrt(3)*nodesNewCons(i,1) + sqrt(3)*l ) < 1e-3 )
        %plot3(nodes(i,1),nodes(i,2),nodes(i,3),'blacko');
        edgeNode2(temp3,1) = i;
        edgeNode2(temp3,2) = nodesNewCons(i,1);
        temp3 = temp3 + 1;
    end
end
aaaa = sortrows(edgeNode2,2);
edgeConsIndex2(:,4) = aaaa(:,1);

temp3 = 1;
for i = 1:rowsConsNode
    if ( abs( nodesNewCons(i,2) + sqrt(3)/2*l ) < 1e-3 )
        %plot3(nodes(i,1),nodes(i,2),nodes(i,3),'blacko');
        edgeNode2(temp3,1) = i;
        edgeNode2(temp3,2) = nodesNewCons(i,1);
        temp3 = temp3 + 1;
    end
end
aaaa = sortrows(edgeNode2,2);
edgeConsIndex2(:,5) = aaaa(:,1);

temp3 = 1;
for i = 1:rowsConsNode
    if ( abs( nodesNewCons(i,2) - sqrt(3)*nodesNewCons(i,1) + sqrt(3)*l ) < 1e-3  )
        %plot3(nodes(i,1),nodes(i,2),nodes(i,3),'blacko');
        edgeNode2(temp3,1) = i;
        edgeNode2(temp3,2) = nodesNewCons(i,1);
        temp3 = temp3 + 1;
    end
end
aaaa = sortrows(edgeNode2,2);
edgeConsIndex2(:,6) = aaaa(:,1);




temp = 1;
bendingIndex = zeros(3,3);

for i = 1:stretchNum
    for j = i+1:stretchNum
        index1 = stretchIndex(i,1);
        index2 = stretchIndex(i,2);
        index3 = stretchIndex(j,1);
        index4 = stretchIndex(j,2);
        
        if index1 == index3 || index1 == index4 || index2 == index3 || index2 == index4
            node1 = nodes(index1,:);
            node2 = nodes(index2,:);
            node3 = nodes(index3,:);
            node4 = nodes(index4,:);
            
            tangent1 = node2 - node1;
            tangent2 = node4 - node3;
            
            tangent1 = tangent1 / norm(tangent1);
            tangent2 = tangent2 / norm(tangent2);
            
            dotProduct = dot(tangent2,tangent1);
            
            if abs(dotProduct) > 0.999
                if index1 == index3
                    bendingIndex(temp, 1) = index2;
                    bendingIndex(temp, 2) = index1;
                    bendingIndex(temp, 3) = index4;
                    
                    temp = temp + 1;
                end
                
                if index2 == index3
                    bendingIndex(temp, 1) = index1;
                    bendingIndex(temp, 2) = index2;
                    bendingIndex(temp, 3) = index4;
                    
                    temp = temp + 1;
                end
                
                if index1 == index4
                    bendingIndex(temp, 1) = index2;
                    bendingIndex(temp, 2) = index1;
                    bendingIndex(temp, 3) = index3;
                    
                    temp = temp + 1;
                end
                
                if index2 == index4
                    bendingIndex(temp, 1) = index1;
                    bendingIndex(temp, 2) = index2;
                    bendingIndex(temp, 3) = index3;
                    
                    temp = temp + 1;
                end
            end
        end
    end
end

% constraintStep = 4;
% [edgeNodeSize,~] = size(edgeConsIndex);
% 
% consIndex = zeros(3,6);
% tempNew = 1;
% for i = 1:constraintStep:edgeNodeSize
%     consIndex(tempNew,:) = edgeConsIndex(i,:);
%     tempNew = tempNew + 1;
% end
% 
% [consNum,~] = size(consIndex);
% for i = 1:consNum
%     for j = 1:6
%         plot3(nodes(consIndex(i,j),1),nodes(consIndex(i,j),2),nodes(consIndex(i,j),3),'blacko');
%         b = num2str(consIndex(i,j));
%         c = cellstr(b);
%         text(nodes(consIndex(i,j),1)+dx, nodes(consIndex(i,j),2)+dy, c);
%     end
% end

[edgeNodeSize,~] = size(edgeConsIndex2);
constraintStep = 1;

consIndex2 = zeros(3,6);
tempNew = 1;
for i = 1:constraintStep:edgeNodeSize
    consIndex2(tempNew,:) = edgeConsIndex2(i,:);
    tempNew = tempNew + 1;
end

[consNum,~] = size(consIndex2);
for i = 1:consNum
    for j = 1:6
        %plot3(nodesNewCons(consIndex2(i,j),1),nodesNewCons(consIndex2(i,j),2),nodesNewCons(consIndex2(i,j),3),'greens');
        b = num2str(consIndex2(i,j));
        c = cellstr(b);
        %text(nodes(consIndex2(i,j),1)+dx, nodes(consIndex2(i,j),2)+dy, c);
    end
end

[rows, cols] = size(consIndex2);
boundaryCons = zeros(rows*cols,1);

temp = 1;
for i = 1:rows
    for j = 1:cols
        boundaryCons(temp) = consIndex2(i,j);
        temp = temp + 1;
    end
end

for i = 1:rows*cols
    plot3(nodesNewCons(boundaryCons(i),1),nodesNewCons(boundaryCons(i),2),nodesNewCons(boundaryCons(i),3),'greens');
    b = num2str(consIndex(boundaryCons(i)));
    c = cellstr(b);
    text(nodesNewCons(boundaryCons(i),1)+dx, nodesNewCons(boundaryCons(i),2)+dy, c);
end


[rows,cols] = size(consIndex);

figure(2);
hold on;

connerIndex = zeros(6,1);
temppppp = 1;
for i = 1:rows
    index = consIndex(i);
    nodelocal = nodes(index,:);
    plot3(nodelocal(1),nodelocal(2),nodelocal(3),'ro');
    b = num2str(index);
    c = cellstr(b);
    text(nodelocal(1)+dx, nodelocal(2)+dy, c);
    
    if (sqrt(nodelocal(1)*nodelocal(1)+nodelocal(2)*nodelocal(2)) > 0.99*l)
        connerIndex(temppppp) = index;
        temppppp = temppppp + 1;
    end
    
    
end


dlmwrite('connerIndex.txt',connerIndex-1,'delimiter',' ');
dlmwrite('consInput.txt',consIndex-1,'delimiter',' ');
dlmwrite('nodesInput.txt',nodes,'delimiter',' ');
dlmwrite('edgeInput.txt',stretchIndex-1,'delimiter',' ');
dlmwrite('bendingInput.txt',bendingIndex-1,'delimiter',' ');

