function [angle1,angle2,all_angles] = find_min_dividing_angle_deg3(Y_nodes_idx,G3,XY3,s_node)
% at any intersection, independent of degree, an edge makes two immediate
% angles; these angles are the abutting angles which needs to be checked
% for acute abutment or orthogonal abutment in the main script; they are
% returned as min_angle and max_angle. The third output lists all angles

%(tnode_1) +------+ (snode_1)-------chain-------(snode_2) +-----+(tnode_2)

chain = [Y_nodes_idx s_node]; % the 1st argument is tnode with deg3
%------------------------------------------------------------------

for i=1:numel(Y_nodes_idx(:,1))
 % listing edges connected to degree 3 node
 edges = [repelem(Y_nodes_idx(i,1),degree(G3,Y_nodes_idx(i,1)))' neighbors(G3,Y_nodes_idx(i,1))];
 % angle of edge from easting
 edges(:,3)=rad2deg(edgeAngle([XY3(edges(:,1),:) XY3(edges(:,2),:) ]));
 % sorting based on angles
 edges = sortrows(edges,3);
 
 
% locating row of chain within edges; stored in 'r'
 for j=1:numel(edges(:,1))
     A(j,1) = isequal(edges(j,1:2),chain);
 end
 r = find(A==1); 
 if isempty(r)
     for j=1:numel(edges(:,1))
       A(j,1) = isequal(edges(j,1:2),fliplr(chain));
     end
   r = find(A==1);  
 end
 
 if ~isempty(edges)
    edges(3,4) = 360-edges(end,3)+edges(1,3); 
    for j=1:numel(edges(:,1))-1 
            edges(j,4)=edges(j+1,3)-edges(j,3);
    end    
 end
 edges(:,3) = edges(:,4);
 edges(:,4) = [];
 Y_nodes_idx(i,2:numel(edges)+1) = reshape(edges',[numel(edges),1])';
 clearvars edges edge2lines 
 
end 

%------------------------------------------------------------------
Y_nodes_idx(:,11:12)=Y_nodes_idx(:,2:3);

% depending upon where the input chain is located (r), the angles that the
% edge makes with other edges is returned as max and min angle
switch r
    case 1
        angle1 = min([Y_nodes_idx(:,4);Y_nodes_idx(:,10)]); 
        angle2 = max([Y_nodes_idx(:,4);Y_nodes_idx(:,10)]);
    case 2
        angle1 = min([Y_nodes_idx(:,4);Y_nodes_idx(:,7)]); 
        angle2 = max([Y_nodes_idx(:,4);Y_nodes_idx(:,7)]);
    case 3 
        angle1 = min([Y_nodes_idx(:,7);Y_nodes_idx(:,10)]); 
        angle2 = max([Y_nodes_idx(:,7);Y_nodes_idx(:,10)]);
end

all_angles = [Y_nodes_idx(:,4) Y_nodes_idx(:,7) Y_nodes_idx(:,10)];
end

