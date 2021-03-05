function [kink_connect, possible_edges] = find_kinked_M_Node_connections(kinked_M_nodes,XY_new,dt2,G3)
% this function parses through the kinked_M_nodes, finds the adjoining
% trielements pertaining to each kinked_M_node, within the triangulation.
% The adjoining trielements are reduced to vertices and a check is
% performed w.r.t the list of kinked_M_nodes vertices, if there are any
% other kinked_M_nodes vertices within the adjoining trielements. The 
% possible connections are stored in a cell array and returned.

X2 = XY_new(:,1:2);

kinked_M_nodes_XY = [XY_new(kinked_M_nodes,1) XY_new(kinked_M_nodes,2)];
% cell array to store the possible kinked connections
kink_connect=cell(numel(kinked_M_nodes_XY(:,1)),1);

for i=1:numel(kinked_M_nodes_XY (:,1))
    disp(i)
    
    % all tri-elements adjacent to the kinked_M_node, vertex attachments
    vertex_attachments = cell2mat(vertexAttachments(dt2,kinked_M_nodes(i,:)))' ;
    
    % extracting points of all the adjacent vertices, reshaping into a
    % column vector, applying unique operation, and then extracting
    % coordinates, removing current kinked_M_node from the P_vertex_attachments
    P_vertex_attachments = setdiff(X2(unique(reshape(dt2.ConnectivityList(vertex_attachments,:),...
        [numel(dt2.ConnectivityList(vertex_attachments,:)),1])),:),kinked_M_nodes_XY(i,:),'rows');
    
    P_vertex_attachments_idx =  setdiff(unique(reshape(dt2.ConnectivityList(vertex_attachments,:),...
        [numel(dt2.ConnectivityList(vertex_attachments,:)),1])), kinked_M_nodes(i,:),'rows');
       
    % check if these points have commonality within the kinked_M_Nodes_XY
    % list
    if ~isempty(intersect(kinked_M_nodes_XY,P_vertex_attachments,'rows'))
        [kink_connect{i,1},ia,~]=intersect(kinked_M_nodes_XY,P_vertex_attachments,'rows'); 
        kink_connect{i,2}=kinked_M_nodes(ia,:);
        kink_connect{i,3}=kinked_M_nodes_XY(i,:);
        kink_connect{i,4}=kinked_M_nodes(i,:);
%         kink_connect{i,5}=Combined_Edges_kinked_idx(i,1);
    end       
    clearvars vertex_attachments P_vertex_attachments P_vertex_attachments_idx ia
end   

% removing emnpty cell arrays; those with no possible kink connections
kink_connect = kink_connect(~cellfun(@isempty, kink_connect(:,1)), :);

j=1;
for i=1:numel(kink_connect(:,1)) 
 if numel(kink_connect{i,1})>1
    no_of_points=numel(kink_connect{i,2});   
    for k=1:no_of_points        
      possible_edges(j,1)=min(kink_connect{i,4},kink_connect{i,2}(k,1)); 
      possible_edges(j,2)=max(kink_connect{i,4},kink_connect{i,2}(k,1));
      j=j+1;
    end       
 else
    possible_edges(j,1) = min(kink_connect{i,2},kink_connect{i,4});
    possible_edges(j,2) = max(kink_connect{i,2},kink_connect{i,4});
    j=j+1;
 end
 clearvars no_of_points
end

% removing duplicates
possible_edges=unique(possible_edges,'rows');

% checking for commonality between possible and existing edges in graph
A=intersect(table2array(G3.Edges),possible_edges,'rows');

% removing the common edges, the remainder can be added to new graph
possible_edges=setdiff(possible_edges,A,'rows');

for i=1:numel(possible_edges(:,1))
  chain  = [XY_new(possible_edges(i,1),1:2) XY_new(possible_edges(i,2),1:2) ];
  possible_edges(i,3) = Lengths2D(chain);
  clearvars chain
end


end

