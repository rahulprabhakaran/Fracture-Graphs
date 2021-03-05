function [possible_edges_I_M_kinked] = find_I_M_possible_edges(dt_5,kinked_M_nodes,I_nodes_idx,I_nodes,XY5,G5)

j=1;
for i=1:numel(I_nodes_idx)
%     disp(i)    
    % all tri-elements adjacent to the kinked_M_node, vertex attachments
    vertex_attachments = cell2mat(vertexAttachments(dt_5,I_nodes_idx(i,:)))' ;
    
    % extracting points of all the adjacent vertices, reshaping into a
    % column vector, calculating graph degree applying unique operation, and then e
 
    P_vertex_attachments = setdiff(XY5(unique(reshape(dt_5.ConnectivityList(vertex_attachments,:),...
        [numel(dt_5.ConnectivityList(vertex_attachments,:)) 1])),1:2),I_nodes(i,1:2),'rows');
    
    P_vertex_attachments(:,3) =  setdiff(unique(reshape(dt_5.ConnectivityList(vertex_attachments,:),...
        [numel(dt_5.ConnectivityList(vertex_attachments,:)),1])), I_nodes_idx(i,:),'rows');
    
    P_vertex_attachments(:,4) =  degree(G5,P_vertex_attachments(:,3));   
    
    % checking for commonality between I_nodes and kinked_M_nodes
    % check is performed if kinked M-node lies on an I-node chain using
    % intersect operation if so, then that kinked M-node must not be considered
    % as a possible edge and is removed using setdiff operation
    if ~isempty(intersect(kinked_M_nodes,P_vertex_attachments(:,3)))
       edge_M_vertices{j,1}=intersect(kinked_M_nodes,P_vertex_attachments(:,3));
       edge_M_vertices{j,2}=I_nodes_idx(i,:);
       if ~isempty(intersect(cell2mat(find_chain(edge_M_vertices{j,2},G5)),edge_M_vertices{j,1}))
         edge_M_vertices{j,1}=setdiff(edge_M_vertices{j,1},intersect(cell2mat(find_chain(edge_M_vertices{j,2},G5)),edge_M_vertices{j,1}));
       end   
       j=j+1;
    end   
     
%     edge_M_vertices = edge_M_vertices(~cellfun(@isempty, edge_M_vertices(:,1)), :);

    clearvars vertex_attachments P_vertex_attachments
end

edge_M_vertices = edge_M_vertices(~cellfun(@isempty, edge_M_vertices(:,1)), :);

j=1;
for m=1:numel(edge_M_vertices(:,1))
   if numel(edge_M_vertices{m,1})>1 
     combinations(1:numel(edge_M_vertices{m,1}),1) = repmat(edge_M_vertices{m,2},numel(edge_M_vertices{m,1}),1);
     combinations(1:numel(edge_M_vertices{m,1}),2) = edge_M_vertices{m,1};
     possible_edges_I_M_kinked(j:j+numel(edge_M_vertices{m,1})-1,:) = combinations;
     j=j+numel(edge_M_vertices{m,1});
   else
     possible_edges_I_M_kinked(j,1)= edge_M_vertices{m,2};
     possible_edges_I_M_kinked(j,2)= edge_M_vertices{m,1};
     j=j+1;
   end    
   clearvars combinations 
end

possible_edges_I_M_kinked(:,3) = Lengths2D( [ XY5(possible_edges_I_M_kinked(:,1),1:2)  XY5(possible_edges_I_M_kinked(:,2),1:2) ]);
[a,ia,iaa] = unique(possible_edges_I_M_kinked(:,1),'rows');
possible_edges_I_M_kinked(:,4)=iaa;
multiple_edges_idx = find(hist(iaa,unique(iaa))>1)';

for m=1:numel(iaa)
   if find(ismember(multiple_edges_idx,m))~=0
       [b,~] = find(possible_edges_I_M_kinked(:,4)==m);
       [c,~] = find(possible_edges_I_M_kinked(b,3)==min(possible_edges_I_M_kinked(b,3)));
       b(c,:)=[];   
       possible_edges_I_M_kinked(b,:)=0;
   end  
   clearvars b c
end    

possible_edges_I_M_kinked=possible_edges_I_M_kinked(possible_edges_I_M_kinked(:,1)~=0,:);

end

