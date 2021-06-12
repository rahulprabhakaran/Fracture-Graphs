function [succeeding_neigh_node_2_degree] = nieghbors_degree_2(G,current_node,preceding_node)
% this function reads a graph and an origin node. The neighboring nodes and
% their degrees are calculated. Only the nodes with a degree of 2 are 
% returned
neigh_nodes = [neighbors(G,current_node) degree(G,neighbors(G,current_node))];

neigh_node_2_degree = neigh_nodes(find(neigh_nodes(:,2)==2),1);
succeeding_neigh_node_2_degree=setdiff(neigh_node_2_degree,preceding_node);

if isempty(succeeding_neigh_node_2_degree)==1
   last_node = setdiff(neigh_nodes(:,1),preceding_node);
%   succeeding_neigh_node_2_degree = last_node;
%    disp(strcat('End of chain. Last Node:   ',num2str(last_node))); 
end  

end

