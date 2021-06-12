function [last_node] = nieghbors_degree_not_2(G,current_node,preceding_node)
% this function reads a graph and two nodes i.e., penultimate node in a
% chain and the preceeding node. The last node in the chain is returned 

neigh_nodes = [neighbors(G,current_node) degree(G,neighbors(G,current_node))];
last_node = setdiff(neigh_nodes(find(neigh_nodes(:,2)~=2),1),preceding_node);

end

