function [chains] = find_chain(Nodes_to_find_chains,G4)
% this function read a vector of graph nodes with degree 1 and finds the
% chains that start from source degree-1 node, traverse threough degree-2
% nodes and terminate in a non-degree 2 node. The chains corresponding to
% each degree-1 node is stored in a cell array called chains and returned

chains={};
    for j=1:length(Nodes_to_find_chains)
      nodeID = Nodes_to_find_chains(j,1);
      neigh_node = neighbors(G4,nodeID);
      edge_path = [nodeID;neigh_node];
      prev_node = nodeID;     
      while degree(G4,neigh_node)==2
         next_neigh = setdiff(neighbors(G4,neigh_node),prev_node); 
         edge_path = [edge_path;next_neigh];
         prev_node = neigh_node;
         neigh_node = next_neigh;
      end
      chains{j,1}=edge_path;
%       clearvars nodeID neigh_node edge_path prev_node next_neigh
    end  
     
end

