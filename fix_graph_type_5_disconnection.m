function [G6] = fix_graph_type_5_disconnection(G5,possible_edges_I_M_kinked)

% this function appends the possible I-M edges to the existing graph edges

% creating copy of graph 
G6=G5;

% adding edges
for i=1:numel(possible_edges_I_M_kinked(:,1))
  disp(i)  
  G6 = addedge(G6,possible_edges_I_M_kinked(i,1),possible_edges_I_M_kinked(i,2));   
end


end

