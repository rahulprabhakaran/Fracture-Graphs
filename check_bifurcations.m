function [out] = check_bifurcations(chain,G3)
% checks if a path has any bifurcations
% G3.Nodes.Name = [];
sub_G3=subgraph(G3,chain);

if check_if_line_graph(sub_G3)==1
   out = 0;
elseif check_if_fork_graph(sub_G3)==1
   out = 1;
else
   out = 2;
end   
    

end

