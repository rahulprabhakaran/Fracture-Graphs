function [out] = check_if_line_graph(sub_G3)
% this function reads a sub-graph and checks whether it is a line
% graph or not. The checks are if the start and end between centrality
% is zero, and if the minimum spanning tree is isomorphic with the
% graph itself

bw_cent = centrality(sub_G3,'betweenness');
T=minspantree(sub_G3);

if  (isisomorphic(T,sub_G3)==1) && (numel(bw_cent(bw_cent==0,1))==2)
   out = 1;
else
   out = 0; 
end    
   
end

