function [sub_G,a_sub_G,xy] = create_subgraph(path_edges,G3,XY3)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

nodes = unique(path_edges,'stable');
sub_G = subgraph(G3,nodes);
xy = XY3(nodes,:);
a_sub_G=adjacency(sub_G);

end

