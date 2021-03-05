function [G3,XY3,removed_V_nodes] = straighten_graph(G2,XY,Combined_Edges)
% this function reads a graph whose edges are curved and have intermediate 
% degree 2 nodes, and returns a graph without these intermediary nodes

G3 = graph();
XY3 = XY;

D2 = degree(G2);
[V_nodes_idx,~] = find(D2(:,1)==2);

for i=1:length(Combined_Edges)
   disp(i)
   Chain=str2num(Combined_Edges{i,4})';
   node_s = Chain(1,1);
   node_t = Chain(end,1);
%    Chain(1:numel(Chain(:,1))-1,2) = Chain(2:numel(Chain(:,1)),1); 
%    Chain(end,:)=[];

   G3 = addedge(G3,node_s,node_t);
%    G3 = rmedge(G3,Chain(:,1),Chain(:,2));
     removed_V_nodes{i,1} = Chain(D2(Chain)==2,1); 
   clearvars Chain
end   

G3 = rmnode(G3,V_nodes_idx);
XY3(V_nodes_idx,:)=[];

end

