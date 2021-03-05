function [poss_stepout] = find_stepout_motif(G3,XY3,stepout_length)
% this function returns a cell array containing all possible stepouts below
% a certain stepout length; function maybe used to 


Edges_G3=table2array(G3.Edges);
if iscell(Edges_G3)
    G3.Nodes.Name = [];
    Edges_G3 = table2array(G3.Edges);
end   

Edges_G3(:,3)=degree(G3,Edges_G3(:,1));
Edges_G3(:,4)=degree(G3,Edges_G3(:,2));
r_3_3 = find(Edges_G3(:,3)==3 & Edges_G3(:,4)==3);
poss_Edges = Edges_G3(r_3_3,:);

poss_Edges(:,5)=edgeLength([XY3(poss_Edges(:,1),1) XY3(poss_Edges(:,1),2) XY3(poss_Edges(:,2),1) XY3(poss_Edges(:,2),2)]);
r_stpout_length = find(poss_Edges(:,5)<stepout_length);

poss_Edges = poss_Edges(r_stpout_length,:);


poss_stepout = cell(0);
for i=1:numel(poss_Edges(:,1))
 poss_stepout{i,1} = poss_Edges(i,1:2);
end


end

