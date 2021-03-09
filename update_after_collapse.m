function [G4,XY4] = update_after_collapse(G4,XY4,Edges_to_add,Edges_to_remove,Nodes_to_remove)
% update graph after performing collapsing a stepout

XY4(Nodes_to_remove(:,1),:)=[];
G4.Nodes.Name = string(1:1:G4.numnodes)';

G4 = rmedge(G4,Edges_to_remove(:,1));
G4 = addedge(G4,string(Edges_to_add(:,1)),string(Edges_to_add(:,2)));
G4 = rmnode(G4,string(Nodes_to_remove(:,1)));

end

