function [G2,XY_new] = fix_graph_type_8_disconnection(G1,small_I_I_I_trielements_AR,XY)
% this function is to rectify, type 8 disconnections, i.e, I-I-I triple
% point type disconnections. The input graph is the input argument and the
% small_I_I_I_trielements. These are the trielements  which are sorted out
% by smallest area. The second criteria is aspect ratio.

trielement_centroids = zeros(numel(small_I_I_I_trielements_AR(:,1)),1);
for i=1:numel(small_I_I_I_trielements_AR(:,1))
  disp(i)  
  trielement_centroids(i,1:2)=centroid(XY(small_I_I_I_trielements_AR(i,1:3),1), XY(small_I_I_I_trielements_AR(i,1:3),2));
end

G2=G1;
XY_new = XY;


for i=1:numel(small_I_I_I_trielements_AR(:,1))
   disp(i)
%   chain1 = Combined_Edges_cleaned_duplicate_removed{small_I_I_I_trielements_AR(i,7),1};
%   chain2 = Combined_Edges_cleaned_duplicate_removed{small_I_I_I_trielements_AR(i,8),1};
%   chain3 = Combined_Edges_cleaned_duplicate_removed{small_I_I_I_trielements_AR(i,9),1};
  G2 = addnode(G2,1);
  G2 = addedge(G2,small_I_I_I_trielements_AR(i,1),numnodes(G2));
  G2 = addedge(G2,small_I_I_I_trielements_AR(i,2),numnodes(G2));
  G2 = addedge(G2,small_I_I_I_trielements_AR(i,3),numnodes(G2));
%   numedges(G3)
%   Combined_Edges_cleaned_duplicate_removed{small_I_I_I_trielements_AR(i,7),1}= ...
%       [Combined_Edges_cleaned_duplicate_removed{small_I_I_I_trielements_AR(i,7),1};numnodes(G2)];
%   Combined_Edges_cleaned_duplicate_removed{small_I_I_I_trielements_AR(i,8),1}= ...
%       [Combined_Edges_cleaned_duplicate_removed{small_I_I_I_trielements_AR(i,8),1};numnodes(G2)];
%   Combined_Edges_cleaned_duplicate_removed{small_I_I_I_trielements_AR(i,9),1}= ...
%       [Combined_Edges_cleaned_duplicate_removed{small_I_I_I_trielements_AR(i,9),1};numnodes(G2)];
  XY_new = [XY_new;trielement_centroids(i,:)];
%   clearvars chain1 chain2 chain3
end    
    
end

