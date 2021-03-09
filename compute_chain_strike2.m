function [length_weighted_strike,Q_info] = compute_chain_strike2(Chain_edges,XY3)
% this function computes length weighted strike for a chain of edges
% if all edges within the chain have unit vectors in Q1-Q3 direction (or
% Q2-Q4 direction) then choosing the minimum is correct. If some edges are
% in Q2-Q4 direction and some in Q1-Q3 then the calculation is inaccurate
% such a situation happens when we have near-horizontal edges.

% Chain_edges = Near_Horz_walks{27,1};
Chain_edges = Chain_edges';
Chain_edges(numel(Chain_edges)-1,2) = Chain_edges(numel(Chain_edges),1);
Chain_edges(end,:) = [];
Chain_edges(1:numel(Chain_edges(:,1))-1,2) = Chain_edges(2:numel(Chain_edges(:,1)),1); 
[Chain_edges(:,3),quadrants] = compute_strike(Chain_edges(:,1:2), XY3);

for i=1:numel(Chain_edges(:,1))    
 Chain_edges(i,4) = distancePoints( [ XY3(Chain_edges(i,1),1) XY3(Chain_edges(i,1),2)], ...
     [XY3(Chain_edges(i,2),1) XY3(Chain_edges(i,2),2)]);
end
cum_length = sum(Chain_edges(:,4));


% check if there are edges in chain that are close to horizontal
if ~isempty(find(Chain_edges(:,3)>150)) && ~isempty(find(Chain_edges(:,3)<30))
  trans = createRotation(deg2rad(60));  % transformation matrix for rotation
  rotated_Chain = transformEdge([XY3(Chain_edges(:,1),:) XY3(Chain_edges(:,2),:)],trans);
  rotated_Chain(:,5) = min(rad2deg(edgeAngle(rotated_Chain(:,1:4))),rad2deg(edgeAngle(reverseEdge(rotated_Chain(:,1:4)))));
  Chain_edges(:,5) = rotated_Chain(:,5).*Chain_edges(:,4);  
  length_weighted_strike = sum(rotated_Chain(:,5).*Chain_edges(:,4)/cum_length) - 60;
else    
  Chain_edges(:,5) = Chain_edges(:,3).*Chain_edges(:,4);
  length_weighted_strike = sum(Chain_edges(:,3).*Chain_edges(:,4)/cum_length);
end  

% in near horizontal chains, there can be some edges which trend in Q2 and 
% Q4; these edges can skew the chain strike calculation, and they are 
% converted into Q1 by subtracting by 180 degrees
r_Q_24 = 0;
r_Q_13 = 0;

 if numel(unique(quadrants))>1    %all(quadrants ~= quadrants(1))   
     r_Q_13 = find(quadrants(:,1)==13);
     r_Q_24 = find(quadrants(:,1)==24);
     
%      if numel(r_Q_13) > numel(r_Q_24) || mean(Chain_edges(r_Q_24,3)) > 160 % 
%          Chain_edges(r_Q_24,3) = 180 - Chain_edges(r_Q_24,3);
%      end     % this condition is set when few near horz edges are in Q24
% 
%      
%      if numel(r_Q_13) < numel(r_Q_24) || mean(Chain_edges(r_Q_13,3)) < 20 % 
%          Chain_edges(r_Q_13,3) = 180 - Chain_edges(r_Q_13,3);
%      end     % this condition is set when few near horz edges are in Q13
%      
     Q_info = [numel(r_Q_13); numel(r_Q_24)];
 else
    Q_info = [numel(find(quadrants(:,1)==13)); numel(find(quadrants(:,1)==24))]; 
end  

end

