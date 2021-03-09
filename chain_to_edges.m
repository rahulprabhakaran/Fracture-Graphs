function [edges] = chain_to_edges(chain)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
[r,c] = size(chain);
if c==1 && r>1 % check for col vector
    edges =chain;
else    
    edges = chain';
end
edges(numel(edges)-1,2) = edges(numel(edges),1);
edges(end,:) = [];
edges(1:numel(edges(:,1))-1,2) = edges(2:numel(edges(:,1)),1); 
end

