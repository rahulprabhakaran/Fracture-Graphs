function [B_s_t] = find_LS_SR_stepout_edges(B,G3)


for i=1:numel(B(:,1))
 B_s_t(i,1) = B{i,1}(1,1);                            % s node 
 B_s_t(i,2) = B{i,1}(1,2);                            % t node
    % neighbours of s node
 B_s_t(i,3:4) = setdiff(neighbors(G3,B_s_t(i,1)),B_s_t(i,2))';
    % neighbours of t node
 B_s_t(i,5:6) = setdiff(neighbors(G3,B_s_t(i,2)),B_s_t(i,1))';   
end


end

