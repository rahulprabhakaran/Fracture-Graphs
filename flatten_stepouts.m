function [G5,XY5,PSST_NSE_flattened,PSST_SE] = flatten_stepouts(poss_stepout, orthogonal_range, strike_threshold, G4, XY4)
% this function reads a possible stepout edge, i.e., an edge which has
% degree 3 topology at both ends and performs a flatten operation that
% realigns the stepout to one pair of edges on either side of the 's'
% and 't' nodes that are favorably aligned; the procedure aids in
% improving connectivity during walk finding

% the flatten operation is performed only if the following criteria is 
% satisfied:
%  -> edges connecting to 's','t' nodes are 1.5 times stepout length
%  -> the stepout does not form a T-motif on either 's' or 't' side

stepout_multiple = 1.5;
orthogonal_min = 90 -  orthogonal_range;
orthogonal_max = 90 +  orthogonal_range;

% one20_deg_min = 120 - one20_deg_range;
% one20_deg_max = 120 + one20_deg_range;

%--------the loop identifies the rows that fit flatten constraints--------%
poss_stepout_s_t = zeros(numel(poss_stepout(:,1),1));
k=1;
for i=1:numel(poss_stepout(:,1)) 
%  disp(i)
 poss_stepout_s_t(i,1) = poss_stepout{i,1}(1,1);                            % s node 
 poss_stepout_s_t(i,2) = poss_stepout{i,1}(1,2);                            % t node

 % dividing angles for s node
 [poss_stepout_s_t(i,3),poss_stepout_s_t(i,4),~] = ...
               find_min_dividing_angle_deg3(poss_stepout_s_t(i,1),G4,XY4,poss_stepout_s_t(i,2));
 
 % dividing angles for t node           
 [poss_stepout_s_t(i,5),poss_stepout_s_t(i,6),~] = ...
               find_min_dividing_angle_deg3(poss_stepout_s_t(i,2),G4,XY4,poss_stepout_s_t(i,1));
 
 % stepout length           
 poss_stepout_s_t(i,7) = compute_chain_length(poss_stepout{i,1}, XY4);
 
 % neighbours of s node
 poss_stepout_s_t(i,8:9) = setdiff(neighbors(G4,poss_stepout_s_t(i,1)),poss_stepout_s_t(i,2))';
 
 % s node neighbour 1 edge length
 poss_stepout_s_t(i,10) = compute_chain_length([poss_stepout_s_t(i,1) poss_stepout_s_t(i,8)], XY4);
 
 % s node neighbour 2 edge length
 poss_stepout_s_t(i,11) = compute_chain_length([poss_stepout_s_t(i,1) poss_stepout_s_t(i,9)], XY4);
 
 % neighbours of t node
 poss_stepout_s_t(i,12:13) = setdiff(neighbors(G4,poss_stepout_s_t(i,2)),poss_stepout_s_t(i,1))';
 
 % t node neighbour 1 edge length
 poss_stepout_s_t(i,14) = compute_chain_length([poss_stepout_s_t(i,2) poss_stepout_s_t(i,12)], XY4);
 

 % t node neighbour 2 edge length
 poss_stepout_s_t(i,15) = compute_chain_length([poss_stepout_s_t(i,2) poss_stepout_s_t(i,13)], XY4);
 
    %isinrange_atleast(poss_stepout_s_t(i,3:6),one20_deg_min, one20_deg_max)==0
if isinrange_atleast(poss_stepout_s_t(i,3:6),orthogonal_min, orthogonal_max)==0 && ...
        poss_stepout_s_t(i,14) > stepout_multiple*poss_stepout_s_t(i,7) && ...
        poss_stepout_s_t(i,15) > stepout_multiple*poss_stepout_s_t(i,7)
    
  rows_to_flatten(k,1) = i;
  k=k+1;
end
         
end

% considering only rows that are to be flattened
poss_stepout_s_t = poss_stepout_s_t(rows_to_flatten,:);
poss_stepout = poss_stepout(rows_to_flatten,:);

%%
PSST = find_LS_SR_stepout_edges(poss_stepout,G4);
edge_counter = [PSST(:,1:2); [PSST(:,1) PSST(:,3)];...
    [PSST(:,1) PSST(:,4)]; [PSST(:,2) PSST(:,5)];...
    [PSST(:,2) PSST(:,6)]];
edge_counter(:,3) = 1;
edge_counter(:,4) = findedge(G4,edge_counter(:,1),edge_counter(:,2));
G4.Nodes.Name = string(1:1:G4.numnodes)';

edge_list =    [  [findedge(G4,PSST(:,1),PSST(:,2))]...  % S-T  edge index
                  [findedge(G4,PSST(:,1),PSST(:,3))]...  % L1-S edge index
                  [findedge(G4,PSST(:,1),PSST(:,4))]...  % L2-S edge index
                  [findedge(G4,PSST(:,2),PSST(:,5))]...  % R1-S edge index
                  [findedge(G4,PSST(:,2),PSST(:,6))]];   % R2-S edge index

% checking for shared edges between S-T and L1,L2, S-T and R1,R2
for i=1:numel(edge_list(:,1))
   % all stepouts that share L-S edges of other stepouts
   share_ST_L1{i,1} =  find(edge_list(:,2)==edge_list(i,1));
   if ~isempty(share_ST_L1{i,1})
    share_ST_L1{i,1} = i;
   end
   share_ST_L2{i,1} =  find(edge_list(:,3)==edge_list(i,1));
   if ~isempty(share_ST_L2{i,1})
    share_ST_L2{i,1} = i;
   end  
   % all stepouts that share R-S edges of other stepouts
   share_ST_R1{i,1} =  find(edge_list(:,4)==edge_list(i,1));
   if ~isempty(share_ST_R1{i,1})
    share_ST_R1{i,1} = i;
   end
   share_ST_R2{i,1} =  find(edge_list(:,5)==edge_list(i,1));
   if ~isempty(share_ST_R2{i,1})
    share_ST_R2{i,1} = i;
   end
   %----------------------------------------------------------------------%
   % all L1_S edges sharing edges with L2_S, R1_S, R2_S
   share_L1S_L2S{i,1} =  find(edge_list(:,3)==edge_list(i,2));
   if ~isempty(share_L1S_L2S{i,1})
    share_L1S_L2S{i,1} = i;
   end
   share_L1S_R1S{i,1} =  find(edge_list(:,4)==edge_list(i,2));
   if ~isempty(share_L1S_R1S{i,1})
    share_L1S_R1S{i,1} = i;
   end
   share_L1S_R2S{i,1} =  find(edge_list(:,5)==edge_list(i,2));
   if ~isempty(share_L1S_R2S{i,1})
    share_L1S_R2S{i,1} = i;
   end
   %----------------------------------------------------------------------%
   share_L2S_R1S{i,1} =  find(edge_list(:,4)==edge_list(i,3));
   if ~isempty(share_L2S_R1S{i,1})
    share_L2S_R1S{i,1} = i;
   end
   share_L2S_R2S{i,1} =  find(edge_list(:,5)==edge_list(i,3));
   if ~isempty(share_L2S_R2S{i,1})
    share_L2S_R2S{i,1} = i;
   end
   %----------------------------------------------------------------------%
   share_R1S_R2S{i,1} =  find(edge_list(:,5)==edge_list(i,4));
   if ~isempty(share_R1S_R2S{i,1})
    share_R1S_R2S{i,1} = i;
   end
end

% extracing rows from edgelist
share_ST_L1 = cell2mat(share_ST_L1);
share_ST_L2 = cell2mat(share_ST_L2);
share_ST_R1 = cell2mat(share_ST_R1);
share_ST_R2 = cell2mat(share_ST_R2);
share_L1S_L2S = cell2mat(share_L1S_L2S);
share_L1S_R1S = cell2mat(share_L1S_R1S);
share_L1S_R2S = cell2mat(share_L1S_R2S);
share_L2S_R1S = cell2mat(share_L2S_R1S);
share_L2S_R2S = cell2mat(share_L2S_R2S);
share_R1S_R2S = cell2mat(share_R1S_R2S);

% rows in poss_stepout that have any type of shared edges (ST, LS or RS)
all_shared_edges_idx = unique([share_ST_L1;share_ST_L2;share_ST_R1;share_ST_R2; ...
   share_L1S_L2S;share_L1S_R1S;share_L1S_R2S;share_L2S_R1S;share_L2S_R2S;share_R1S_R2S]);


% all stepouts that share edges with L-S, R-S edges of other stepouts   
%all_shared_edges_idx = unique([share_ST_L1;share_ST_L2;share_ST_R1;share_ST_R2]);


% rows in poss_stepout that have no shared edges
no_shared_edges = setdiff(1:1:numel(poss_stepout(:,1)),all_shared_edges_idx)';

% stepout lists split into those with and without shared edges
% the ones with shared edges may need dynamic graph manipulation within the loop
% the ones without shared edges can be manipulated in one vectorized operation
PSST_SE = poss_stepout(all_shared_edges_idx,1);
PSST_NSE = poss_stepout(no_shared_edges,1);
PSST_S_T_SE = poss_stepout_s_t(all_shared_edges_idx,:);
PSST_S_T_NSE = poss_stepout_s_t(no_shared_edges,:);

clearvars share_L1S_L2S share_L1S_R1S share_L1S_R2S share_L2S_R1S share_L2S_R2S share_R1S_R2S share_ST_L1 share_ST_L2 share_ST_R1 share_ST_R2
%%
%-finding combinations of edges on both sides of stepout, no shared edges-%
%
k=1;
for i=1:numel(PSST_S_T_NSE(:,1))
%    disp(i) 
   stepout_strike =  compute_chain_strike2(PSST_S_T_NSE(i,1:2), XY4);
   % all possible combos of edges on left and right & including stepout
   combs = [[repmat(PSST_S_T_NSE(i,8),[2 1]);repmat(PSST_S_T_NSE(i,9),[2 1])] ...
       repmat(PSST_S_T_NSE(i,1:2),[4 1]) repmat(PSST_S_T_NSE(i,12:13)',[2 1])];
   
   combs_idx = [1    3; ...   % L1-S  R1-S  combination index structure      
                1    4; ...   % L1-S  R2-S
                2    3; ...   % L2-S  R1-S
                2    4];      % L2-S  R2-S
   
   % combos on left and right L1-S, L2-S, R1-S, R2-S for walk checks
   combs2 = { [PSST_S_T_NSE(i,8) PSST_S_T_NSE(i,1)]; ...
                [PSST_S_T_NSE(i,9) PSST_S_T_NSE(i,1)]; ...
                 [PSST_S_T_NSE(i,2) PSST_S_T_NSE(i,12)]; ...
                  [PSST_S_T_NSE(i,2) PSST_S_T_NSE(i,13)] };
              
   combs2 = { combs(1,1:2); combs(3,1:2); combs(3,3:4); combs(4,3:4)};
        
   % finding walks for each of L1-S, L2-S, R1-S, R2-S, length of each walk     
   for m=1:4
      combs2{m,2} = find_walk([combs2{m,1}],strike_threshold,G4,XY4);            
      combs2{m,3} = numel(combs2{m,2});
   end
   
   % finding the row index of the longest walk on left
   if ~isempty(find( cell2mat(combs2(1:2,3))>2))
     [~,r_max_L_S_walk] =  max(cell2mat(combs2(1:2,3)));  % row of longest walk
   else
     r_max_L_S_walk = [];  
   end    
   
   % finding the row index of the longest walk on right
   if ~isempty(find( cell2mat(combs2(3:4,3))>2))
     [~,r_max_R_S_walk] =  max(cell2mat(combs2(3:4,3)));  % row of longest walk
     r_max_R_S_walk = r_max_R_S_walk+2;                   % adding 2 to get cumulative row number
   else
     r_max_R_S_walk = [];
   end    
       
   % difference in strike for each combination of left and right
   diff = abs(compute_strike(combs(:,1:2),XY4) - compute_strike(combs(:,3:4),XY4));
   
   % condition when connected edge combo strikes are within threshold
   if ~isempty(find(diff<strike_threshold))
%      if numel(find(diff<strike_threshold))<3  % choosing only those with 1 or 2 possibilities
%       r_conns{k,1}=numel(find(diff<strike_threshold));   % number of possible combos
%       r_conns{k,2}=combs(find(diff(:,1)==min(diff)),:);  % choosing one with min difference
%       r_conns{k,3}=i;                                    % row numbers to extract
%       k=k+1;
%      end 
     %------------------------------------------------------------------% 
     if numel(find(diff<strike_threshold))<3 % choosing only those with 1 or 2 possibilities
       % when there is a definite longest walk on both sides
       if ~isempty(r_max_L_S_walk) && ~isempty(r_max_R_S_walk)
          r_conns{k,1}=numel(find(diff<strike_threshold));   % number of possible combos
          r_conns{k,2} = [combs2{r_max_L_S_walk,1} combs2{r_max_R_S_walk,1}];
          r_conns{k,3}=i;                                    % row numbers to extract
          k=k+1;  
       end
       % no walk on right, walk on left
       if ~isempty(r_max_L_S_walk) && isempty(r_max_R_S_walk)
          r_conns{k,1}=numel(find(diff<strike_threshold));   % number of possible combos
          r_tmp = find(combs_idx(:,1)==r_max_L_S_walk);      % rows in comb_idx selected for diff
          r_conns{k,2} = combs(find(diff(r_tmp,1)==min(diff(r_tmp,1))),:);
          r_conns{k,3}=i;                                    % row numbers to extract
          k=k+1;  
       end 
       % no walk on left, walk on right
       if isempty(r_max_L_S_walk) && ~isempty(r_max_R_S_walk)
          r_conns{k,1}=numel(find(diff<strike_threshold));   % number of possible combos 
          r_tmp = find(combs_idx(:,2)==r_max_R_S_walk);      % rows in comb_idx selected for diff
          r_conns{k,2} = combs(find(diff(r_tmp,1)==min(diff(r_tmp,1))),:); 
          r_conns{k,3}=i;                                    % row numbers to extract
          k=k+1;  
       end
       % no walk either on left or right
       if isempty(r_max_L_S_walk) && isempty(r_max_R_S_walk)
          r_conns{k,1}=numel(find(diff<strike_threshold));   % number of possible combos
          r_conns{k,2}=combs(find(diff(:,1)==min(diff)),:);  % choosing combo with min difference
          r_conns{k,3}=i;                                    % row numbers to extract
          k=k+1;  
       end    
     end    
   end 
end

PSST_NSE_combinable = PSST_NSE(cell2mat(r_conns(:,3)),:);
% r_4_conns = find(cell2mat(r_conns(:,1))==4);
% r_3_conns = find(cell2mat(r_conns(:,1))==3);
r_2_conns = find(cell2mat(r_conns(:,1))==2);
% r_1_conns = find(cell2mat(r_conns(:,1))==1);
% poss_stepout_4_conns = poss_stepout(cell2mat(r_conns(r_4_conns,3)),:);
% poss_stepout_3_conns = poss_stepout(cell2mat(r_conns(r_3_conns,3)),:);
% poss_stepout_1_conns = poss_stepout(cell2mat(r_conns(r_1_conns,3)),:);
PSST_NSE_2_conns = PSST_NSE(cell2mat(r_conns(r_2_conns,3)),:);

%%
%--------flattening the stepouts------------------------------------------%
% for each stepout, 3 edges added & removed, 1 node added & 1 removed; this
% loop makes a list of all these edges and nodes
Edges_to_add = []; Edges_to_remove = []; Nodes_to_add =[]; Nodes_to_remove =[];
k=1;
flattened_rows=[];
for i=1:numel(PSST_NSE_combinable(:,1))
%    tic 
%    disp(i) 
   edge_comb = r_conns{i,2};

   % converting 2nd edge diverging out of left of stepout to a line &
   % finding projection point; G3 is used as reference because graph G4 and
   % spatial referencing matrix is dynamically varying as loop progresses
   edge_to_project = [setdiff(neighbors(G4,edge_comb(1,2)),[edge_comb(1,3) edge_comb(1,1)]) edge_comb(1,2)];
   line_to_project = edgeToLine([XY4( edge_to_project(1,1),1:2) XY4( edge_to_project(1,2),1:2)]);
   proj_point = intersectLineEdge(line_to_project, [XY4(edge_comb(1,1),1:2) XY4(edge_comb(1,3),1:2)]);
   
   % when the intersectLinEdge returns [NaN NaN], no projection is possible
   if isequal(isnan(proj_point),[0 0])
       Nodes_to_add = [Nodes_to_add; [G4.numnodes+k proj_point]];
       Edges_to_add = [Edges_to_add ;[edge_comb(1,1) G4.numnodes+k; edge_comb(1,3) G4.numnodes+k; edge_to_project(1,1) G4.numnodes+k]];
       Edges_to_remove = [Edges_to_remove; [edge_to_project(1,1) edge_comb(1,2); edge_comb(1,1) edge_comb(1,2); edge_comb(1,2) edge_comb(1,3)]];
       Nodes_to_remove = [Nodes_to_remove ; edge_comb(1,2)];
       flattened_rows(k,1)=i;
       k=k+1;
   end
   clearvars proj_point edge_comb edge_to_project line_to_project
%    toc
end   
    
PSST_NSE_flattened = PSST_NSE_combinable(flattened_rows,:);
clearvars r_conns r_2_conns flattened_rows

%% flattening for stepouts with shared edges

%-finding combinations of edges on both sides of stepout, no shared edges-%
%
% k=1;
% for i=1:numel(PSST_S_T_SE(:,1))
%    disp(i) 
%    stepout_strike =  compute_chain_strike2(PSST_S_T_SE(i,1:2), XY4);
%    % all possible combos of edges on left and right & including stepout
%    combs = [[repmat(PSST_S_T_SE(i,8),[2 1]);repmat(PSST_S_T_SE(i,9),[2 1])] ...
%        repmat(PSST_S_T_SE(i,1:2),[4 1]) repmat(PSST_S_T_SE(i,12:13)',[2 1])];
%    
%    combs_idx = [1    3; ...   % L1-S  R1-S  combination index structure      
%                 1    4; ...   % L1-S  R2-S
%                 2    3; ...   % L2-S  R1-S
%                 2    4];      % L2-S  R2-S
%    
%    % combos on left and right L1-S, L2-S, R1-S, R2-S for walk checks
%    combs2 = { [PSST_S_T_SE(i,8) PSST_S_T_SE(i,1)]; ...
%                 [PSST_S_T_SE(i,9) PSST_S_T_SE(i,1)]; ...
%                  [PSST_S_T_SE(i,2) PSST_S_T_SE(i,12)]; ...
%                   [PSST_S_T_SE(i,2) PSST_S_T_SE(i,13)] };
%               
%    combs2 = { combs(1,1:2); combs(3,1:2); combs(3,3:4); combs(4,3:4)};
%         
%    % finding walks for each of L1-S, L2-S, R1-S, R2-S, length of each walk     
%    for m=1:4
%       combs2{m,2} = find_walk([combs2{m,1}],strike_threshold,G4,XY4);            
%       combs2{m,3} = numel(combs2{m,2});
%    end
%    
%    % finding the row index of the longest walk on left
%    if ~isempty(find( cell2mat(combs2(1:2,3))>2))
%      [~,r_max_L_S_walk] =  max(cell2mat(combs2(1:2,3)));  % row of longest walk
%    else
%      r_max_L_S_walk = [];  
%    end    
%    
%    % finding the row index of the longest walk on right
%    if ~isempty(find( cell2mat(combs2(3:4,3))>2))
%      [~,r_max_R_S_walk] =  max(cell2mat(combs2(3:4,3)));  % row of longest walk
%      r_max_R_S_walk = r_max_R_S_walk+2;                   % adding 2 to get cumulative row number
%    else
%      r_max_R_S_walk = [];
%    end    
%        
%    % difference in strike for each combination of left and right
%    diff = abs(compute_strike(combs(:,1:2),XY4) - compute_strike(combs(:,3:4),XY4));
%    
%    % condition when connected edge combo strikes are within threshold
%    if ~isempty(find(diff<strike_threshold))
% %      if numel(find(diff<strike_threshold))<3  % choosing only those with 1 or 2 possibilities
% %       r_conns{k,1}=numel(find(diff<strike_threshold));   % number of possible combos
% %       r_conns{k,2}=combs(find(diff(:,1)==min(diff)),:);  % choosing one with min difference
% %       r_conns{k,3}=i;                                    % row numbers to extract
% %       k=k+1;
% %      end 
%      %------------------------------------------------------------------% 
%      if numel(find(diff<strike_threshold))<3 % choosing only those with 1 or 2 possibilities
%        % when there is a definite longest walk on both sides
%        if ~isempty(r_max_L_S_walk) && ~isempty(r_max_R_S_walk)
%           r_conns{k,1}=numel(find(diff<strike_threshold));   % number of possible combos
%           r_conns{k,2} = [combs2{r_max_L_S_walk,1} combs2{r_max_R_S_walk,1}];
%           r_conns{k,3}=i;                                    % row numbers to extract
%           k=k+1;  
%        end
%        % no walk on right, walk on left
%        if ~isempty(r_max_L_S_walk) && isempty(r_max_R_S_walk)
%           r_conns{k,1}=numel(find(diff<strike_threshold));   % number of possible combos
%           r_tmp = find(combs_idx(:,1)==r_max_L_S_walk);      % rows in comb_idx selected for diff
%           r_conns{k,2} = combs(find(diff(r_tmp,1)==min(diff(r_tmp,1))),:);
%           r_conns{k,3}=i;                                    % row numbers to extract
%           k=k+1;  
%        end 
%        % no walk on left, walk on right
%        if isempty(r_max_L_S_walk) && ~isempty(r_max_R_S_walk)
%           r_conns{k,1}=numel(find(diff<strike_threshold));   % number of possible combos 
%           r_tmp = find(combs_idx(:,2)==r_max_R_S_walk);      % rows in comb_idx selected for diff
%           r_conns{k,2} = combs(find(diff(r_tmp,1)==min(diff(r_tmp,1))),:); 
%           r_conns{k,3}=i;                                    % row numbers to extract
%           k=k+1;  
%        end
%        % no walk either on left or right
%        if isempty(r_max_L_S_walk) && isempty(r_max_R_S_walk)
%           r_conns{k,1}=numel(find(diff<strike_threshold));   % number of possible combos
%           r_conns{k,2}=combs(find(diff(:,1)==min(diff)),:);  % choosing combo with min difference
%           r_conns{k,3}=i;                                    % row numbers to extract
%           k=k+1;  
%        end    
%      end    
%    end 
% end
% 
% PSST_SE_combinable = PSST_SE(cell2mat(r_conns(:,3)),:);
% % r_4_conns = find(cell2mat(r_conns(:,1))==4);
% % r_3_conns = find(cell2mat(r_conns(:,1))==3);
% r_2_conns = find(cell2mat(r_conns(:,1))==2);
% % r_1_conns = find(cell2mat(r_conns(:,1))==1);
% % poss_stepout_4_conns = poss_stepout(cell2mat(r_conns(r_4_conns,3)),:);
% % poss_stepout_3_conns = poss_stepout(cell2mat(r_conns(r_3_conns,3)),:);
% % poss_stepout_1_conns = poss_stepout(cell2mat(r_conns(r_1_conns,3)),:);
% PSST_SE_2_conns = PSST_SE(cell2mat(r_conns(r_2_conns,3)),:);
% 
% %%
% %--------flattening the stepouts with shared edges------------------------%
% % for each stepout, 3 edges added & removed, 1 node added & 1 removed; this
% % loop makes a list of all these edges and nodes
% 
% k=1;
% flattened_rows=[];
% for i=1:numel(PSST_SE_combinable(:,1))
%    tic 
%    disp(i) 
%    edge_comb = r_conns{i,2};
% 
%    % converting 2nd edge diverging out of left of stepout to a line &
%    % finding projection point; G3 is used as reference because graph G4 and
%    % spatial referencing matrix is dynamically varying as loop progresses
%    edge_to_project = [setdiff(neighbors(G4,edge_comb(1,2)),[edge_comb(1,3) edge_comb(1,1)]) edge_comb(1,2)];
%    line_to_project = edgeToLine([XY4( edge_to_project(1,1),1:2) XY4( edge_to_project(1,2),1:2)]);
%    proj_point = intersectLineEdge(line_to_project, [XY4(edge_comb(1,1),1:2) XY4(edge_comb(1,3),1:2)]);
%    
%    % when the intersectLinEdge returns [NaN NaN], no projection is possible
%    if isequal(isnan(proj_point),[0 0])
%        Nodes_to_add = [Nodes_to_add; [G4.numnodes+k proj_point]];
%        Edges_to_add = [Edges_to_add ;[edge_comb(1,1) G4.numnodes+k; edge_comb(1,3) G4.numnodes+k; edge_to_project(1,1) G4.numnodes+k]];
%        Edges_to_remove = [Edges_to_remove; [edge_to_project(1,1) edge_comb(1,2); edge_comb(1,1) edge_comb(1,2); edge_comb(1,2) edge_comb(1,3)]];
%        Nodes_to_remove = [Nodes_to_remove ; edge_comb(1,2)];
%        flattened_rows(k,1)=i;
%        k=k+1;
%    end
%    clearvars proj_point edge_comb edge_to_project line_to_project
%    toc
% end   
%     
% PSST_SE_flattened = PSST_SE_combinable(flattened_rows,:);

%% manipulating graph for the stepouts
G4.Nodes.Name = string(1:1:G4.numnodes)'; %% assigning node names
G5=G4;      % duplicating the graph
XY5 = XY4;  % duplicating spatial referencing matrix

XY5 = [XY5; Nodes_to_add(:,2:3)];
XY5(Nodes_to_remove(:,1),:)=[];

Edges_to_add = unique(Edges_to_add,'rows');
Edges_to_remove = unique(Edges_to_remove,'rows');

% Edges_to_add = sort(Edges_to_add,2);
% Edges_to_remove = sort(Edges_to_remove,2);

G5 = addnode(G5,string(Nodes_to_add(:,1)));
G5 = addedge(G5,string(Edges_to_add(:,1)),string(Edges_to_add(:,2)));
G5 = rmedge(G5,string(Edges_to_remove(:,1)),string(Edges_to_remove(:,2)));
G5 = rmnode(G5,string(Nodes_to_remove(:,1)));

end

