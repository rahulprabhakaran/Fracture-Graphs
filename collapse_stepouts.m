function [G4,XY4,uncollapsable_stepout,PSSTSE,poss_stepout_no_shared_edges] = collapse_stepouts(poss_stepout,G3,XY3)
% outputs: A -> uncollapsable_stepout
%          B -> poss_stepout_shared_edges
%          C -> poss_stepout_no_shared_edges

% stepouts that are very small, i.e., around 2 pixel lengths
% can be collapsed; for each stepout s-t, three edges are removed,
% 2 edges are added, and 1 node removed

% the ones with no intersection points are likely to be twisted stepouts or
% those with T-motifs; the ones with twists need to be fixed in QGIS, the
% others can be re-processed with the flattening stepout function

%%
G4 = G3; 
XY4=XY3;
G4.Nodes.Name = [];

k=1;
 Edges_to_add = []; Edges_to_remove = []; Nodes_to_remove =[];

% listing the S_T, L1, L2, R1, R2 nodes for all stepouts
for i=1:numel(poss_stepout(:,1))
 disp(i)   
 poss_stepout_s_t(i,1) = poss_stepout{i,1}(1,1);                            % s node 
 poss_stepout_s_t(i,2) = poss_stepout{i,1}(1,2);                            % t node
    % neighbours of s node
 poss_stepout_s_t(i,3:4) = setdiff(neighbors(G3,poss_stepout_s_t(i,1)),poss_stepout_s_t(i,2))';
    % neighbours of t node
 poss_stepout_s_t(i,5:6) = setdiff(neighbors(G3,poss_stepout_s_t(i,2)),poss_stepout_s_t(i,1))';   
end

% listing the 5 edges that comprise stepout motif in graph edge indices
edge_counter = [poss_stepout_s_t(:,1:2); [poss_stepout_s_t(:,1) poss_stepout_s_t(:,3)];...
    [poss_stepout_s_t(:,1) poss_stepout_s_t(:,4)]; [poss_stepout_s_t(:,2) poss_stepout_s_t(:,5)];...
    [poss_stepout_s_t(:,2) poss_stepout_s_t(:,6)]];
edge_counter(:,3) = 1;
edge_counter(:,4) = findedge(G3,edge_counter(:,1),edge_counter(:,2));
G3.Nodes.Name = string(1:1:G3.numnodes)';

edge_list =    [  [findedge(G3,poss_stepout_s_t(:,1),poss_stepout_s_t(:,2))]...  % S-T  edge index
                  [findedge(G3,poss_stepout_s_t(:,1),poss_stepout_s_t(:,3))]...  % L1-S edge index
                  [findedge(G3,poss_stepout_s_t(:,1),poss_stepout_s_t(:,4))]...  % L2-S edge index
                  [findedge(G3,poss_stepout_s_t(:,2),poss_stepout_s_t(:,5))]...  % R1-S edge index
                  [findedge(G3,poss_stepout_s_t(:,2),poss_stepout_s_t(:,6))]];   % R2-S edge index

% checking for shared edges between S-T and L1,L2, S-T and R1,R2
for i=1:numel(edge_list(:,1))
   disp(i) 
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

poss_stepout_no_shared_edges = poss_stepout(no_shared_edges,1);
poss_stepout_s_t_no_shared_edges = poss_stepout_s_t(no_shared_edges,:);

%% loop that collapses the stepouts that do not overlap
for i=1:numel(poss_stepout_no_shared_edges(:,1))
    tic
    disp(i);
    
    combs_idx = [3 5 4 6; ... % L1-R1 intersects L2-R2  combination index structure      
                 3 6 4 5];    % L1-R2 intersects L2-R1
    E1 = [XY3(poss_stepout_s_t_no_shared_edges(i,combs_idx(:,1)),:) XY3(poss_stepout_s_t_no_shared_edges(i,combs_idx(:,2)),:)];
    E2 = [XY3(poss_stepout_s_t_no_shared_edges(i,combs_idx(:,3)),:) XY3(poss_stepout_s_t_no_shared_edges(i,combs_idx(:,4)),:)];
    intersect_node = intersectEdges(E1,E2);
    intersect_node = intersect_node(~isnan(intersect_node))';
     
    if ~isempty(intersect_node)
        
        % checking edge counter for edge removal history
        
        [r11,~] = findedge(G3,poss_stepout_s_t_no_shared_edges(i,1),poss_stepout_s_t_no_shared_edges(i,2));
         r12 = find(edge_counter(:,4)==r11);
        
        [r21,~] = findedge(G3,poss_stepout_s_t_no_shared_edges(i,1),poss_stepout_s_t_no_shared_edges(i,3));
         r22 = find(edge_counter(:,4)==r21);
        
        [r31,~] = findedge(G3,poss_stepout_s_t_no_shared_edges(i,1),poss_stepout_s_t_no_shared_edges(i,4));
         r32 = find(edge_counter(:,4)==r31);
        
         % appending edges to Edges_to_remove when counter is 1
        if edge_counter(r12,3)~=0
          edge_counter(r12,3)=0; 
          Edges_to_remove = [Edges_to_remove; r11];
        end
        if edge_counter(r22,3)~=0
         edge_counter(r22,3)=0;
         Edges_to_remove = [Edges_to_remove; r21];
        end
        if edge_counter(r32,2)~=0
         edge_counter(r32,3)=0;
         Edges_to_remove = [Edges_to_remove; r31];   
        end    

        Nodes_to_remove = [Nodes_to_remove; poss_stepout_s_t_no_shared_edges(i,1)];
        Edges_to_add = [Edges_to_add; [poss_stepout_s_t_no_shared_edges(i,3) poss_stepout_s_t_no_shared_edges(i,2)]; ...
                                      [poss_stepout_s_t_no_shared_edges(i,4) poss_stepout_s_t_no_shared_edges(i,2)] ];                                           
    else
       no_intersection(k,1)=i;
       k=k+1;
    end    

    toc
end

%% loop that collapses the stepouts that share S-T edges

PSSTSE = poss_stepout(all_shared_edges_idx,1);
PSSTSE_ST = find_LS_SR_stepout_edges(PSSTSE,G3);
% poss_stepout_s_t_shared_edges = poss_stepout_s_t(all_shared_edges_idx,:);

% identifying rows in which stepouts are linked to each other with one
% intersection node shared by S-T edges; these are stored in rows_linked_stepouts 
if ~isempty(PSSTSE_ST)
    rows_linked_stepouts=[];
    k=1;
    rows_shared_edges =[];
    for i=1:numel(PSSTSE_ST(:,1))
     disp(i)   
     if ~isempty(find(PSSTSE_ST(:,1)==PSSTSE_ST(i,2))) || ~isempty(find(PSSTSE_ST(:,2)==PSSTSE_ST(i,1))) ==1
        if ~isempty(find(PSSTSE_ST(:,1)==PSSTSE_ST(i,2)))
         rows_linked_stepouts{k,1} = find(PSSTSE_ST(:,1)==PSSTSE_ST(i,2)); 
         rows_linked_stepouts{k,2} = i;
         k=k+1;
        end 
    %  elseif numel(find(PSSTSE_ST(:,1)==PSSTSE_ST(i,1)))>1
    %      A{1,1}=find(PSSTSE_ST(:,1)==PSSTSE_ST(i,1))  ;
    %      if numel(cell2mat(A))==2
    %          rows_linked_stepouts{k,1} = A{1,1}(1,1); 
    %          rows_linked_stepouts{k,2} = A{1,1}(2,1);
    %          k=k+1;
    %      end
     clearvars A    
     end
    end
end
% it is possible that three stepouts share s-t nodes; in such a case we 
% need to find the rows where these cases occur and if an intersec

% following loop is a diagnostic to check how many stepouts coincide and
% whether 
% q=1;
% for j=1:numel(rows_linked_stepouts(:,1))
%     if numel(rows_linked_stepouts{j,1})>1
%        r_more_than_2_stpouts(q,1)=j;
%        more_than_2_stpouts{q,1} = rows_linked_stepouts{j,1};
%        more_than_2_stpouts{q,2} = rows_linked_stepouts{j,2};
%        more_than_2_stpouts{q,3} = intersect([PSSTSE_ST( more_than_2_stpouts{q,1}(1,1),1:2) ...
%                                              PSSTSE_ST( more_than_2_stpouts{q,1}(2,1),1:2)], ...
%                                              PSSTSE_ST( more_than_2_stpouts{q,2},1:2));
%        q=q+1;
%     end    
% end

%% loop that collapses stepouts sharing S-T edges

Edges_to_add2 =[];
Edges_to_remove2 =[];
Nodes_to_remove2=[];
if exist('rows_linked_stepouts','var')
    if ~isempty(rows_linked_stepouts)

        if numel(PSSTSE_ST(:,1))== numel(rows_linked_stepouts(:,1))*2
           disp('All stepouts are pair connected ')         
        end    

        for i=1:numel(rows_linked_stepouts(:,1)) 
            disp(i)

            A=[PSSTSE_ST(rows_linked_stepouts{i,1},:);PSSTSE_ST(rows_linked_stepouts{i,2},:)];
            stpout_nodes = unique([A(:,1) A(:,2)]); 
            all_nodes = unique(A);
            collapse_centre = intersect(A(1,1:2),A(2:end,1:2));

            B = [];
            for j = 1:numel(A(:,1))
             [a,b]=ndgrid(A(j,1),A(j,3));             
                          B = [B;[a(:),b(:)]];  

             [a,b]=ndgrid(A(j,1),A(j,4));             
                          B = [B;[a(:),b(:)]];  

             [a,b]=ndgrid(A(j,2),A(j,5));             
                          B = [B;[a(:),b(:)]]; 

             [a,b]=ndgrid(A(j,2),A(j,6));             
                          B = [B;[a(:),b(:)]]; 
            end
            B = sort(B,2);
            B = unique(B,'rows');

            Edges_to_add2 = [Edges_to_add2; [setdiff(all_nodes,stpout_nodes) ...
                repmat(collapse_centre, [numel(setdiff(all_nodes,stpout_nodes)') 1]) ]];

            Nodes_to_remove2 = [Nodes_to_remove2; setdiff( unique([PSSTSE_ST(rows_linked_stepouts{i,1},1:2);...
                PSSTSE_ST(rows_linked_stepouts{i,2},1:2)]),collapse_centre)];

            Edges_to_remove2 = [Edges_to_remove2; B];
            clearvars A

            %-----------------------------------------------------------------% 
            % following code was previous version which only worked for two 
            % stepouts sharing edges; this is modified for any number of
            % shares

    %         rows_linked_stepouts = unique(cell2mat(rows_linked_stepouts),'rows');

    %         PSSTSE_pairs = PSSTSE_ST(unique(rows_linked_stepouts),:);

    %         collapse_centre = intersect(PSSTSE_ST(rows_linked_stepouts{i,1},1:2),...
    %             PSSTSE_ST(rows_linked_stepouts{i,2},1:2));
    %         Nodes_to_remove = [Nodes_to_remove; setdiff( unique([PSSTSE_ST(rows_linked_stepouts{i,1},1:2);...
    %             PSSTSE_ST(rows_linked_stepouts{i,2},1:2)]),collapse_centre)];

    %         Edges_to_remove = [Edges_to_remove; PSSTSE_ST(rows_linked_stepouts{i,1},1:2); ...
    %                 [repmat(PSSTSE_ST(rows_linked_stepouts{i,1},1),[2 1]) PSSTSE_ST(rows_linked_stepouts{i,1},3:4)'];...
    %                 [repmat(PSSTSE_ST(rows_linked_stepouts{i,1},2),[2 1]) PSSTSE_ST(rows_linked_stepouts{i,1},5:6)'];  ...
    %                 PSSTSE_ST(rows_linked_stepouts{i,2},1:2); ...
    %                 [repmat(PSSTSE_ST(rows_linked_stepouts{i,2},1),[2 1]) PSSTSE_ST(rows_linked_stepouts{i,2},3:4)'];  ...
    %                 [repmat(PSSTSE_ST(rows_linked_stepouts{i,2},2),[2 1]) PSSTSE_ST(rows_linked_stepouts{i,2},5:6)'] ];

    %         stpout_nodes = unique([PSSTSE_ST(rows_linked_stepouts{i,1},1:2) PSSTSE_ST(rows_linked_stepouts{i,2},1:2)]);
    %         all_nodes = unique([PSSTSE_ST(rows_linked_stepouts{i,1},:) PSSTSE_ST(rows_linked_stepouts{i,2},:)]);
    %         Edges_to_add = [Edges_to_add; [setdiff(all_nodes,stpout_nodes)' repmat(collapse_centre, [numel(setdiff(all_nodes,stpout_nodes)') 1])]]; 

            %  ----------------use the following code snippet to check graph manipulations    
        %     figure(1)
        %     sub_G3 = subgraph(G3, unique([PSSTSE(linked_stepouts(i,1),:) PSSTSE(linked_stepouts(i,2),:)]));
        %       plot(sub_G3)
        %       PSSTSE_pairs(linked_stepouts(i,1),1:2)
        %     PSSTSE_pairs(linked_stepouts(i,2),1:2)
        %     %%
    %         Edges_to_remove = unique(Edges_to_remove,'rows');
    %         sub_G3 = rmedge(sub_G3,string(Edges_to_remove(:,1)),string(Edges_to_remove(:,2)));
    %         sub_G3 = addedge(sub_G3,string(Edges_to_add(:,1)),string(Edges_to_add(:,2)));
    %         sub_G3 = rmnode(sub_G3,string(Nodes_to_remove));
    %         figure(2)
    %         plot(sub_G3)
    %        clearvars sub_G3

        %     if check_if_cycle_graph(sub_G3)==1
        %       disp('cycle exists within stepout pairs')
        %       cycle_pair(m,1:2)=[linked_stepouts(i,1) linked_stepouts(i,2)];
        %       m=m+1;
        %     end 
        end 
    end
end
%%
% appending edges,nodes manipulations for the stepouts with shared edges

if ~isempty(Edges_to_remove2)   % if 
 A = findedge(G4, Edges_to_remove2(:,1),Edges_to_remove2(:,2));
 Edges_to_remove = [Edges_to_remove; A];
end
Edges_to_add = [Edges_to_add; Edges_to_add2];
Nodes_to_remove2 = [Nodes_to_remove; Nodes_to_remove2];

if exist('no_intersection','var')
  uncollapsable_stepout = poss_stepout(no_intersection,:);
else
  uncollapsable_stepout = 'NaN';
end    

G4.Nodes.Name = string(1:1:G4.numnodes)';

[G4,XY4] = update_after_collapse(G4,XY4,Edges_to_add,Edges_to_remove,Nodes_to_remove);


end

