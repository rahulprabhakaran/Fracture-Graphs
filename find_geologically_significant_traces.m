function [unique_walks] = find_geologically_significant_traces(G3,XY3,strike_threshold)
% this function takes as input a graph and its spatial positioning matrix
% along with an input strike threshold, and converts it into a cell array
% of node sequences that are geologically consistent

Edges_G3=table2array(G3.Edges);   % listing all graph edges
if iscell(Edges_G3)==1 
   G3.Nodes.Name = [];
   Edges_G3=table2array(G3.Edges);
end    
Edges_G3(:,3)=compute_strike(Edges_G3, XY3); % normalized strike of edges
Edges_G3_copy=Edges_G3;

edge_counter = table2array(G3.Edges); 
edge_counter(:,3) = 1; % the 3rd col is a counter set initially to 1
G3.Nodes.Name = string(1:1:G3.numnodes)';

components=[];
k=1;
m=1;

for i=1:numel(Edges_G3_copy(:,1))
   tic 
%   i
    disp(i)
   chain = sort(Edges_G3_copy(i,1:2),2);
   node_s = Edges_G3_copy(i,1);  
   walks{i,1} = find_possible_edges3(G3,node_s,chain,XY3,strike_threshold);  
   node_t = Edges_G3_copy(i,2);
   chain2 = fliplr(chain);
   walks{i,2} = find_possible_edges3(G3,node_t,chain2,XY3,strike_threshold); 
   walks{i,3} = [walks{i,1} fliplr(walks{i,2})];
   [~,u,~]=unique(walks{i,3}, 'first') ;
   walks{i,4} = walks{i,3}(sort(u));
   
   %---------------------------------------------------------------------%
   % checking edge-counter, extracting only 'ON' edges and saving as a
   % subgraph
   Edges = chain_to_edges(walks{i,4});
   edge_idx1 =  setdiff(findedge(G3,Edges(:,1),Edges(:,2)),0); % prevent zero idx
   edge_idx2 = find(edge_counter(edge_idx1,3)==1);  % finding edges in "on" condition
   Edges = Edges(edge_idx2,:);                      % extracting only "on" edges
   edge_idx1 = edge_idx1(edge_idx2,1);              % idx of "on" edges
   path_graph = create_subgraph(unique(Edges,'stable')',G3,XY3);
   % checks: if subgraph is not a line graph and is not empty
   if path_graph.numnodes~=0 && path_graph.numedges~=0 && check_if_line_graph(path_graph)==0
       bif{m,1} = path_graph;
       bif{m,2} = i;
       m=m+1;
   elseif numel(unique(conncomp(path_graph)))>1 % check if there are more than one components   
       components{k,1}=path_graph;
       components{k,2}=i;
       k=k+1;
   else
       walks{i,5} = edges_to_chain(path_graph);    
   end    
   edge_counter(edge_idx1,3) = 0;   
   walks{i,6} = num2str(walks{i,5}); 
   %---------------------------------------------------------------------%
   
   % walks{i,6} = unique(cellfun(@str2num,table2array(path_graph.Edges)),'stable')';   
   % this method of assigning walk node order was wrong a few times
%    if numel(walks{i,6})==2
%        walks{i,6} = walks{i,6}';
%    end  
%    walks{i,7} = num2str(walks{i,5});
     
%    if numel(walks{i,5})==numel(unique(walks{i,5}))
%        walks{i,9} = 0;
%    else
%        walks{i,9} = 1;
%    end    
%    if ~isequal(walks{i,5},walks{i,6}) && ~isempty(walks{i,5})
%        walks{i,10} = 1;
%    else
%        walks{i,10} = 0;
%    end    
   clearvars node_s node_t chain Edges edge_idx1 edge_idx2 edges path_graph
   toc
end   

%% some bifurcations have components within them; identifying these
if exist('bif','var')
    rows_comp = [];
    for i=1:numel(bif(:,1))
      comp = unique(conncomp(bif{i,1}));  
      if numel(comp)>1
         rows_comp = [rows_comp; i]; 
      end   
    end
    rows_bif=setdiff((1:1:numel(bif(:,1)))',rows_comp);
    components = bif(rows_comp,:);
    bif = bif(rows_bif,:);
    if exist('bif_idx','var')
     bif_idx = bif_idx(rows_bif,:);
    end
    clearvars comp rows_comp rows_bif
end


%% converting components to walks
if ~isempty(components)
    paths_to_add ={};
    for i=1:numel(components(:,1))
       g1 =  components{i,1};
       n_comp = numel(unique(conncomp(g1)));
       for j=1:n_comp
           full_path = str2double(table2array(g1.Nodes))';
           comp{j,1}=full_path(1,find(conncomp(g1)==j));
           g2 = subgraph(G3,comp{j,1});
           if check_if_line_graph(g2)==1
             paths_to_add = [paths_to_add; edges_to_chain(g2)];
           elseif check_if_fork_graph(g2)==1
                  [~,bif_node] =  check_if_fork_graph(g2);
                  paths_to_add = [paths_to_add ; find_subgraph_bifurcations1(g2,bif_node)];
           else
                  paths_to_add = [paths_to_add ; find_subgraph_bifurcations3(g2,XY3)];
           end 
           clearvars g2 bif_node 
       end  
       clearvars g1 n_comp full_path
    end
end

%% sub_graph rows with bifurcations need to be split; these are very few
bifurcations = [];



if exist('bif','var')
    
%[bif_idx,ia] = setdiff(cell2mat(bif(:,2)),cell2mat(components(:,2)));
[bif_idx] = cell2mat(bif(:,2));         % when there are no components

    for i=1:numel(bif_idx)
      disp(i)  
    %   sub_G3  = bif{ia(i,1),1};      
      sub_G3 =  bif{i,1};                   % when there are no components

      if check_if_fork_graph(sub_G3)==1
          [~,bif_node] =  check_if_fork_graph(sub_G3);
          bifurcations = [bifurcations ; find_subgraph_bifurcations1(sub_G3,bif_node)];
      else
          bifurcations = [bifurcations ; find_subgraph_bifurcations3(sub_G3,XY3)];   
      end    
      clearvars sub_G3  bif_node 
    end    
end

%%

if exist('bif','var')
    walks(bif_idx,:)={[]};                   % removing paths with bifurcations
end
if ~isempty(components)
    walks(cell2mat(components(:,2)),:)={[]}; % removing paths with multiple components
end
%%

unique_walks = walks(:,5);
unique_walks = unique_walks(~cellfun('isempty',unique_walks));

%%
if ~isempty(bifurcations)
    unique_walks = [unique_walks; bifurcations];
end
if exist('paths_to_add','var')
 %if ~isempty(bifurcations) && ~isempty(paths_to_add) 
    unique_walks = [unique_walks; paths_to_add; bifurcations];
 %end
end
unique_walks(:,2) = cellfun(@num2str,unique_walks(:,1),'UniformOutput',false);

%% identifying duplicates
for i=1:length(unique_walks)
 disp(i) 
 try
 rows_to_remove{i,1} = convertCharsToStrings(num2str(sort([find(ismember(unique_walks(:,2),num2str(str2num(unique_walks{i,2})))==1); ...
     find( ismember(unique_walks(:,2),num2str(fliplr(str2num(unique_walks{i,2}))))==1) ]) )) ; 
 end
end

%% retaining only unique walks
if exist('rows_to_remove','var')
    rows_to_remove = rows_to_remove(~cellfun('isempty',rows_to_remove));
    rows_to_remove = cellstr(rows_to_remove);
    [~, ia1, ~] = unique(rows_to_remove);
    ia1=sort(ia1);
    unique_walks = unique_walks(ia1,:);
    clearvars ia1 rows_to_remove
end

%% checking for bifurcations within the walks
rows_bif = []; k=1;
for i=1:numel(unique_walks(:,1))
    disp(i)
    out = check_bifurcations(unique_walks{i,1},G3);
    if out~=0       
       rows_bif(k,1)=i;
       k=k+1;
    end        
end

%%
bif_walks=[];
for i=1:numel(rows_bif)
    sub_G3=subgraph(G3,unique_walks{rows_bif(i,1),1});
    if check_if_fork_graph(sub_G3)==1
      [~,bif_node] =  check_if_fork_graph(sub_G3);
      bif_walks = [bif_walks ; find_subgraph_bifurcations1(sub_G3,bif_node)];
    else
      bif_walks = [bif_walks ; find_subgraph_bifurcations3(sub_G3,XY3)];
    end  
    clearvars sub_G3 bif_node
end    

unique_walks(rows_bif,:) = {[]};
unique_walks = [unique_walks(:,1); bif_walks];
unique_walks = unique_walks(~cellfun('isempty',unique_walks));

for i=1:numel(unique_walks)
   [m,n] = size(unique_walks{i,1});
   if m>n
      unique_walks{i,1}= unique_walks{i,1}';
   end    
end    

end

