function [Combined_Edges_cleaned_duplicate_removed] = resolve_artificial_fragmentation(G2,XY)
%UNTITLED2 Summary of this function goes here
    % this section is a 3 tier nested loop (not too slow). Parses through each
    % non-degree 2 node within the graph, finds possible chain connections  
    % from the non-degree 2 node through degree-2 nodes and each node with 
    % degree 2 (Loop 2) and then checks how far degree 2
    D2 = degree(G2);
    XY(:,3)=D2;     % add degree to the node list
    Node_neighbors_to_check  = find(D2~=2);

    Combined_Edges_chains1={};
    for j=1:length(Node_neighbors_to_check)
       % disp(j)
        nodeID = Node_neighbors_to_check(j,1);

        % queue of all nodes connected directly and indirectly to selected node
        % bfs_queue = bfsearch(G2,nodeID); 
        % immediate neighbours of selected node within graph and degrees
        neigh_nodes = [neighbors(G2,nodeID), degree(G2,neighbors(G2,nodeID))];

        % calculating degree of all nodes within the BFS queue
        % v(:,2)=degree(G2,v(:,1));
        % all nodes which have degree of 2 since we need to go deeper 
        neigh_nodes_deg_2 = neigh_nodes(find(neigh_nodes(:,2)==2),:);
        % all nodes which have degree not equal to 2 since they are 
        % possible_paths = neigh_nodes(find(neigh_nodes(:,2)~=2),:);

        origin_node = nodeID;
        [m,n]=size(neigh_nodes_deg_2);
        for i=1:m
          current_node = neigh_nodes_deg_2(i,1);  
          preceding_node = origin_node;
          edge_path = [origin_node;current_node];
          counter = 1;
          % note: if the graph contains a closed sub-graph, the while loop will
          % not end. Stop, identify the loop and correct it within the
          % shapefile
          while isempty(nieghbors_degree_2(G2,current_node,preceding_node))~=1 
            succeeding_node = nieghbors_degree_2(G2,current_node,preceding_node);
            preceding_node = current_node;
            current_node = succeeding_node;
            edge_path = [edge_path; current_node]; 
            counter = counter + 1;
          end
          last_node = nieghbors_degree_not_2(G2,current_node,preceding_node);
          edge_path=[edge_path;last_node];
          Combined_Edges_chains1{j,i}=edge_path;
        end   
        clearvars edge_path neigh_nodes neigh_nodes_deg_2
    end

    clearvars origin_node preceding_node current_node last_node counter neigh_nodes nodeID succeeding_node neigh_nodes_deg_2

    %
    Combined_Edges_chains2={};
    for j=1:length(Node_neighbors_to_check)
      %disp(j)
      nodeID = Node_neighbors_to_check(j,1);
      neigh_nodes = [neighbors(G2,nodeID), degree(G2,neighbors(G2,nodeID))];
      neigh_nodes_deg_not_2 = neigh_nodes(find(neigh_nodes(:,2)~=2),:);
      origin_node = nodeID;
      [m,n]=size(neigh_nodes_deg_not_2);
          for i=1:m
            current_node = neigh_nodes_deg_not_2(i,1);  
            preceding_node = origin_node;
            edge_path = [origin_node;current_node];
    %       counter = 1;
          % note: if the graph contains a closed sub-graph, the while loop will
          % not end. Stop, identify the loop and correct it within the
          % shapefile
    %       while isempty(nieghbors_degree_not_2(G2,current_node,preceding_node))~=1 
    %         succeeding_node = nieghbors_degree_not_2(G2,current_node,preceding_node);
    %         preceding_node = current_node;
    %         current_node = succeeding_node;
    %         edge_path = [edge_path; current_node]; 
    %         counter = counter + 1;
    %       end
    %       last_node = nieghbors_degree_not_2(G2,current_node,preceding_node);
    %       edge_path=[edge_path;last_node];
            Combined_Edges_chains2{j,i}=edge_path;
    %         clearvars edge_path
        end 
    end
    clearvars current_node preceding_node edge_path

    % reshaping the Combined_Edges_chains1 cell array to remove empty entries. 
    [m,n] = size(Combined_Edges_chains1);
    Combined_Edges_chains1_reshaped = reshape(Combined_Edges_chains1',m*n,1);
    empty_cells=find(cellfun(@isempty,Combined_Edges_chains1_reshaped));
    Combined_Edges_chains1_reshaped(empty_cells,:)=[];

    clearvars empty_cells m n Combined_Edges

    % adding the degree of the nodes to Combined_Edges_reshaped

    for i=1:length(Combined_Edges_chains1_reshaped)
     Combined_Edges_chains1_reshaped{i,2} = degree(G2,Combined_Edges_chains1_reshaped{i,1});
    end    
    clearvars i


    % reshaping the Combined_Edges_chains1 cell array to remove empty entries. 
    [m,n] = size(Combined_Edges_chains2);

    Combined_Edges_chains2_reshaped = reshape(Combined_Edges_chains2',m*n,1);
    empty_cells=find(cellfun(@isempty,Combined_Edges_chains2_reshaped));

    Combined_Edges_chains2_reshaped(empty_cells,:)=[];

    clearvars empty_cells m n Combined_Edges

    % adding the degree of the nodes to Combined_Edges_reshaped

    for i=1:length(Combined_Edges_chains2_reshaped)
     Combined_Edges_chains2_reshaped{i,2} = degree(G2,Combined_Edges_chains2_reshaped{i,1});
    end    
    clearvars i

    %
    % j=1;
    % for i=1:length(Combined_Edges_chains2_reshaped)
    %    chain2_degree =  Combined_Edges_chains2_reshaped{i,2}';
    %    if isequal(chain2_degree,[1 1])
    %       I_node_cells(j,1)=i;
    %       j=j+1;
    %    end    
    % end    
    % 
    % % not sure why we need to remove 1 - 1 edges, 
    % if exist('I_node_cells','var')
    %  Combined_Edges_chains2_reshaped(I_node_cells,:)=[];
    % end
    %
    % Combined_Edges_reshaped = [Combined_Edges_chains1_reshaped;Combined_Edges_chains2_reshaped];

    % Trying to get some statistics on edge types. Cell arrays with numeric
    % data cannot be searched. Conversion to character array is required for
    % search, or compare operations

    for i =1:length(Combined_Edges_chains1_reshaped)
       % disp(i)
        Combined_Edges_chains1_reshaped{i,3} = num2str(Combined_Edges_chains1_reshaped{i,2}');
        Combined_Edges_chains1_reshaped{i,4} = num2str(Combined_Edges_chains1_reshaped{i,1}');
    end
    clearvars i
    topology_combs_chains1 = unique(Combined_Edges_chains1_reshaped(:,3));
    % Combined_Edges_chains1 = unique(Combined_Edges_chains1_reshaped(:,4));


    % because of conversion to character arrays, unique on the character cell
    % array entries does not account for reversed entries '1 3' ~= '3 1' 
    % clearvars repeated_topo_type
    j=1;
    k=1;
    for i=1:length(topology_combs_chains1)
      topo_type= topology_combs_chains1{i,1};
      flip_topo_type = fliplr(topo_type);
      if ~isequal(flip_topo_type,topo_type)
         repeated_topo_type(j,1)=find(strcmp(topology_combs_chains1(:,1),topo_type));
         repeated_topo_type(j,2)=find(strcmp(topology_combs_chains1(:,1),flip_topo_type)); 
         j=j+1;
      else
         rows_to_retain(k,1)=i; % only rows where topology is mirrored
         k=k+1;
      end    
      clearvars topo_type flip_topo_type
    end
    clearvars i j k

    % finding row indices of the particular topology combination and the total 
    % number of incidences pertaining to each type
    for i=1:length(topology_combs_chains1)
        %disp(i)
        topology_combs_chains1{i,2}= find(strcmp(Combined_Edges_chains1_reshaped(:,3),topology_combs_chains1(i,1)));
        topology_combs_chains1{i,3}= num2str(length(topology_combs_chains1{i,2}));
    end 
    clearvars i


    % since the repetitions are mirrored, the first half may be retained or 
    % deleted. The first half id removed from topology_combs_chains1
    % rows_to_retain = [rows_to_retain; repeated_topo_type(1:floor(length(repeated_topo_type)/2),1)];
    % topology_combs_chains1_retained = topology_combs_chains1(rows_to_retain,:);

    % the order of repetitions is not structured, hence the repeated topological
    % edges are reduced
    repeated_topo_type_shrinked = repeated_topo_type;
    j=1;
    for i=1:length(repeated_topo_type)
     %disp(i)   
     b = repeated_topo_type_shrinked(i,2);
     c = find(repeated_topo_type_shrinked(:,1)==b);
     if ~isempty(repeated_topo_type_shrinked(c,1))
      repeated_topo_type_shrinked(i,:)=0;
      rows_to_shrink(j,1)=i; 
      j=j+1;
     end 
     clearvars a b c
    end

    repeated_topo_type_shrinked(repeated_topo_type_shrinked(:,1)==0,:)=[];
    rows_to_retain = [rows_to_retain; repeated_topo_type_shrinked(:,1)];
    topology_combs_chains1_retained = topology_combs_chains1(rows_to_retain(:,1),:);

    % the unique topological chains are stored in Combined_Edges_chains2_cleaned
    % the size of this cell array is smaller than Combined_Edges_chains2_reshaped
    % because of removal of duplicates
    Combined_Edges_chains1_cleaned={};
    for i=1:length(topology_combs_chains1_retained)
      %  disp(i)
        rows_to_add = topology_combs_chains1_retained{i,2};
        Combined_Edges_chains1_cleaned = [Combined_Edges_chains1_cleaned;Combined_Edges_chains1_reshaped(rows_to_add,:)];
        clearvars rows_to_add 
    end   

    %

    % for chains 2, (no embedded 2 degree nodes)
    for i =1:length(Combined_Edges_chains2_reshaped)
       % disp(i)
        Combined_Edges_chains2_reshaped{i,3} = num2str(Combined_Edges_chains2_reshaped{i,2}');
        Combined_Edges_chains2_reshaped{i,4} = num2str(Combined_Edges_chains2_reshaped{i,1}');
    end
    topology_combs_chains2 = unique(Combined_Edges_chains2_reshaped(:,3));
    Combined_Edges_chains2 = unique(Combined_Edges_chains2_reshaped(:,4));

    % because of conversion to character arrays, unique on the character cell
    % array entries does not account for reversed entries '1 3' ~= '3 1' 
    clearvars repeated_topo_type rows_to_retain
    j=1;
    k=1;
    for i=1:length(topology_combs_chains2)
      topo_type= topology_combs_chains2{i,1};
      flip_topo_type = fliplr(topo_type);
      if ~isequal(flip_topo_type,topo_type)
         repeated_topo_type(j,1)=find(strcmp(topology_combs_chains2(:,1),topo_type));
         repeated_topo_type(j,2)=find(strcmp(topology_combs_chains2(:,1),flip_topo_type)); 
         j=j+1;
      else
         rows_to_retain(k,1)=i;
         k=k+1;
      end    
      clearvars topo_type flip_topo_type
    end

    % finding row indices of the particular topology combination and the total 
    % number of incidences pertaining to each type
    for i=1:length(topology_combs_chains2)
        %disp(i)
        topology_combs_chains2{i,2}= find(strcmp(Combined_Edges_chains2_reshaped(:,3),topology_combs_chains2(i,1)));
        topology_combs_chains2{i,3}= num2str(length(topology_combs_chains2{i,2}));
    end 
    clearvars i

    % since the repetitions are mirrored, the first half may be retained or 
    % deleted. The first half id removed from topology_combs_chains1
    % rows_to_retain = [rows_to_retain; repeated_topo_type(1:floor(length(repeated_topo_type)/2),1)];
    % topology_combs_chains2_retained = topology_combs_chains2(rows_to_retain,:);


    % the order of repetitions is not structured, hence the repeated topological
    % edges are reduced
    clearvars rows_to_shrink
    repeated_topo_type_shrinked = repeated_topo_type;
    j=1;
    for i=1:length(repeated_topo_type)
     %disp(i)   
     b = repeated_topo_type_shrinked(i,2);
     c = find(repeated_topo_type_shrinked(:,1)==b);
     if ~isempty(repeated_topo_type_shrinked(c,1))
      repeated_topo_type_shrinked(i,:)=0;
      rows_to_shrink(j,1)=i; 
      j=j+1;
     end 
     clearvars a b c
    end

    repeated_topo_type_shrinked(repeated_topo_type_shrinked(:,1)==0,:)=[];
    rows_to_retain = [rows_to_retain; repeated_topo_type_shrinked(:,1)];
    topology_combs_chains2_retained = topology_combs_chains2(rows_to_retain(:,1),:);

    % the unique topological chains are stored in Combined_Edges_chains2_cleaned
    % the size of this cell array is smaller than Combined_Edges_chains2_reshaped
    % because of removal of duplicates
    Combined_Edges_chains2_cleaned={};
    for i=1:length(topology_combs_chains2_retained(:,1))
      %  disp(i)
        rows_to_add = topology_combs_chains2_retained{i,2};
        Combined_Edges_chains2_cleaned = [Combined_Edges_chains2_cleaned;Combined_Edges_chains2_reshaped(rows_to_add,:)];
        clearvars rows_to_add 
    end    

    %
    Combined_Edges_cleaned = [Combined_Edges_chains1_cleaned; Combined_Edges_chains2_cleaned];   

    % Combined_Edges_reshaped = [Combined_Edges_chains1_reshaped; Combined_Edges_chains2_reshaped];



    % creating the shapefile structure 
%     for i=1:length(Combined_Edges_cleaned)
%      disp(i)   
%      Chain =  unique(Combined_Edges_cleaned{i,1},'stable');   
%      [Linked_Edges(i).Geometry] = 'Line';
%      [Linked_Edges(i).X] = XY_georef(Chain,1);
%      [Linked_Edges(i).Y] = XY_georef(Chain,2);
%      clearvars Chain
%     end


    % writing the shapefile
%     outfolder='E:\PhD\Journal_Papers\Bristol_Channel_RWTH_Aachen\Lilstock_Benches\20_2_tiles\20_2_Comb_Ridges\Bin_Ridges\Shapefiles_Tiled\Bench_2\';
%     shapewrite(Linked_Edges,strcat(outfolder,InFileListShort{bigloop_i}));
    
    
%      clearvars A1 A2 B Combined_Edges_chains1 Combined_Edges_chains2 Combined_Edges_chains1_cleaned Combined_Edges_chains2_cleaned ...
%         Combined_Edges_chains1_reshaped Combined_Edges_chains2_reshaped Combined_Edges_cleaned D2 Edges_G2 EndNodes G1 G2 i j k Linked_Edges ...
%         neigh_nodes neigh_nodes_deg_not_2 Node_neighbors_to_check nodeID origin_node repeated_topo_type repeated_topo_type_shrinked ...
%         rows_to_retain rows_to_shrink topology_combs_chains1 topology_combs_chains2 topology_combs_chains1_retained topology_combs_chains2_retained ...
%         xmin XY XY_georef ymin
%     toc

    % test code block ...the order of repetitions is not structured, hence the repeated topological
    % edges are reduced
    % repeated_topo_type_shrinked = repeated_topo_type;
    % j=1;
    % for i=1:length(repeated_topo_type)
    %  disp(i)   
    %  b = repeated_topo_type_shrinked(i,2);
    %  c = find(repeated_topo_type_shrinked(:,1)==b);
    %  if ~isempty(repeated_topo_type_shrinked(c,1))
    %   repeated_topo_type_shrinked(i,:)=0;
    %   rows_to_shrink(j,1)=i; 
    %   j=j+1;
    %  end 
    %  clearvars a b c
    % end
    % 
    % repeated_topo_type_shrinked(repeated_topo_type_shrinked(:,1)==0,:)=[];

 %
%  n=8;
%  count=0;
%  for i=1:n
%      for j=1:n
%         for k=1:n
%             if i<j && j<k
%                 count=count+1;
%                 disp(i)
%                 disp(j)
%                 disp(k)
%                 disp(count)
%             end
%         end
%      end
%  end     
%             



% there seems to be repetitions within the cleaned edges, this needs to
% be checked and duplicates removed. Not clear is this is an artefact of
% resaving the clipped vector file from QGIS and the reloading it

% creating an N1N2 matrix which contains row indices of each topology type
% idea is to separately go through all topology types, convert to numeric
% array, flip and then apply unique rows. Cannot directly perform this on
% the cleaned edges as it is a cell array and because sizes change
N1N2={};
N1N2{1,1}=1;
N1N2{1,2}=length(topology_combs_chains1_retained{1,2});
for i=2:length(topology_combs_chains1_retained(:,1))
    disp(i)
    no_of_rows_topol = length(topology_combs_chains1_retained{i,2});
    N1N2{i,2}=N1N2{i-1,2}+no_of_rows_topol;
    N1N2{i,1}=N1N2{i-1,2}+1;
    clearvars no_of_rows_topol
end    
clearvars i

% parsing through chains of each topology type, left-right flip of chain in 
% ascending order, and performing sort 
rows_of_chains1 =[];
for i=1:length(N1N2)
  disp(i)  
  k=1;  
  for j=N1N2{i,1}:N1N2{i,2}
   topology_type_chains(k,:) = Combined_Edges_chains1_cleaned{j,1}';
   k=k+1;
  end
  [m,n] = size(topology_type_chains);
  for j=1:m
    disp(j)  
    if topology_type_chains(j,1)>topology_type_chains(j,end)
        topology_type_chains(j,:)=fliplr(topology_type_chains(j,:)); 
    end    
  end 
  [~,idx] = unique(topology_type_chains,'rows');
  temp_rows = (N1N2{i,1}:N1N2{i,2})';
  rows_of_chains1 = [rows_of_chains1; temp_rows(idx,:)];
  clearvars topology_type_chains idx temp_rows
end    

Combined_Edges_chains1_cleaned_duplicate_removed = Combined_Edges_chains1_cleaned(rows_of_chains1,:);

% 
clearvars no_of_rows_topol

% creating an N1N2 matrix which contains row indices of each topology type
% idea is to separately go through all topology types, convert to numeric
% array, flip and then apply unique rows. Cannot directly perform this on
% the cleaned edges as it is a cell array and because sizes change
N1N2={};
N1N2{1,1}=1;
N1N2{1,2}=length(topology_combs_chains2_retained{1,2});    
for i=2:length(topology_combs_chains2_retained(:,1))
    disp(i)
    no_of_rows_topol = length(topology_combs_chains2_retained{i,2});
    N1N2{i,2}=N1N2{i-1,2}+no_of_rows_topol;
    N1N2{i,1}=N1N2{i-1,2}+1;
    clearvars no_of_rows_topol
end    
clearvars i

% parsing through chains of each topology type, left-right flip of chain in 
% ascending order, and performing sort 
clearvars topology_type_chains idx temp_rows
rows_of_chains2 =[];
for i=1:length(N1N2)
  disp(i)  
  k=1;  
  for j=N1N2{i,1}:N1N2{i,2}
   topology_type_chains(k,:) = Combined_Edges_chains2_cleaned{j,1}';
   k=k+1;
  end
  [m,n] = size(topology_type_chains);
  for j=1:m
    disp(j)  
    if topology_type_chains(j,1)>topology_type_chains(j,end)
        topology_type_chains(j,:)=fliplr(topology_type_chains(j,:)); 
    end    
  end 
  [~,idx] = unique(topology_type_chains,'rows');
  temp_rows = (N1N2{i,1}:N1N2{i,2})';
  rows_of_chains2 = [rows_of_chains2; temp_rows(idx,:)];
  clearvars topology_type_chains idx temp_rows
end    

Combined_Edges_chains2_cleaned_duplicate_removed = Combined_Edges_chains2_cleaned(rows_of_chains2,:);

%
Combined_Edges_cleaned_duplicate_removed = [Combined_Edges_chains1_cleaned_duplicate_removed; Combined_Edges_chains2_cleaned_duplicate_removed];   


end

