%%------------------------------------------------------------------------%
% this script is the code supplement to the manuscript, "Prabhakaran et al 
% (2021), Large-scale fracture network patterns: Insights from automated 
% mapping in the Lilstock (Bristol Channel) limestone outcrop" submitted to
% Journal of Structural Geology and contains the implementation of several 
% graph-based fracture network functions. 

% Copyright Rahul Prabhakaran, 2021

% Please report any bugs, errors, corrections to r.prabhakaran@tudelft.nl

% Requires -> MATLAB Mapping Toolbox. 
% Uses -> Geom2D Toolbox

% ----------------------Contents-----------------------------------------%
% 1. Loading a shapefile and converting it to a graph
% 2. Fixing three types of common topological discontinuities
%  2.1 Type-1 discontinuity fix
%  2.2 Type-2 discontinuity fix
%  2.3 Type-3 discontinuity fix
% 3. Resolving artificial fragmentation in shapefile attribute tables
% 4. Resolving Stepouts within fracture graphs
%  4.1 Resolving stepouts using merge operation
%  4.2 Resolving stepouts using flatten operation
% 5. Straightening fracture segments (removing degree-2 nodes)
% 6. Converting fracture graphs to geologically-significant fracture traces
% 7. Converting primals graphs to dual graph representations

% functions list:
% Shape_to_Data, build_shp_struct, chain_to_edges, collapse_stepouts, 
% compute_chain_length, create_GraphEdges, find_I_M_possible_edges,
% find_geologically_significant_traces, find_kinked_M_Node_connections,
% find_kinked_M_Nodes, find_small_I_I_I_trielements_AR, find_stepout_motif,
% fix_graph_type_5_disconnection, fix_graph_type_8_disconnection,
% flatten_stepouts, resolve_artificial_fragmentation, straighten_graph

clc
clear all
close all
format compact

%% 1. Loading a shapefile and converting it to a graph

%    2D Fracture networks are often saved as shapefiles in which the 
%    geometry is represented in the form of polylines. The shapefiles can
%    be converted to a graph representation using graph objects. This section
%    depicts how a shapefile can be loaded and converted into a graph
%    structure.

    % we load a shapefile from the open-source fracture network dataset:
    % https://data.4tu.nl/articles/dataset/Deterministic_fracture_network_models_from_the_Potiguar_basin_Brazil/12708305
    % using the shaperead function from the Mapping toolbox
    B = shaperead('E:\PhD\Journal_Papers\Manuscript_JSG\code_supplement\apodi4.shp');

    % the data is read in the form of F = [x1 y1 x2 y2] using the 
    % Shape_to_Data() function
    [F, xmin, xmax, ymin, ymax] = Shape_to_Data(B);

    % checking for and removing duplicates
    Fractures = unique(F,'rows'); 

    % the following function converts the geometry into a set of edges
    % specified in the EndNodes variable with corresponding spatial
    % positioning matrix XY
    [EndNodes,xy] = create_GraphEdges(Fractures, xmin, ymin);
    
    % 
    for i=1:numel(EndNodes(:,1))
     disp(i)
     if numel(EndNodes{i,1}{1,1})>1
         EndNodes{i,1}{1,1} = min (EndNodes{i,1}{1,1});
     end    
     EndNodes2(i,1)=EndNodes{i,1}{1,1};
     if numel(EndNodes{i,1}{1,2})>1
         EndNodes{i,1}{1,2} = min (EndNodes{i,1}{1,2});
     end 
     EndNodes2(i,2)=EndNodes{i,1}{1,2};
    end
    
    % creating and empty graph and assigning edges to it
    g=graph(); 
    g =  addedge(g,EndNodes2(:,1),EndNodes2(:,2));
    
    % once we have the graph we can compute its adjacency matrix
    a = adjacency(g);
    
    % combination of adjacency and spatial positioning is sufficient to
    % display the graph in a planar, spatial form
    figure(1)
    gplot(a,xy)
    
    % the adjacency matrix is a sparse-representation of edge connections
    % and can be visualized using:
    figure(2)
    spy(a)
    
    % the graph can be drawn without spatial positioning in different
    % layouts, the following plot depicts the number of components 
    figure(3)
    plot(g)
    
    figure(4)
    plot(g, 'layout', 'force', 'UseGravity', true)
    
    % the topology of a graph is calculated from the degree
    d = degree(g);
    
    % and the topology range visualized using a histogram of node degrees
    figure(5)
    histogram(d)
    
    % once the graph has been created, we can simply save the graph object
    % 'g' and the spatial positioning matrix 'xy' 
 
 %%   2.1 Automatic detection can lead to certain topological discontinuities
 %    that needs to be resolved so that network metrics can be applied.
 %    The first type of topological discontinuity is when the a fracture
 %    with a degree-3 node is incorrectly represented as a degree-2 node
 %    in close proximity to a degree-1 node. In this section, we implement
 %    a function to resolve this using delaunay triangulation
 
      clear all; close all; clc;    % clearing all workspace variables
      
      % to illustrate, an example graph cut from a shapefile is loaded
      % with the corresponding spatial positioning matrix
      load('E:\PhD\Journal_Papers\Manuscript_JSG\code_supplement\g_topo_1');
      load('E:\PhD\Journal_Papers\Manuscript_JSG\code_supplement\xy_topo_1')
      
      % creating a delaunay triangulation around nodes
      dt = delaunayTriangulation(xy);
      tri = dt.ConnectivityList;
      xy1 = dt.Points;
       
      figure(1)
      triplot(dt,'c')
      hold on
      a = adjacency(g);  
      gplot(a,xy,'k')
      pbaspect([1 1 1])
      
      
      d = degree(g);
      angle_1 = 15;   %20
      angle_2 = 165;  %160 
      [M_nodes_idx,~] = find(d(:,1)==2);
      M_nodes=xy(M_nodes_idx,:);
      [kinked_M_nodes] = find_kinked_M_Nodes(xy,M_nodes_idx,angle_1,angle_2,g);
      
      [I_nodes_idx,~] = find(d(:,1)==1);
      I_nodes=xy(I_nodes_idx,:);
      [possible_edges_I_M_kinked] = find_I_M_possible_edges(dt,kinked_M_nodes,I_nodes_idx,I_nodes,xy,g);

      length_threshold = 0.1;
      possible_edges_I_M_kinked = possible_edges_I_M_kinked(possible_edges_I_M_kinked(:,3)<length_threshold,:);
      
      g2 = fix_graph_type_5_disconnection(g,possible_edges_I_M_kinked);
      xy2 = xy;
      
      a2 = adjacency(g2);
      
      figure(2)
      subplot(1,3,1)
      gplot(a,xy,'k')
      xlim([5.75 6.05]); ylim([5.6 5.9])
      pbaspect([1 1 1])
      set(gca,'xtick',[]); set(gca,'ytick',[]);
      %xticks([5.75 6.05])
      %yticks([5.6 5.9])
      
      subplot(1,3,2)
      triplot(dt,'c')
      hold on
      gplot(a,xy,'k')
      hold on
      scatter(xy(possible_edges_I_M_kinked(:,1),1),xy(possible_edges_I_M_kinked(:,1),2),'r','Filled')
      hold on
      scatter(xy(possible_edges_I_M_kinked(:,2),1),xy(possible_edges_I_M_kinked(:,2),2),'r','Filled')
      xlim([5.75 6.05]); ylim([5.6 5.9])
      pbaspect([1 1 1])
      set(gca,'xtick',[]); set(gca,'ytick',[]);
      %xticks([5.75 6.05])
      %yticks([5.6 5.9])
      
      subplot(1,3,3)
      gplot(a2,xy2,'k');
      xlim([5.75 6.05]); ylim([5.6 5.9])
      pbaspect([1 1 1])
      set(gca,'xtick',[]); set(gca,'ytick',[]);
      %xticks([5.75 6.05])
      %yticks([5.6 5.9])
 
 %%   2.2 The second type of topological discontinuity is when three degree-1
 %    nodes lie in close proximity. Again the delaunay triangulation is
 %    applied and based on thresholds set for tri-element size and aspect-
 %    ratio, these discontinuities can be identified and fixed. 

      clear all; close all; clc;    % clearing all workspace variables
      
      % loading the same example
      load('E:\PhD\Journal_Papers\Manuscript_JSG\code_supplement\g_topo_1');
      load('E:\PhD\Journal_Papers\Manuscript_JSG\code_supplement\xy_topo_1')
      
      % creating a delaunay triangulation around nodes
      dt = delaunayTriangulation(xy);
      tri = dt.ConnectivityList;
      xy1 = dt.Points;
       
      figure(1)
      triplot(dt,'c')
      hold on
      a = adjacency(g);
      gplot(a,xy,'k')
      pbaspect([1 1 1])
      
      d = degree(g);
      trielement_list = dt.ConnectivityList;
      trielement_list(:,4)= d(trielement_list(:,1));
      trielement_list(:,5)= d(trielement_list(:,2));
      trielement_list(:,6)= d(trielement_list(:,3));
      
      
      % setting a threshold for tri-element areas and aspect ratios to 
      % identify I-I-I elements
      area_threshold = 0.000129; 
      aspect_ratio_threshold = 1.3;
      
      % identifying the I-I-I elements
      [small_I_I_I_trielements_AR]=find_small_I_I_I_trielements_AR...
          (trielement_list,area_threshold,aspect_ratio_threshold,xy1);
      
      % visualizing the I-I-I trielements sorted for aspect ratio and area
      figure(1)
      triplot(dt,'c')
      hold on
      gplot(a,xy1,'b')
      hold on
      triplot(tri(small_I_I_I_trielements_AR(:,4),:),xy1(:, 1), xy1(:, 2),'k','LineWidth',2);
      pbaspect([1 1 1])
      
      % resolving the topological discontinuity 
      [g2,xy2] = fix_graph_type_8_disconnection(g,small_I_I_I_trielements_AR,xy1);
      a2 = adjacency(g2);
      
      figure(2)
      gplot(a2,xy2);
      pbaspect([1 1 1])
      
      % zooming into a region in the updated graph where a topological fix 
      % has taken place and comparing with original graph
      figure(3)
      
      subplot(1,3,1)
      gplot(a,xy,'k');
      xlim([8.3 8.6]); ylim([6.2 6.5])
      pbaspect([1 1 1])
      set(gca,'xtick',[]); set(gca,'ytick',[]);
      %xticks([8.3 8.6])
      %yticks([6.2 6.5])
      
      subplot(1,3,2)
      triplot(dt,'c')
      hold on
      gplot(a,xy1,'k')
      hold on
      triplot(tri(small_I_I_I_trielements_AR(:,4),:),xy1(:, 1), xy1(:, 2),'r','LineWidth',2);
      xlim([8.3 8.6]); ylim([6.2 6.5])
      pbaspect([1 1 1])
      set(gca,'xtick',[]); set(gca,'ytick',[]);
      %xticks([8.3 8.6])
      %yticks([6.2 6.5])
      
      subplot(1,3,3)
      gplot(a2,xy2,'k');
      xlim([8.3 8.6]); ylim([6.2 6.5])
      pbaspect([1 1 1])
      set(gca,'xtick',[]); set(gca,'ytick',[]);
      %xticks([8.3 8.6])
      %yticks([6.2 6.5])

 %%   2.3 The third type of topological discontinuity is when two degree-2
 %    nodes lie in close proximity with sharp angles. Again the delaunay 
 %    triangulation is applied and based on thresholds set for tri-element 
 %    size and aspect-ratio, these discontinuities can be identified and
 %    fixed
 
 clear all; close all; clc;    % clearing all workspace variables
      
      % loading the same example
      load('E:\PhD\Journal_Papers\Manuscript_JSG\code_supplement\g_topo_1');
      load('E:\PhD\Journal_Papers\Manuscript_JSG\code_supplement\xy_topo_1')
      a= adjacency(g);
      
      % creating a delaunay triangulation around nodes
      dt = delaunayTriangulation(xy);
      tri = dt.ConnectivityList;
      xy1 = dt.Points;

          % setting angle thresholds for connection and obtaining kinked M-nodes,
    % calculating possible edges for M-M connection

    angle_1 = 25;   %20
    angle_2 = 160;  %160 

    clearvars M_nodes_idx M_nodes
    d = degree(g);
    [M_nodes_idx,~] = find(d(:,1)==2);
    M_nodes=xy(M_nodes_idx,:);

    [kinked_M_nodes] = find_kinked_M_Nodes(xy,M_nodes_idx,...
        angle_1,angle_2,g);

    [kink_connect, possible_edges_M_M] = find_kinked_M_Node_connections...
        (kinked_M_nodes,xy,dt,g);

    % applying a length threshold for M-M connections
    M_M_length_threshold=0.075; 
    possible_edges_M_M= possible_edges_M_M(possible_edges_M_M(:,3)<M_M_length_threshold,:);
    possible_edges_M_M(:,3)=[];


    % appending the possible M-M edges and creating a new graph

    g2 = fix_graph_type_5_disconnection(g,possible_edges_M_M);
    xy2 = xy; % the spatial matrix remains the same due to edge addition
      
    a2=adjacency(g2);
    
    figure(3)
    gplot(a2,xy2)
      
    figure(2)
      subplot(1,3,1)
      gplot(a,xy,'k')
      xlim([7.45 7.75]); ylim([6.55 6.85])
      pbaspect([1 1 1])
      set(gca,'xtick',[]); set(gca,'ytick',[]);
      %xticks([7.45 7.75])
      %yticks([6.55 6.85])
      
      subplot(1,3,2)
      triplot(dt,'c')
      hold on
      gplot(a,xy,'k')
      hold on
      scatter(xy(possible_edges_M_M(:,1),1),xy(possible_edges_M_M(:,1),2),'r','Filled')
      hold on
      scatter(xy(possible_edges_M_M(:,2),1),xy(possible_edges_M_M(:,2),2),'r','Filled')
      xlim([7.45 7.75]); ylim([6.55 6.85])
      pbaspect([1 1 1])
      set(gca,'xtick',[]); set(gca,'ytick',[]);
      %xticks([7.45 7.75])
      %yticks([6.55 6.85])
      
      subplot(1,3,3)
      gplot(a2,xy2,'k');
      xlim([7.45 7.75]); ylim([6.55 6.85])
      pbaspect([1 1 1])
      set(gca,'xtick',[]); set(gca,'ytick',[]);
      %xticks([7.45 7.75])
      %yticks([6.55 6.85])
 
%%   3. Resolving artificial fragmentation
    % when using automated tracing, the vectorized traces may be topologically
    % connected but fragmented. This effect is due to the image processing
    % steps that results in shapefiles of fragmented edges. This can be 
    % resolved easily by converting the shapefile to a graph and then
    % grouping edges based on degree-2 connection. We illustrate this using
    % an example:
    
    close all; clear all; clc

    % we load a shapefile from the open-source fracture network dataset:
    % https://data.4tu.nl/articles/dataset/Deterministic_fracture_network_models_from_the_Potiguar_basin_Brazil/12708305
    % using the shaperead function from the Mapping toolbox
    B = shaperead('E:\PhD\Journal_Papers\Manuscript_JSG\code_supplement\apodi4.shp');

    % the data is read in the form of F = [x1 y1 x2 y2] using the 
    % Shape_to_Data() function
    [F, xmin, xmax, ymin, ymax] = Shape_to_Data(B);

    % checking for and removing duplicates
    Fractures = unique(F,'rows'); 

    % the following function converts the geometry into a set of edges
    % specified in the EndNodes variable with corresponding spatial
    % positioning matrix XY
    [EndNodes,xy] = create_GraphEdges(Fractures, xmin, ymin);
    
    % 
    for i=1:numel(EndNodes(:,1))
     disp(i)
     if numel(EndNodes{i,1}{1,1})>1
         EndNodes{i,1}{1,1} = min (EndNodes{i,1}{1,1});
     end    
     EndNodes2(i,1)=EndNodes{i,1}{1,1};
     if numel(EndNodes{i,1}{1,2})>1
         EndNodes{i,1}{1,2} = min (EndNodes{i,1}{1,2});
     end 
     EndNodes2(i,2)=EndNodes{i,1}{1,2};
    end
    
    % creating and empty graph and assigning edges to it
    g=graph(); 
    g =  addedge(g,EndNodes2(:,1),EndNodes2(:,2));
    
    % the loaded shapefile was traced manually. Now converting the graph 
    % back into a shapefile, with all entries in the shapefile as individual 
    % edges. This is the typical situation when an automatically traced 
    % shapefile is converted into a graph.
    e = table2array(g.Edges);
    for i=1:numel(e(:,1))
       edge2shp{i,1} = e(i,:); 
    end    
    xmin = min(xy(:,1));
    ymin = min(xy(:,2));
    
    % writing this as a shapefile so that it can be viewed in QGIS
    outfilename = 'E:\PhD\Journal_Papers\Manuscript_JSG\code_supplement\Apodi4_disconnected.shp';
    fractures_shape = build_shp_struct(edge2shp,xy,xmin,ymin,1);
    shapewrite(fractures_shape,outfilename);
    
    % it can be seen that the edges are all fragmented. we would like to have 
    % edges linked so that they represent continuous segments
    
    % the following function does just that; grouping edges between 
    % non-degree 2 nodes 
    [Shape_Edges] = resolve_artificial_fragmentation(g,xy);
    
    % the output is saved as a shapefile to view the effect of the combination
    outfilename = 'E:\PhD\Journal_Papers\Manuscript_JSG\code_supplement\Apodi4_connected.shp';
    fractures_shape = build_shp_struct(Shape_Edges,xy,xmin,ymin,1);
    shapewrite(fractures_shape,outfilename);
        
%%   4.1 Resolving stepouts using merge operation
      % in this example, small stepouts based on a threshold are merged
      % or collapsed. Merging a stepout reduces two degree-3 nodes into
      % a degree-4 node

      clear all; close all; clc;    % clearing all workspace variables
      
      % loading a new example from Bristol Channel
      load('E:\PhD\Journal_Papers\Simplicity_Profiles\25_20_g\25_20_g');
      load('E:\PhD\Journal_Papers\Simplicity_Profiles\25_20_g\25_20_xy')
      
      a= adjacency(g);    % adjacency matrix
      g.Edges.Weight=[];  % removing edge weights
      
      % setting a stepout threshold for merging
      collapse_stepout_length = 0.02;   
      
      % all possible stepouts as per the threshold
      poss_stepout = find_stepout_motif(g,xy,collapse_stepout_length);
      g.Nodes.Name = string(1:1:g.numnodes)';
      
      % performing the collapse (or merge operation) and generating a new
      % graph with the merged stepouts
      [g2,xy2,~,~,~] = collapse_stepouts(poss_stepout,g,xy);
      if ismultigraph(g2)
         g2 = simplify(g2);
      end
      
      a2 = adjacency(g2);  % adjacency matrix of new graph
      poss_stp_outs = cell2mat(poss_stepout);
      
      
      % plot to depict a merged stepout
      subplot(1,2,1)
      gplot(a,xy,'k') 
      hold on
      scatter(xy(poss_stp_outs(:,1),1),xy(poss_stp_outs(:,1),2),10,'r','Filled')
      hold on
      scatter(xy(poss_stp_outs(:,2),1),xy(poss_stp_outs(:,2),2),10,'r','Filled')
      xlim([25.3 25.8]); ylim([21.7 22.2])
      pbaspect([1 1 1])
      set(gca,'xtick',[]); set(gca,'ytick',[]);
      
      subplot(1,2,2)
      gplot(a2,xy2,'k')
      xlim([25.3 25.8]); ylim([21.7 22.2])
      pbaspect([1 1 1])
      set(gca,'xtick',[]); set(gca,'ytick',[]);
      

%%   4.2 Resolving stepouts using flatten operation
      % in this example, stepouts based on length threshold and based on
      % two additional inputs of angle and orthogonal range are flattened
      % so that continuity is maintained in one direction
      
      clear all; close all; clc;    % clearing all workspace variables
      
      % loading a new example from Bristol Channel
      load('E:\PhD\Journal_Papers\Manuscript_JSG\code_supplement\g_flatten');
      load('E:\PhD\Journal_Papers\Manuscript_JSG\code_supplement\xy_flatten')
      
      a= adjacency(g);    % adjacency matrix
     
      % inputs to flatten stepouts
      collapse_stepout_length = 0.05; 
      orthogonal_range = 15;
      strike_threshold = 20;
      
      % identifying the stepout motif
      poss_stepout = find_stepout_motif(g,xy,collapse_stepout_length);
      
      % flattening the possible stepouts based on inputs into a new graph
      [g2,xy2,~,~] = flatten_stepouts(poss_stepout, orthogonal_range, strike_threshold, g, xy);
      
      a2 = adjacency(g2);  % adjacency matrix of new graph
      poss_stp_outs = cell2mat(poss_stepout);
      
      % plotting an example of a flattened stepout
      figure(1)
      subplot(1,2,1)
      gplot(a,xy,'k') 
      hold on
      scatter(xy(poss_stp_outs(:,1),1),xy(poss_stp_outs(:,1),2),15,'r','Filled')
      hold on
      scatter(xy(poss_stp_outs(:,2),1),xy(poss_stp_outs(:,2),2),15,'r','Filled')
      xlim([140.3 140.7]); ylim([35.2 35.5])
      pbaspect([1 1 1])
      set(gca,'xtick',[]); set(gca,'ytick',[]);
      
      subplot(1,2,2)
      gplot(a2,xy2,'k')
      xlim([140.3 140.7]); ylim([35.2 35.5])
      pbaspect([1 1 1])
      set(gca,'xtick',[]); set(gca,'ytick',[]);

%%  5. Straightening fracture segments
    % when converting a shapefile to a graph, there are a large number of
    % degree-2 nodes. These represent the natural sinuousity of fracture
    % traces. Within the graph representation, they increase the size of
    % the graph and on some occasion we would like to remove these degree-2
    % nodes while retaining the overall graph structure
    
     clear all; close all; clc;    % clearing all workspace variables
     
    % loading the Apodi example as a graph
    load('E:\PhD\Journal_Papers\Manuscript_JSG\code_supplement\g_apodi');
    load('E:\PhD\Journal_Papers\Manuscript_JSG\code_supplement\xy_apodi')
    
    % computing adjacency matrix
    a = adjacency(g)  ;
     
    % we need the data structure that holds the combined edges
    [Combined_Edges] = resolve_artificial_fragmentation(g,xy);
    
    % straightening the graph using this information into a new graph
    [g2,xy2,removed_V_nodes] = straighten_graph(g,xy,Combined_Edges);
    
    % adjacency matrix of the straightened graph
    a2 = adjacency(g2);
    
    % plotting to visualize the effect of straightening
    figure(1)
    subplot(1,2,1)
    gplot(a,xy);
    xlim([86 96]); ylim([95 105])
    pbaspect([1 1 1])
    set(gca,'xtick',[]); set(gca,'ytick',[]);
    %xticks([7.45 7.75])
    %yticks([6.55 6.85])
    
    subplot(1,2,2)
    gplot(a2,xy2);
    xlim([86 96]); ylim([95 105])
    pbaspect([1 1 1])
    set(gca,'xtick',[]); set(gca,'ytick',[]);
    
    % the straightening operation gets rid of the degree-2 nodes. The 
    % change can be depicted in a degree distribution comparison
    
    d = degree(g);        % degree of original graph
    d2 = degree(g2);      % degree of straightened graph   
    
    figure(2)
    subplot(1,2,1)
    histogram(d); 
    title('degree dist - original graph')
    pbaspect([1 1 1])
    
    subplot(1,2,2)
    histogram(d2);
    title('degree dist - straightened graph')
    pbaspect([1 1 1])
    
    % the straightened graph can be converted and saved as a shapefile
    e2 = table2array(g2.Edges);
    for i=1:numel(e2(:,1))
       edge2shp{i,1} = e2(i,:); 
    end    
    xmin = min(xy2(:,1));
    ymin = min(xy2(:,2));
    
    % writing this as a shapefile so that it can be viewed in QGIS
    outfilename = 'E:\PhD\Journal_Papers\Manuscript_JSG\code_supplement\Apodi4_straightened.shp';
    fractures_shape = build_shp_struct(edge2shp,xy,xmin,ymin,1);
    shapewrite(fractures_shape,outfilename);
    
      
%%   6. Converting fracture graphs to geologically-significant fracture traces     
     % in this section, we deal with conversion of graph edges into
     % geologically significant fractures. A fracture from tip-to-tip can
     % be identified from a series of graph edges provided they are
     % continuous in terms of edge angle. An input threshold strike angle
     % is set and the function identifies a set of fracture edges that are
     % within the threshold.
     
     clear all; close all; clc;    % clearing all workspace variables
     
     % an example graph is loaded
     load('E:\PhD\Journal_Papers\Manuscript_JSG\code_supplement\25_20_g.mat')
     load('E:\PhD\Journal_Papers\Manuscript_JSG\code_supplement\25_20_xy.mat')
     
     g.Edges.Weight =[];  % removing edge weights
     [unique_walks] = find_geologically_significant_traces(g,xy,20);
     
     % in order to depict the continuity of traces, we plot the fractures
     % by length
     for i=1:numel(unique_walks(:,1))
        unique_walks{i,2} =  compute_chain_length(unique_walks{i,1},xy);
     end    
     
     % creating a colorbar to depict different lengths of the geologically
     % significant traces
     c1 = parula(numel(unique_walks(:,1)));
     % c = flipud(c);
     unique_walks = sortrows(unique_walks,2);
     
     % computing lengths of the edges
     e = table2array(g.Edges);
     for i=1:numel(e(:,1))
        e(i,3) =  compute_chain_length(e(i,1:2),xy);
     end
     e = sortrows(e,3);
     
     % creating a similar colorbar for the edge lengths
     c2 = parula(numel(e(:,1)));
     
     figure(1)
     subplot(1,2,1)
     for i=1:numel(unique_walks(:,1))
        disp(i)
        chain = unique_walks{i,1};
        plot(xy(chain,1),xy(chain,2),'Color',c1(i,:),'LineWidth',4)
        hold on
     end    
     a=adjacency(g);
     gplot(a,xy,'k')
     title('Fractures (tip-to-tip) colored by length')
     pbaspect([1 1 1])
     
     subplot(1,2,2)
      for i=1:numel(e(:,1))
        disp(i)
        chain = e(i,1:2);
        plot(xy(chain,1),xy(chain,2),'Color',c2(i,:),'LineWidth',4)
        hold on
     end    
     a=adjacency(g);
     gplot(a,xy,'k')
     title('Fracture segments colored by length')
     pbaspect([1 1 1])

     % the difference can also be depicted as a length histogram
     figure(2)
     histogram(cell2mat(unique_walks(:,end)),'FaceColor','#7E2F8E');
     hold on
     histogram(e(:,end),'FaceColor','#D95319')
     legend('Fracture Lengths(m)', 'Segment Lengths(m)')
     xlabel('Length(m)')
     ylabel('Frequency')
     
%%   7. Converting primals graphs to dual graph representations
    % in this section, a primal graph is converted into a dual graph

     clear all; close all; clc;    % clearing all workspace variables
     
     % an example graph is loaded
     load('E:\PhD\Journal_Papers\Manuscript_JSG\code_supplement\25_20_g.mat')
     load('E:\PhD\Journal_Papers\Manuscript_JSG\code_supplement\25_20_xy.mat')
     
     g.Edges.Weight =[];  % removing edge weights
     a = adjacency(g);
     
     % to compute the dual, the list of geologically significant traces is
     % needed to be computed first
     [unique_walks] = find_geologically_significant_traces(g,xy,20);
     
     Edges_p = table2array(g.Edges);  % edges in primal graph
     Edges_p = Edges_p(:,1:2);         % removing weights of edges
     for i=1:numel(unique_walks)
       disp(i) 
       chain = unique_walks{i,1};
       E = chain_to_edges(chain);
       r = findedge(g,E(:,1),E(:,2));
       if nnz(r)==numel(r)       
        Edges_p(r,3)=i;              % assign a frac_id to a set of edges
       else
        r(r==0)=[];
        Edges_p(r,3)=i;
       end    
     end    
     % Edges_p = sortrows(Edges_p,3);

     n_dual_graph_nodes = numel(unique_walks);

     % creating a sparse adjacency matrix
     a_d = sparse(n_dual_graph_nodes,n_dual_graph_nodes);

     % comparing each fracture with each other; if there are common nodes, then
     % fractures are intersecting and there is an edge between the fractures
     for i=1:n_dual_graph_nodes
       disp(i) 
       walk1 = unique_walks{i,1}; 
       [m,n] = size(walk1);
       if m<n
        walk1 = walk1';
       end
       % finding nodes that neighbor current fracture / walk nodes
       N = [];
       for j=1:numel(walk1)
          N = [N; neighbors(g, walk1(j,1))]; 
       end    

       N = setdiff(N, walk1); % subtracting walk nodes from the neighbor array

       % finding frac ids for each of N
       F = [];       % rows of Edges_p that contain neighbors to the walk
       for j=1:numel(N)
           %F = [F; find(Edges_p(:,1)==N(j,1)); find(Edges_p(:,2)==N(j,1))];
           [r,~]=find(Edges_p(:,1:2)==N(j,1));
           F = [F;r];
       end    
       F = unique(F);

       % walk indices that correspond to rows of Edges_p neighboring walk
       frac_ids = unique(Edges_p(F,3)); 
       frac_ids = setdiff(frac_ids,i);

       % removing indices 
       if ~isempty(find(frac_ids==0))
           frac_ids = setdiff(frac_ids,0);
       end    

       % final check for intersections between current walk and walks that are
       % to be added as edges in dual graph
       r=[];
       for j=1:numel(frac_ids)
           if isempty(intersect(walk1,unique_walks{frac_ids(j,1)}))
               r = [r;frac_ids(j,1)];
           end    
       end 
       frac_ids = setdiff(frac_ids,r);

       % constructing the adjacency matrix
       a_d(repmat(i,numel(frac_ids),1), frac_ids) = 1;
       a_d(frac_ids,repmat(i,numel(frac_ids),1)) = 1;
       clearvars walk N F frac_ids r
    end    

    % creating the dual graph
    g_d = graph(a_d);
    d_d = degree(g_d);
 
    % visualizing the difference between primals and duals
    figure(1)
    subplot(1,2,1)
    gplot(a,xy,'k');
    title('Primal graph')
    pbaspect([1 1 1])
    subplot(1,2,2)
    p = plot(g_d,'Layout','force','UseGravity',true);
    title('Dual graph')
    pbaspect([1 1 1])
    p.MarkerSize = d_d./5;
    p.EdgeColor = 'k';%[0 0.4470 0.7410];
    p.NodeColor = 'm';
    
    % plotting the node degree distribution of primals and duals
    d=degree(g);
    figure(2)
    subplot(1,2,1) 
    histogram(d)
    title('Primal graph')
    xlabel('Node degree dist.')
    pbaspect([1 1 1])
    subplot(1,2,2)
    histogram(d_d)
    title('Dual graph')
    xlabel('Node degree dist.')
    subplot(1,2,2)
    pbaspect([1 1 1])