# Fracture-Graphs

This repository contains MATLAB code and functions that forms the Code supplement to the manuscript titled, "Prabhakaran et al (2021), Large-scale natural fracture network patterns: Insights
from automated mapping in the Lilstock (Bristol Channel) limestone outcrops". The Code_Supplement.m script showcases the use of graph-based methods to read-in shapefiles and manipulate them using graph algorithms. Functions from the Geom2D Toolbox by David Legland (see https://nl.mathworks.com/matlabcentral/fileexchange/7844-geom2d). MATLAB Mapping toolbox is required to read and write shapefiles. The Data folder contains *.mat & *.shp files which are used as examples to illustrate the working of the code. The data may be downloaded and paths set accordingly by the user.

-------------CONTENTS------------------------------
1. Loading a shapefile and converting it to a graph: 2D fracture networks are often saved in the ESRI shapefile format in the form of polylines. These are converted to graph structures with edges and nodes. Spatial positioning information corresponding to each node makes it a spatial graph
2. Fixing three types of common topological discontinuites: Automatic detection lead to certain topological discontinuities that needs to be resolved so that network metrics can be applied.
2.1 Type-1 discontinuity fix: fracture with a degree-3 node is incorrectly represented as a degree-2 node in close proximity to a degree-1 node
2.2 Type-2 discontinuity fix: three degree-1 nodes lie in close proximity
2.3 Type-3 discontinuity fix: two degree-2 nodes lie in close proximity with sharp angles
3. Resolving artificial fragmentation in shapefile attribute tables: Vectorized traces may be topologically connected but fragmented and this artificial fragmentation is resolved.
4. Resolving Stepouts within fracture graphs
4.1 Resolving stepouts using merge operation: Merging a stepout reduces two degree-3 nodes into a degree-4 node 
4.2 Resolving stepouts using flatten operation: Stepouts based on length threshold and based on two additional inputs of angle and orthogonal range are flattened
5. Straightening fracture segments (removing degree-2 nodes): when converting a shapefile to a graph, there are a large number of degree-2 nodes. These represent the natural   sinuousity of fracture traces. Within the graph representation, they increase the size of the graph and on some occasion we would like to remove these degree-2nodes while retaining the overall graph structure
6. Converting fracture graphs to geologically-significant fracture traces: A fracture from tip-to-tip can be identified from a series of graph edges provided they are continuous in terms of edge angle. An input threshold strike angle is set and the function identifies a set of fracture edges that are within the threshold
7. Converting primals graphs to dual graph representations: In the primal graph representation, intersections between the fractures and fracture tips are regarded as nodes with the edges being fracture segments between the nodes. This primal representation can be converted into a dual graph where a fracture from tip-to-tip represent nodes and intersection between fractures are edges

functions list:
build_shp_struct, chain_to_edges, check_bifurcations, check_if_line_graph, collapse_stepouts, compute_chain_length, compute_chain_strike2, compute_strike, create_GraphEdges, create_subgraph, edges_to_chain, find_chain, find_geologically_significant_traces, find_I_M_possible_edges, find_kinked_M_Node_connections, find_kinked_M_Nodes, find_small_I_I_I_trielements_AR, find_stepout_motif, find_walk, fix_graph_type_5_disconnection, fix_graph_type_8_disconnection,
flatten_stepouts, isinrange_atleast, Lengths2D, resolve_artificial_fragmentation, Shape_to_Data, straighten_graph, tri_aspect_ratio, update_after_collapse
