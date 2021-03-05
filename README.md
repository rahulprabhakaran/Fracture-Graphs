# Fracture-Graphs

This repository contains MATLAB code and functions that forms the Code supplement to the manuscript titled, "Prabhakaran et al (2021), Large-scale natural fracture network patterns: Insights
from automated mapping in the Lilstock (Bristol Channel) limestone outcrops". The Code_Supplement.m script showcases the use of graph-based methods to read-in shapefiles and manipulate them using graph algorithms.

The code uses certain functions from the Geom2D Toolbox by David Legland (see https://nl.mathworks.com/matlabcentral/fileexchange/7844-geom2d).

-------------CONTENTS------------------------------
1. Loading a shapefile and converting it to a graph
2. Fixing three types of common topological discontinuites
2.1 Type-1 discontinuity fix
2.2 Type-2 discontinuity fix
2.3 Type-3 discontinuity fix
3. Resolving artificial fragmentation in shapefile attribute tables
4. Resolving Stepouts within fracture graphs
4.1 Resolving stepouts using merge operation
4.2 Resolving stepouts using flatten operation
5. Straightening fracture segments (removing degree-2 nodes)
6. Converting fracture graphs to geologically-significant fracture traces
7. Converting primals graphs to dual graph representations

functions list:
Shape_to_Data, build_shp_struct, chain_to_edges, collapse_stepouts,compute_chain_length, create_GraphEdges, find_I_M_possible_edges,find_geologically_significant_traces, find_kinked_M_Node_connections,find_kinked_M_Nodes, find_small_I_I_I_trielements_AR, find_stepout_motif,fix_graph_type_5_disconnection, fix_graph_type_8_disconnection,
flatten_stepouts, resolve_artificial_fragmentation, straighten_graph
