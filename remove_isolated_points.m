function [XY,G] = remove_isolated_points(XY0,G0)
% this function reads a graph and identifies isolated points. Presence of
% isolated and discontinuous points can cause errors in further graph
% processing. Such nodes are identifiable by the following two checks:
%    -> node neighbor check reveals that neighbor is node itself
%    -> edgelist shows that edge is connected to itself

    j=1;
    EdgeList = table2array(G0.Edges);
    discontinuous_nodes = [];
    for i=1:numnodes(G0)
       disp(i) 
       current_node = i;
       if neighbors(G0, current_node)==current_node
           discontinuous_nodes(j,1) = i;
           edge1 = find(EdgeList(:,1)==current_node);
           disp(['Node:', num2str(i), ' is a discontinuous point']) 
           disp(['Edge:', num2str(edge1), ' will be removed from graph']);
           j=j+1;
       end    
       clearvars current_node 
    end    
    if ~isempty(discontinuous_nodes)
        G = rmnode(G0,discontinuous_nodes);
        XY = XY0;
        XY(discontinuous_nodes,:) = [];
    else    
        G = G0; XY = XY0;
    end    
end

