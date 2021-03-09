function walk = find_possible_edges3(G3,node_s,chain,XY3,strike_threshold)
% modification of find_possible_edges2; since the strike is normalized from
% [0,180], this can lead to loss of connectivity when the chain encounters
% edges that are near horizontal, i.e., 0-10 deg and 170-180 deg. The
% change w.r.t older function is that now the strike range to be checked is 
% enlarged for these edge cases. The second change was made to rectify the
% problem of walks that travelled backwards. The acute angle check between
% chain and possible edges is now done prior to checking whether possible
% edges are within range.

% the acute angle check is done between chain [node_s t] and possible edges
% [node_s t]. This ordering is crucial as the lineAngle function from 
% geom2D toolbox gives a different result when the nodes are not specified
% correctly.

% initializing to enter the while loop
possible_edges_final = 1;

% first edge of the walk is the input edge
walk = chain;

% counter
j=1;

while ~isempty(possible_edges_final)
     % calculating strike of current chain
     [chain_strike,quadrants] = compute_strike(chain, XY3);
     
     max_strike_range = chain_strike + strike_threshold;
     min_strike_range = chain_strike - strike_threshold;
              
     % enumerating the possible number of edges for current chain
     possible_edges = [neighbors(G3, node_s) repelem(node_s,numel(neighbors(G3, node_s)))'];
       
     % removing the starting edge from the possible edges list  
     possible_edges = setdiff(possible_edges,chain,'rows');
     possible_edges = setdiff(possible_edges,fliplr(chain),'rows');
     [m,~] = size(possible_edges);
     
     % strike for all possible chains
     possible_edges(:,3) = compute_strike(possible_edges, XY3);
     
     % replicating the origin chain 'm' times and converting to lines  
     % A1 -> line starting from node_s
     A1 = edgeToLine([repelem(XY3(node_s,1),m)' repelem(XY3(node_s,2),m)' ...
         repelem(XY3(setdiff(chain,node_s),1),m)' repelem(XY3(setdiff(chain,node_s),2),m)']);
     
     % A2 -> line ending in node_s
     A2 = edgeToLine([ repelem(XY3(setdiff(chain,node_s),1),m)' repelem(XY3(setdiff(chain,node_s),2),m)'...
         repelem(XY3(node_s,1),m)' repelem(XY3(node_s,2),m)']);    
     % angle_A = [min(rad2deg(lineAngle(A1)),rad2deg(lineAngle(A2))) max(rad2deg(lineAngle(A1)),rad2deg(lineAngle(A2)))];
     
     % B1 -> possible edges lines ending in node_s
     B1 = edgeToLine([XY3(possible_edges(:,1),1) XY3(possible_edges(:,1),2) ...
         XY3(possible_edges(:,2),1) XY3(possible_edges(:,2),2)]);

     % B2 -> possible edges lines starting in node_s
     B2 = edgeToLine([XY3(possible_edges(:,2),1) XY3(possible_edges(:,2),2) ...
         XY3(possible_edges(:,1),1) XY3(possible_edges(:,1),2)]);
     % angle_B = [min(rad2deg(lineAngle(B1)),rad2deg(lineAngle(B2))) max(rad2deg(lineAngle(B1)),rad2deg(lineAngle(B2)))];
     
     % when chain strike is close to horizontal, then the max_strike_range
     % is greater than 180, and edges with very low angles also need to
     % be considered as possible edges; the strike of possible neighboring
     % edges from 0 to max_strike_range-180 is increased above 180
     if max_strike_range > 180 
        [r1,~] = find(possible_edges(:,3) < max_strike_range-180);
        possible_edges(r1,3) = possible_edges(r1,3)+180;
     end  
     
     % when chain strike is close to horizontal, then min_strike_range can
     % become less than 0, and edges with very high angles also need to
     % be considered as possible edges; the strike of possible neighboring
     % edges from 180 + min_strike_range is reduced by 180
     if min_strike_range <0
        [r2,~] = find(possible_edges(:,3) > 180 + min_strike_range);
        possible_edges(r2,3) = 180 - possible_edges(r2,3);
     end    
            
     % angle difference between origin chain and possible edges, always
     % positive and between 0 to 180
     % lineAngle returns the directed angle between two lines; in
     % in counter-clockwise direction.

     possible_edges(:,4) = abs(rad2deg(lineAngle(A1)) - rad2deg(lineAngle(B2)));

     %      possible_edges(:,4) = min(rad2deg(lineAngle(B,A)),rad2deg(lineAngle(A,B))); 
     
     % filtering for only those edges that diverge forward from
     % origin node; r_s1 -> edges that diverge away
     [r_s1,~]= find(possible_edges(:,4)>90);
     possible_edges = possible_edges(r_s1,:);
     clearvars r_s1
      
     % finding edges within the range of strike; 'r' row indices of edges
     % within range of strike
     [r,~]=find( (min_strike_range <= possible_edges(:,3)) & ...
                  (possible_edges(:,3) <= max_strike_range)); 
      
       
     if ~isempty(r) 
         
         % filtering out edges that are in range of strike threshold
         possible_edges = possible_edges(r,:);
        
         % condition when there are more than one choice
         if numel(r) > 1 
             
             % filtering for only those edges that diverge forward from
             % origin node; r_s1 -> edges that diverge away
%              [r_s1,~]= find(possible_edges(:,4)>90);
%              possible_edges = possible_edges(r_s1,:);

             
             % filtering edges that have least deviation from chain
             possible_edges(:,5) = abs(possible_edges(:,3)-chain_strike);
             [~,r_s2]= min(possible_edges(:,5));
             
             
             % retaining only one edge after all eliminations
             possible_edges_final = possible_edges(r_s2,1:2);
             clearvars r_s2
           
         % condition when there is only one possible choice     
         else
             possible_edges_final = possible_edges(:,1:2);
         end  
         clearvars r
         % updating source node and chain  
         tmp = node_s;
         node_s = setdiff(possible_edges_final,node_s);
         chain = [node_s tmp];    
         walk = [node_s walk];
         
         % this hard exit is so that chain does not repeat endlessly
         if length(walk)~=length(unique(walk))
           walk(1,1)=0;  
           walk = nonzeros(walk)';
           possible_edges_final = [];
         end   
         
     else
         possible_edges_final = [];
     end
     
     j=j+1;    
end 
end

