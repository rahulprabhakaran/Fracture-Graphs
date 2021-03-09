function [chosen_strike,quadrants] = compute_strike(Edges_G3, XY3)
% this function calculates strike of an edge w.r.t east counterclockwise
% in degrees and normalizes the angle from [0, 180]

% Caution: when strike is used to compute possible connections, be mindful
% of the edge case. i.e,  0 +/- strike threshold and 180+/- strike threshold
% for e.g., when an edge has strike of 2 degrees, and strike threshold is
% 20 degrees, then search for neighbors has to be from 2-12 and 172-180
% to be 

 strike(:,1) = rad2deg(edgeAngle([ XY3(Edges_G3(:,1),1) XY3(Edges_G3(:,1),2) ...
        XY3(Edges_G3(:,2),1) XY3(Edges_G3(:,2),2)]));
 strike(:,2) = rad2deg(edgeAngle([ XY3(Edges_G3(:,2),1) XY3(Edges_G3(:,2),2) ...
        XY3(Edges_G3(:,1),1) XY3(Edges_G3(:,1),2)]));
 strike = sort(strike,2);
 
 quadrants = zeros(numel(strike(:,1)),1);  % listing of quadrants that edges lie


 for i=1:numel(strike(:,1))
     if strike(i,1) < 90
        quadrants(i,1) = 13;  % quadrants 1 and 2
     else
        quadrants(i,1) = 24;  % quadrants 3 and 4
     end    
 end
 
 chosen_strike = min(strike(:,1),strike(:,2));
 


% strike1 = rad2deg(edgeAngle([ XY3(Edges_G3(:,1),1) XY3(Edges_G3(:,1),2) ...
%     XY3(Edges_G3(:,2),1) XY3(Edges_G3(:,2),2)]));
% strike2 = rad2deg(edgeAngle([ XY3(Edges_G3(:,2),1) XY3(Edges_G3(:,2),2) ...
%     XY3(Edges_G3(:,1),1) XY3(Edges_G3(:,1),2)]));
% 
% chosen_strike = min(strike1,strike2);


end

