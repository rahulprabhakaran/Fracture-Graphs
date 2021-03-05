function [kinked_M_nodes] = find_kinked_M_Nodes(XY3,M_nodes_idx,angle_1,angle_2,G3)
% within the chains of the graph, there are degree-2 nodes where there is
% a sharp kink. These nodes and the chains corresponding to them need to
% be identified

% not all of the Combined_Edges contain degree-2 nodes. Identifying the 
% index of these chains

% j=1;
% for i=numel(Combined_Edges(:,1))
%    chain =  Combined_Edges{i,3};
%    if ~isempty(find(contains(chain,'2')))
%        idx(j,1)=i;
%        j=j+1;
%    end  
%    clearvars chain
% end    

% only those chains with degree 2 nodes are considered
% Combined_Edges_with_2_deg_nodes = Combined_Edges(idx,:);
% 
% for i=1:numel(Combined_Edges_with_2_deg_nodes(:,1))
%  disp(i)   
%  Chain = Combined_Edges_with_2_deg_nodes{i,1};
%  Chain_XY = XY_new(Chain,:);
%  Angles_Chain = Angles2D( [Chain_XY(1,:)  Chain_XY(2,:)],1);
%  Combined_Edges_reshaped_set2{i,2} = degree(G2,Combined_Edges_reshaped_set2{i,1});
%  Combined_Edges_reshaped_set2{i,3}=Angles_Chain;
%  clearvars Chain Angles_Chain 
% end 

j=1;
for i=1:numel(M_nodes_idx)
 
 disp(i)   
 current_node = M_nodes_idx(i,1);
 if neighbors(G3, current_node)~=current_node
     n1n2 = neighbors(G3,current_node);

     branch1 = [XY3(current_node,1:2) XY3(n1n2(1,:),1:2)];
     branch2 = [XY3(n1n2(2,:),1:2) XY3(current_node,1:2)];
     branch1_flip = [XY3(n1n2(1,:),1:2) XY3(current_node,1:2)];
     branch2_flip = [XY3(current_node,1:2) XY3(n1n2(2,:),1:2)];
     ang_branch1_min = min(rad2deg(edgeAngle(branch1)),rad2deg(edgeAngle(branch1_flip)) );
     ang_branch1_max = max(rad2deg(edgeAngle(branch1)),rad2deg(edgeAngle(branch1_flip)) );
     ang_branch2_min = min(rad2deg(edgeAngle(branch2)),rad2deg(edgeAngle(branch2_flip)) );
     ang_branch2_max = max(rad2deg(edgeAngle(branch2)),rad2deg(edgeAngle(branch2_flip)) );
    %  ang_branch1 = Angles2D(branch1,1);
    %  ang_branch2 = Angles2D(branch2,1);
     angle_diff_min = abs(ang_branch1_min - ang_branch2_min) ;
    %  angle_diff_max = abs(ang_branch1_max - ang_branch2_max) ;

     if angle_diff_min>angle_1 && angle_diff_min < angle_2    
    %      if ~isempty(find(contains(Combined_Edges(:,4),num2str([n1n2(1,:);current_node;n1n2(2,:)]'),'IgnoreCase',true)))
             kinked_M_nodes(j,1) = current_node;
    %          Combined_Edges_kinked_idx(j,1)=find(contains(Combined_Edges(:,4),num2str([n1n2(1,:);current_node;n1n2(2,:)]'),'IgnoreCase',true));
    %      end
         j=j+1;
     end  

    %  Chain_XY(1,:) = branch2(1,1:2); 
    %  Chain_XY(2,:) = branch2(1,3:4); 
    %  Chain_XY(3,:) = branch1(1,3:4);
    %  figure(3)
    %  plot(Chain_XY(:,1),Chain_XY(:,2),'*r');
 else
    disp('Node is a discontinuous point')    
 end    
     clearvars n1n2 current_node branch1 branch2 ang_branch1 ang_branch2 angle_diff_min
end    

% Combined_Edges_kinked_idx(Combined_Edges_kinked_idx==0)=[];
kinked_M_nodes(kinked_M_nodes==0)=[];

end

