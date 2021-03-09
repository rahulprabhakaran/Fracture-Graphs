function [chain] = edges_to_chain(sub_G3)
% convert the edges of a line graph into a chain with correct order 
% input argument is the sub-graph; the edges may be in disorder row
% or columnwise
  if ~isempty(sub_G3.Edges) && ~isempty(sub_G3.Nodes)
    % extracting edges from the sub-graph
    E = str2double(table2array(sub_G3.Edges));
    
    % creating an anonymized sub-graph using same edges
    sub_G3_anon = graph();
    sub_G3_anon = addedge(sub_G3_anon,string(E(:,1)),string(E(:,2)));  
    
    % edges of anonymized graph with node degrees
    E1 = str2double(table2array(sub_G3_anon.Edges));
    E1(:,3) = degree(sub_G3_anon,string(E1(:,1)));
    E1(:,4) = degree(sub_G3_anon,string(E1(:,2)));
    
    %--------------------------------------------------------------------%
    % the first edge must contain deg 1 node in 1st column
    [r1,~] = find(E1(:,3)==1);   % finding 's' or 't' nodes in edge starts
    if ~isempty(r1)              
      if numel(r1)>1             % when two choices, just choose minimum
        r1 = min(r1);
      end
      % interchanging row with 's' or 't' edge as 1st row
      tmp = E1(1,:);           
      E1(1,:) = E1(r1,:);
      E1(r1,:)= tmp;     
    else                          
      [r1,~] = find(E1(:,4)==1); % finding 's' or 't' nodes in edge tails
      if numel(r1)>1
        r1 = min(r1);            % when two choices, just choose minimum
      end
      % interchanging row with 's' or 't' edge as 1st row
      E1(r1,1:2) = fliplr(E1(r1,1:2)); % flipping identified row
      tmp = E1(1,:);
      E1(1,:) = E1(r1,:);              % replacing 1st row with flipped 
      E1(r1,:)= tmp;                   % 1st row re-assigned
    end
    E1(:,3:4)=[];               % removing cols containing node degrees
    
    %--------------------------------------------------------------------%
    % the first edge must contain deg 1 node in 1st column
    chain(:,1)=E1(1,:)';
    E1(1,:)=0;
    for i=1:numel(E1(:,1)-1)
    %  disp(i)
        [r,c] = find(E1==chain(end,1));  
        chain = [chain; E1(r,setdiff([1 2],c))];
        E1(r,:)=0;
        clearvars r c
    end    
    chain = chain';
  else
    chain =[];  
  end  
    
end

