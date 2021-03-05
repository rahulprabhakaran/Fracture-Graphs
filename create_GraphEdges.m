function [EndNodes,XY] = create_GraphEdges(Fractures, xmin, ymin)
% This function generates an EdgeTable and NodeTable (in table) format
% which are inputs to create graph. An XY matrix (N x 2) is also returned   
% which can be used to reference the graph. The Fractures matrix (N x 4)
% consists of the fracture segment start and end coordinates
% Rahul Prabhakaran (2019)
% 
Fractures(:,1)=Fractures(:,1)-xmin;
Fractures(:,3)=Fractures(:,3)-xmin;
Fractures(:,2)=Fractures(:,2)-ymin;
Fractures(:,4)=Fractures(:,4)-ymin;

% Node matrix for all nodes
% tic
All_Nodes = [ [Fractures(:,1:2) Fractures(:,5)]; [Fractures(:,3:4) Fractures(:,5)] ];
% All_Nodes = [Fractures(:,1:2); Fractures(:,3:4)];
% toc


% removing all duplicates
% tic
[All_Nodes_uniq,~,~] = unique(All_Nodes,'rows');
[z,~]=size(All_Nodes_uniq);
All_Nodes_Concat{z,1} ={};
Fracture_storage{z,1} ={};
string_length{z,1} ={};
% toc

tol = 10E15;

% Loop 1: this loop creates a cell array with string entries that concatenates the X
% and Y coordinate so that it is a unique identifier....1.53 secs
% tic
for i=1:z
   disp('Loop1')      
   disp(i)
   All_Nodes_Concat{i,1} = strcat(num2str(round(All_Nodes_uniq(i,1),10)*tol),num2str(round(All_Nodes_uniq(i,2),10)*tol));
   Fracture_storage{i,1} = All_Nodes_uniq(i,3);
   string_length{i,1} = [length(num2str(round(All_Nodes_uniq(i,1),10)*tol)) length(num2str(round(All_Nodes_uniq(i,2),10)*tol))];
   %clearvars A 
end
% toc

% removing duplicates in concatenated cell array
[All_Nodes_Concat_uniq,id,~] = unique(All_Nodes_Concat);

% assigning the polyline ID to the All_Nodes_Concat_uniq cell array
for i=1:length(All_Nodes_Concat_uniq)
 Fracture_storage_uniq{i,1} = Fracture_storage{id(i,1),1};
end

% Loop 2: creating a string length array. We will use this later for the XY matrix
% tic
for i=1:length(All_Nodes_Concat_uniq)
   disp('Loop2')      
   disp(i)
   string_length_uniq(i,:) = string_length{id(i),1};
end
% toc

% creating a node table. Not necessary to construct graph. only for comparison 
NodeTable = array2table(All_Nodes_Concat_uniq);

% Loop 3: this loop creates a cell array of string entries that concatenates the X
% Y coordinate of the edge start and edge end nodes
% tic
[z,~]=size(Fractures);
Edges_Concat{z,2}={};
for i=1:z
  disp('Loop3')      
  disp(i)
  Edges_Concat{i,1} = strcat(num2str(round(Fractures(i,1),10)*tol),num2str(round(Fractures(i,2),10)*tol));
  Edges_Concat{i,2} = strcat(num2str(round(Fractures(i,3),10)*tol),num2str(round(Fractures(i,4),10)*tol));
end
% toc

% Loop 4: now the edges can be numbered by comparing with the unique list of
% nodes. Loop speeded up around 18 times using find(contains(....
%  tic

% find(contains( )) did not work on some datasets...converting the character
% cell arrays to string aarys and using find(contains( )) on the string
% arrays
All_Nodes_Concat_uniq_str = convertCharsToStrings(All_Nodes_Concat_uniq);


[z,~]=size(Fractures);
EndNodes{z,1}={};
for i=1:z   
    disp('Loop4')      
    disp(i)
% %     tic
%      A = strfind(All_Nodes_Concat_uniq,Edges_Concat{i,1});
% %     toc         
% %     tic
%      B = strfind(All_Nodes_Concat_uniq,Edges_Concat{i,2});
% %     toc    
% %     tic
%      EndNodes{i,1} = {find(~cellfun(@isempty,B)),find(~cellfun(@isempty,A))};     
% %     toc     
%      clearvars A B    
%     tic
    %idx_A = find(contains(All_Nodes_Concat_uniq,Edges_Concat{i,1}));
%     toc  
%     tic
    %idx_B = find(contains(All_Nodes_Concat_uniq,Edges_Concat{i,2}));
%     toc    
%     tic
    %EndNodes{i,1} = {idx_B,idx_A};     
%     toc  
%       EndNodes{i,1} = {find(contains(All_Nodes_Concat_uniq,Edges_Concat{i,2})),find(contains(All_Nodes_Concat_uniq,Edges_Concat{i,1}))};
      EndNodes{i,1} = {find(contains(All_Nodes_Concat_uniq_str,convertCharsToStrings(Edges_Concat{i,2}))),...
          find(contains(All_Nodes_Concat_uniq_str,convertCharsToStrings(Edges_Concat{i,1})))};
      EndNodes{i,2} = Fractures(i,5); 
%      C2 = contains(All_Nodes_Concat_uniq,Edges_Concat{i,2});
%      C1 = contains(All_Nodes_Concat_uniq,Edges_Concat{i,1});
%      [nonZeroRows2, ~] = find(C2 ~= 0);
%      [nonZeroRows1, ~] = find(C1 ~= 0);
%      EndNodes{i,1} = {nonZeroRows2,nonZeroRows1};
    
    %clearvars idx_A idx_B 
end
%    toc

% Loop 5: creating the XY matrix
% tic
[z,~]=size(All_Nodes_Concat_uniq);
XY=zeros(z,2);
for i=1:z
  disp('Loop5')      
  disp(i)
  concat_string= All_Nodes_Concat_uniq{i,1};  
  if string_length_uniq(i,1) + string_length_uniq(i,2) ~= length(concat_string)
   print("Mismatch in string size.")   
  else
 
   XY(i,1) = str2double(concat_string(1:string_length_uniq(i,1)));


   XY(i,2) = str2double(concat_string(string_length_uniq(i,1)+1:end));

   %XY(i,:)
  end
  clearvars  concat_string 
%   disp(i)
end
XY = XY./tol;
% toc



% converting the cell array into a table
%EndNodes = cell2table(EndNodes);
%EndNodes = sortrows(EndNodes);

%All_Nodes_uniq = unique(All_Nodes,'rows');
%NodeTable = [1:1:length(All_Nodes_uniq)]';
%NodeTable = array2table(NodeTable);


end

