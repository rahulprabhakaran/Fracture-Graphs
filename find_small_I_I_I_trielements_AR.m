function [small_I_I_I_trielements_AR] = find_small_I_I_I_trielements_AR(trielement_list,area_threshold,...
    aspect_ratio_threshold,XY)
% this function extracts I_I_I trielements subject to a minimum area
% condition and an aspect ratio criterion and returns a list 

j=1;
for i=1:numel(trielement_list(:,1))
 disp(i)   
 if (trielement_list(i,4)==1 && trielement_list(i,5)==1 && trielement_list(i,6)==1)==1
     idx_I_I_I(j,1)=i;
     I_I_I_trielements(j,1:3) = trielement_list(i,1:3);
     I_I_I_trielements(j,4) = i; 
     I_I_I_trielements(j,5)=polyarea(XY(I_I_I_trielements(j,1:3),1),XY(I_I_I_trielements(j,1:3),2));
     j=j+1;
 end     
end


small_I_I_I_trielements_idx= find(I_I_I_trielements(:,5)<area_threshold);
small_I_I_I_trielements = I_I_I_trielements(small_I_I_I_trielements_idx,:);


for i=1:numel(small_I_I_I_trielements(:,1))
 small_I_I_I_trielements(i,6)= tri_aspect_ratio(XY(small_I_I_I_trielements(i,1:3),1), XY(small_I_I_I_trielements(i,1:3),2)); 
end

% adjusting for aspect ratio, only selecting those trielements with aspec
small_I_I_I_trielements_AR= small_I_I_I_trielements(find(small_I_I_I_trielements(:,6)<aspect_ratio_threshold),:);
%  [small_I_I_I_trielements_AR] = find_I_I_I_edges(Combined_Edges_cleaned_duplicate_removed,small_I_I_I_trielements_AR);

end

