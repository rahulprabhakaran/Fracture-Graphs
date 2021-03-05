function [Path_Shape] = build_shp_struct(cell_array,XY3,xmin,ymin,flag)
% this function builds a shape structure for polyline fractures or
% points

% flag = 1, for 'Line'
% flag = 2, for 'Point'
% flag = 3, for 'Polygon'

if flag == 1
    for i=1:numel(cell_array(:,1))    
     disp(i)   
     [Path_Shape(i).Geometry] = 'Line';
     [Path_Shape(i).X] = XY3(cell_array{i,1},1) + xmin;
     [Path_Shape(i).Y] = XY3(cell_array{i,1},2) + ymin;
     [Path_Shape(i).ID] = i;
     % [Path_Shape3(i).Azi] = path_no_overlaps{i,6};
    end
end

if flag == 2
 for i=1:numel(cell_array(:,1))    
  disp(i)   
  [Path_Shape(i).Geometry] = 'Point';
  [Path_Shape(i).X] = XY3(cell_array(i,1),1) + xmin ;
  [Path_Shape(i).Y] = XY3(cell_array(i,1),2) + ymin;
  [Path_Shape(i).ID] = i;
 end
end

if flag == 3
    for i=1:numel(cell_array(:,1))    
     disp(i)   
     [Path_Shape(i).Geometry] = 'Polygon';
     [Path_Shape(i).X] = cell_array{i,1}(:,1) + xmin;
     [Path_Shape(i).Y] = cell_array{i,1}(:,2) + ymin;
     [Path_Shape(i).Area] = cell_array{i,2}(:,1);
     [Path_Shape(i).ID] = i;
     % [Path_Shape3(i).Azi] = path_no_overlaps{i,6};
    end
end

end

