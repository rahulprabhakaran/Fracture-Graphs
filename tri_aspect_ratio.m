function [AR] = tri_aspect_ratio(x,y)
% this function computes the aspect ratio of a trielement
% inputs are coordinates of the three corners
a = sqrt((x(1,1)-x(2,1))^2 + (y(1,1)-y(2,1))^2);
b = sqrt((x(2,1)-x(3,1))^2 + (y(2,1)-y(3,1))^2);
c = sqrt((x(3,1)-x(1,1))^2 + (y(3,1)-y(1,1))^2);

s = (a+b+c)/2;

AR = a*b*c / 8 / (s-a) / (s-b) / (s-c);
end

