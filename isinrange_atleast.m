function [out] = isinrange_atleast(value, minval, maxval)
% checks if a value, or atleast one value in a row vector of values lies within an interval 
% returns '0' or '1' depending on whether atleast one value lie within the interval

chk_vector1 = value > minval;
chk_vector2 = value < maxval;
chk_vector3 = chk_vector1.*chk_vector2;

if nnz(chk_vector3)==0
    out = 0;
else
    out = 1;
end    

end

