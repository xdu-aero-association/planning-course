function [c_pts1] = get_extended_c_pts(c_pts,p,init_type)
% extend control points to make up for front/end kth order continuity
% p: degree
c_f = ones(p-1,2).*c_pts(1,:);
c_b = ones(p-1,2).*c_pts(end,:);
if init_type == 0 % pos constraint only
    c_pts1 = [c_f;c_pts;c_b];
elseif init_type == 1 % full start & end pos
    c_pts1 = [c_pts;c_b];
else % full start & full end
    c_pts1 = c_pts;
end 

end

