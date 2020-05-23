function [c_pts_opt] = trim_c_pts_opt(c_pts_opt,p,init_type)
% get rid of extended control points to make sure no densely points at
% start/end
% p: degree
if init_type == 0 % pos constraint only
    c_pts_opt = c_pts_opt(p:size(c_pts_opt,1)-(p-1),:);
elseif init_type == 1 % full start & end pos
    c_pts_opt = c_pts_opt(1:size(c_pts_opt,1)-(p-1),:);
else % full start & full end
end 
end

