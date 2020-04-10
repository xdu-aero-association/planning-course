function [dist] = cal_dist(p1,p2)
%CAL_DIST 此处显示有关此函数的摘要
%   此处显示详细说明
dp = p1 - p2;
dist = sqrt(dp(1)^2+dp(2)^2);
end

