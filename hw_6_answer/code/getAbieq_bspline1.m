function [Aieq, bieq] = getAbieq_bspline1(n, p, corridor_range, v_max, a_max)
% generate bspline equality constraints matrix
% n: n+1 control points
% p: degree of bspline
% note we only need to optimize middle n+1-2*p control points
% the first & last p control points are determined by start & end
% conditions
% first calculate knots
knots = deboor_knot(p,n,2); % we need each curve seg to be clamped at both ends
dknots = knots(2:max(size(knots))-1);
ddknots = dknots(2:max(size(dknots))-1);
dp = p-1;
ddp = dp-1;
dn = n-1;
ddn = dn-1;
C = getC(n,p,knots);
dC = getC(dn,dp,dknots);
ddC = getC(ddn,ddp,ddknots);
%#####################################################
% STEP 3.2.1 p constraint
Aieq_p = [eye(n+1);-eye(n+1)];
Aieq_p_upper = eye(n+1);
Aieq_p_upper = Aieq_p_upper(p+1:n-p+1,:);
Aieq_p_lower = -eye(n+1);
Aieq_p_lower = Aieq_p_lower(p+1:n-p+1,:);
Aieq_p = [Aieq_p_upper;Aieq_p_lower];
bieq_p = zeros((n+1-2*p)*2,1);
for i = 1:n+1-2*p
    bieq_p(i) = corridor_range(i+p,2); % max limit
    bieq_p(i+n+1-2*p) = -corridor_range(i+p,1); % min limit
end
%#####################################################
% STEP 3.2.2 v constraint  
Aieq_v = [C;-C];
Aieq_v_upper = C;
Aieq_v_upper = Aieq_v_upper(dp+1:dn-dp+1,:);
Aieq_v_lower = -C;
Aieq_v_lower = Aieq_v_lower(dp+1:dn-dp+1,:);
Aieq_v = [Aieq_v_upper;Aieq_v_lower];
b_half_size = dn+1-2*dp;
bieq_v = zeros(b_half_size*2,1);
for i = 1:b_half_size
    bieq_v(i) = v_max; % max limit
    bieq_v(i+b_half_size) = v_max; % min limit
end
%#####################################################
% STEP 3.2.3 a constraint   
Aieq_a = [dC*C;-dC*C];
Aieq_a_upper = dC*C;
Aieq_a_upper = Aieq_a_upper(ddp+1:ddn-ddp+1,:);
Aieq_a_lower = -dC*C;
Aieq_a_lower = Aieq_a_lower(ddp+1:ddn-ddp+1,:);
Aieq_a = [Aieq_a_upper;Aieq_a_lower];
b_half_size = ddn+1-2*ddp;
bieq_a = zeros(b_half_size*2,1);
for i = 1:b_half_size
    bieq_a(i) = a_max; % max limit
    bieq_a(i+b_half_size) = a_max; % min limit
end
%#####################################################
% combine all components to form Aieq and bieq   
Aieq = [Aieq_p; Aieq_v; Aieq_a];
bieq = [bieq_p; bieq_v; bieq_a];


end