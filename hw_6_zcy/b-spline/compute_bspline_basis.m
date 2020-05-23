function [N] = compute_bspline_basis(degree,n,knot,t)
% compute b-spline basis for Ni,p(u) for each control point Pi, i=[0,n]
% input:
% n: n+1 control points
% p: degree
% knots: m+1 knots, m = n+p+1
% t: knot t in knot span [knot(k),knot(k+1))
N = zeros(n+1,1);
m = max(size(knot))-1;
p = degree;
% first rule out special cases
if t == knot(1)
    N(1) = 1;
    return
elseif t == knot(m+1)
    N(n+1) = 1;
    return
end
% now t is btwn knot(1) and knot(m+1)
% let t be in knot span [knot(k),knot(k+1))
% find out k
k = 1;
for kk = 1:m+1
    if knot(kk) > t
        k = max(kk-1,1);
        break;
    end
end
N(k) = 1;
% I added some if conditions to include no clamped curve
for d = 1:p % degree d goes from 1 to p
    if k > d
        N(k-d) = (knot(k+1)-t)/(knot(k+1)-knot(k-d+1))*N(k-d+1); % right(south-west corner) term only
    end
    if k >= d
        for i = k-d+1:k-1 % internal terms
            N(i) = (t-knot(i))/(knot(i+d)-knot(i))*N(i)+(knot(i+d+1)-t)/(knot(i+d+1)-knot(i+1))*N(i+1);
        end
    end
    if k <= m-d
        N(k) = (t-knot(k))/(knot(k+d)-knot(k))*N(k); % left (north-west corner) term only
    end
end

end

