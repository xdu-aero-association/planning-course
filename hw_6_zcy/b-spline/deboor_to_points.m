function pts = deboor_to_points(degree, n, c_pts, knot, dense, type)
% calculate points on b-spline
% degree: degree of b-spline
% n for n+1 c_pts
% c_pts: control points
% knot: array of knot positions
% dense: sample num on each curve segment
% type: 0 for deboor algorithm, 1 for standard def.

cnt = 0;
pts = zeros((n-degree+1)*dense,2);
for i = degree+1:n+1
    if knot(i+1)>knot(i)
        for ii = 0:dense-1
            t = knot(i) + ii*(knot(i+1)-knot(i))/dense;
            if type == 0
                pts(cnt+1,:) = deboor(degree, c_pts, knot, t, i);
            else
                N_set = compute_bspline_basis(degree,n,knot,t);
                pts(cnt+1,1) = sum(N_set.*c_pts(:,1));
                pts(cnt+1,2) = sum(N_set.*c_pts(:,2));
            end
            cnt = cnt+1;
        end
    end
    
end
fprintf(1,'pts calculation done\n')
end

