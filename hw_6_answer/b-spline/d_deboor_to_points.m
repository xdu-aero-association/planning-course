function d_pts = d_deboor_to_points(degree, n, c_pts, knot, dense)
% calculate 1st order deri points on b-spline
% degree: degree of orginal b-spline
% n for original n+1 c_pts
% c_pts: original control points
% knot: array of original knot positions
% dense: sample num on each curve segment
% derivative control points & knots
[d_c_pts, dknot] = get_deri_c_pts(degree, n, c_pts, knot);

cnt = 0;
d_pts = zeros((n-degree+1)*dense,2);
for i = degree:n
    if dknot(i+1)>dknot(i)
        for ii = 0:dense-1
            t = dknot(i) + ii*(dknot(i+1)-dknot(i))/dense;
            N_set = compute_bspline_basis(degree-1,n-1,dknot,t);
            d_pts(cnt+1,1) = sum(N_set.*d_c_pts(:,1));
            d_pts(cnt+1,2) = sum(N_set.*d_c_pts(:,2));
            cnt = cnt+1;
        end
    end
    
end
fprintf(1,'d_pts calculation done\n')
end

