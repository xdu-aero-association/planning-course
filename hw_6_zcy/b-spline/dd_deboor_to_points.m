function dd_pts = dd_deboor_to_points(degree, n, c_pts, knot, dense)
% calculate 1st order deri points on b-spline
% degree: degree of orginal b-spline
% n for original n+1 c_pts
% c_pts: original control points
% knot: array of original knot positions
% dense: sample num on each curve segment
% derivative control points & knots
[d_c_pts, dknot] = get_deri_c_pts(degree, n, c_pts, knot);
[dd_c_pts, ddknot] = get_deri_c_pts(degree-1, n-1, d_c_pts, dknot);

cnt = 0;
dd_pts = zeros((n-degree+1)*dense,2);
for i = degree-1:n-1
    if ddknot(i+1)>ddknot(i)
        for ii = 0:dense-1
            t = ddknot(i) + ii*(ddknot(i+1)-ddknot(i))/dense;
            N_set = compute_bspline_basis(degree-2,n-2,ddknot,t);
            dd_pts(cnt+1,1) = sum(N_set.*dd_c_pts(:,1));
            dd_pts(cnt+1,2) = sum(N_set.*dd_c_pts(:,2));
            cnt = cnt+1;
        end
    end
    
end
fprintf(1,'d_pts calculation done\n')
end

