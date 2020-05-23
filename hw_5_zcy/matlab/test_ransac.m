path_points = [X_n;Y_n].';
x_size = size(path_points,1);
indices = linspace(0,x_size-1,x_size).';
% first fit x = f(t), t = indices
poly_x = indices;
poly_y = path_points(:,1);
[order_x, coeff_x] = my_ransac(poly_x,poly_y,'x',0);
% then fit y = f(t), t = indices
poly_y = path_points(:,2);
[order_y, coeff_y] = my_ransac(poly_x,poly_y,'y',0);
% new points & original comparison
p_x_new = polyval(coeff_x,indices);
p_y_new = polyval(coeff_y,indices);
% plot
figure
plot(path_points(:,1),path_points(:,2),'.-')
hold on
plot(p_x_new,p_y_new,'.-')
legend('original','new')
title('new points & original comparison')
