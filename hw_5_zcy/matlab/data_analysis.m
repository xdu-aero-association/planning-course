function [mse,rmse,var_y,r_squared] = data_analysis(poly_x, poly_y, n, coeff_in)
% Calculate mse,rmse,var_y,r_squared comparing with self(observer)
% mse: Mean Squared Error
% rmse: Root Mean Squared Error
% var_y: Variation of y
% r_squared: 
% straight line polyfit & rmse
if nargin  < 4
    coeff = polyfit(poly_x, poly_y, n);
else
    coeff = coeff_in;
end
sq_sum = 0;
var_y = 0;
y_mean = mean(poly_y);
num_data = size(poly_x, 1);
for i = 1:num_data
    y_ideal = cal_y_value(coeff, poly_x(i), n);
    sq_sum = sq_sum + (poly_y(i) - y_ideal)^2;
    var_y = var_y + (poly_y(i) - y_mean)^2;
end
mse = sq_sum / num_data;
rmse = sqrt(mse);
var_y = var_y / num_data;
r_squared = 1 - mse / var_y;
end

function [y] = cal_y_value(coeff, x, n)
y = coeff(1);
for i = 2:n+1
    y = y * x + coeff(i);
end
end