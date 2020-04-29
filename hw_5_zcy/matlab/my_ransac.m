function [order, coeff_out] = my_ransac(poly_x,poly_y,string,ArB_exist_flag)
% Use RANSAC to get curve fit order
% Params
thres_dist = 0.3; % (m)
max_iter = 25;
fitted_rate = 0;
ideal_fitted_rate = 0.50;
init_fitted_rate = 0.5;
order_list = [1,2,3,4,5,6,7,8,9]; % too high order leads to unstable numeric results!!!!!

% transform poly points in to its start point coord frame
% fprintf(1,'ArB exist status:%d\n',ArB_exist_flag)
% [poly_x, poly_y] = Transform2origin(poly_x, poly_y,ArB_exist_flag);
% if ArB_exist_flag ~= 1
%     assignin('base','x_start',poly_x(1))
%     assignin('base','x_end',poly_x(end))
% end
% resample points
% [poly_x_ransac,poly_y_ransac] = resample_points(poly_x,poly_y);
% poly_x = poly_x_ransac;
% poly_y = poly_y_ransac;

max_num = length(poly_x);
assert((max_num >= 2), 'Not enough points!!!');
rand_sel_num = max(ceil(max_num * init_fitted_rate), 2);
ideal_num = max(ceil(max_num * ideal_fitted_rate), 2);
best_fitted_rate = zeros(length(order_list),1);
list_len_set = zeros(length(order_list),1);
coeff_set = zeros(length(order_list), order_list(end)+1);

for j = 1:length(order_list)
    best_r_squared = 0;
    best_mse = 0;
    first_flag = 0;
    n = order_list(j);
    list_len = 0;
    iter = 0;
    best_coeff = zeros(1, n+1);
    while iter < max_iter
        idx_list = randperm(max_num, rand_sel_num).';
        idx_list_remain = get_remain_idx(idx_list, max_num);
        coeff_temp = polyfit(poly_x(idx_list),poly_y(idx_list),n);
        for i = 1:length(idx_list_remain)
            x_data = poly_x(idx_list_remain(i));
            y_data = poly_y(idx_list_remain(i));
            y_model = cal_y_value(coeff_temp, x_data, n);
            dist = cal_dist(x_data, y_data, x_data, y_model);
            if dist < thres_dist
                idx_list = [idx_list;idx_list_remain(i)];
            end
        end
        list_len = length(idx_list);
        if length(idx_list) > ideal_num
            % this implies that we may have found a good model
            % now test how good it is
            better_coeff = polyfit(poly_x(idx_list),poly_y(idx_list),n);
%             assignin('base','better_coeff',better_coeff);
            [mse,rmse,var_y,r_squared] = data_analysis(poly_x(idx_list),poly_y(idx_list), n);
%             assignin('base','r_squared',r_squared);
%             if best_r_squared < r_squared
%                 best_r_squared = r_squared;
%                 best_coeff = better_coeff;
%             end
            if first_flag == 0
                first_flag = 1;
                best_mse = mse;
                best_coeff = better_coeff;
            elseif best_mse > mse
                best_mse = mse;
                best_coeff = better_coeff;
            end
        end
        iter = iter + 1;
    end
    list_len_set(j) = list_len;
%     best_fitted_rate(j) = best_r_squared;
    best_fitted_rate(j) = best_mse;
    coeff_set(j,1:n+1) = best_coeff;
end
assignin('base','best_fitted_rate',best_fitted_rate); 
assignin('base','list_len_set',list_len_set); 
assignin('base','coeff_set',coeff_set); 
% [v,max_idx] = max(best_fitted_rate);
[v,max_idx] = min(best_fitted_rate);
[v,max_num_idx] = max(list_len_set);
% max_idx = min(max_idx, max_num_idx);
order = order_list(max_idx);
% set order directly
% order = 3;
coeff_out = coeff_set(order, 1:order+1);
% if nargout == 0
figure
plot(poly_x, poly_y, '.-')
y_fitted = zeros(length(poly_x),1);
for i = 1:length(poly_x)
    y_fitted(i) = cal_y_value(coeff_out, poly_x(i), order);
end
hold on
plot(poly_x, y_fitted, '.-')
axis equal
title(['ransac result comparison, order: ', num2str(order)])
hold on 
plot(poly_x(1),poly_y(1),'*')
legend(['original ' string], 'fitted','start')

% end

end

function idx_list_remain = get_remain_idx(idx_list, len)
    idx_list = insert_sort(idx_list);
%     assignin('base','idx_list',idx_list);
    idx_list_remain = zeros(len - length(idx_list), 1);
    i = 1;
    j = 1;
    idx = 1;
    while idx <= len && i <= length(idx_list) && j <= length(idx_list_remain)
        if idx_list(i) ~= idx
            idx_list_remain(j) = idx;
            j = j + 1;
            idx = idx + 1;
        else
            i = i + 1;
            idx = idx + 1;
        end
    end
    while idx <= len && j <= length(idx_list_remain)
        idx_list_remain(j) = idx;
        j = j + 1;
        idx = idx + 1;
    end
    assignin('base','idx_list_remain',idx_list_remain);
end

function [y] = cal_y_value(coeff, x, n)
y = coeff(1);
% assignin('base','coeff',coeff)
% assignin('base','n_temp',coeff)
for i = 2:n+1
    y = y * x + coeff(i);
end
end

function [x] = insert_sort(x)
len = length(x);
for j = 2:len
    i = j - 1;
    key = x(j);
    while i > 0 && x(i) > key
        x(i + 1) = x(i);
        i = i - 1;
    end
    x(i + 1) = key;
end
end

function [dist] = cal_dist(x1, y1, x2, y2)
    dist = sqrt((x1-x2)^2 + (y1-y2)^2);
end

function [x,y] = Transform2origin(x,y,ArB_exist_flag)
    theta = atan((y(end)-y(1)) / (x(end)-x(1)));
    ArB = [cos(theta) -sin(theta);sin(theta) cos(theta)];
    AtB = [x(1);y(1)];
    for i = 1:length(x)
        point_new = ArB^-1 * ([x(i);y(i)] - AtB);
        x(i) = point_new(1);
        y(i) = point_new(2);
    end
    if ArB_exist_flag ~= 1
        assignin('base','ArB',ArB);
        assignin('base','AtB',AtB);
    end
end
