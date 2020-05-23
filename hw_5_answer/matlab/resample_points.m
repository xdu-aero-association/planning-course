function [x,y] = resample_points(x_in,y_in)
% resample points at fixed interval
current_idx = 1;
min_dist = 0.1;
num_data = length(x_in);
assert((num_data >= 2), 'Not enough points!!!');
resampled_flags = zeros(num_data,1);
resampled_flags(current_idx) = 1;
next_idx = current_idx + 1;
resampled_num = 1;
while next_idx < num_data
    if cal_dist(x_in(current_idx),y_in(current_idx),x_in(next_idx),y_in(next_idx)) < min_dist
        next_idx = next_idx + 1;
    else
        current_idx = next_idx;
        next_idx = next_idx + 1;
        resampled_flags(current_idx) = 1;
        resampled_num = resampled_num + 1;
    end
end
resampled_flags(next_idx) = 1;
resampled_num = resampled_num + 1;
x = zeros(resampled_num,1);
y = x;
j = 1;
for i = 1:num_data
    if resampled_flags(i) == 1
        x(j) = x_in(i);
        y(j) = y_in(i);
        j = j + 1;
    end
end

end

function [dist] = cal_dist(x1, y1, x2, y2)
    dist = sqrt((x1-x2)^2 + (y1-y2)^2);
end