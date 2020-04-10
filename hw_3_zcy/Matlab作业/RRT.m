%***************************************
%Author: Chaoqun Wang
%Date: 2019-10-15
%***************************************
%% 流程初始化
clear all; close all;
x_I=1; y_I=1;           % 设置初始点
x_G=700; y_G=700;       % 设置目标点
Thr=50;                 %设置目标点阈值
Delta= 30;              % 设置扩展步长
FLAG_REWIRE = 1;        % RRT* 剪枝开关
rewire_radius = 2*Delta;% RRT*
FLAG_INFORMED = 1;      % Informed RRT* (generate random new point in a Ellipse that its const_dist = current path dist)
flag_path_found = 0;    % informed RRT* can only perform when a path is found
goal_index = 0;         % informed RRT* has to know last point to goal
flag_informed_done = 0; % in case we are not generating random point inside ellipse
max_iter_in = 2000;     % max iter number for inform random
prob_toward_goal = 0.8; % if p > prob_toward_goal select goal as new point
% Explanation:
% the const dist of ellipse = PA+PB, P is a point on curve, A and B are its
% two two focal points
%% 建树初始化
T.v(1).x = x_I;         % T是我们要做的树，v是节点，这里先把起始点加入到T里面来
T.v(1).y = y_I; 
T.v(1).xPrev = x_I;     % 起始节点的父节点仍然是其本身
T.v(1).yPrev = y_I;
T.v(1).dist=0;          %从父节点到该节点的距离，这里可取欧氏距离
T.v(1).indPrev = 0;     %
%% 开始构建树——作业部分
figure(1);
ImpRgb=imread('newmap.png');
Imp=rgb2gray(ImpRgb);
imshow(Imp)
xL=size(Imp,1);%地图x轴长度
yL=size(Imp,2);%地图y轴长度
hold on
plot(x_I, y_I, 'ro', 'MarkerSize',10, 'MarkerFaceColor','r');
plot(x_G, y_G, 'go', 'MarkerSize',10, 'MarkerFaceColor','g');% 绘制起点和目标点
count=1;
record_path_length = 0;
for iter = 1:3000
    x_rand=[];
    %Step 1: 在地图中随机采样一个点x_rand
    %提示：用（x_rand(1),x_rand(2)）表示环境中采样点的坐标
    x_rand = zeros(1,2);
    x_rand(1) = floor(xL * rand()) + 1;
    x_rand(2) = floor(yL * rand()) + 1;
    % Informed RRT* random sample in ellipse
    p_rand = rand();
    if p_rand > prob_toward_goal
        x_rand = [x_G y_G];
    elseif FLAG_INFORMED == 1 && flag_path_found == 1
        iter_in = 0;
        while 1
            x_rand = zeros(1,2);
            x_rand(1) = floor(xL * rand()) + 1;
            x_rand(2) = floor(yL * rand()) + 1;
            last_point = [T.v(goal_index).x T.v(goal_index).y];
            dist_temp = cal_dist(last_point,[x_G y_G]);
            const_dist = T.v(goal_index).dist + dist_temp;
            record_path_length = const_dist;
            a = const_dist/2; % half long axis
            c = cal_dist([x_I y_I],[x_G y_G])/2; % focal of ellipse
            b = sqrt(a^2-c^2); % half short axis
            theta = atan2(y_G-y_I,x_G-x_I);
            coord_local = [cos(theta) -sin(theta) x_I;sin(theta) cos(theta) y_I;0 0 1]^-1 * [x_rand.';1];
            coord_local(1) = coord_local(1) - c;
            if coord_local(1)^2/a^2 + coord_local(2)^2/b^2 <= 1
                break
            end
            if iter_in > max_iter_in
                flag_informed_done = 1;
%                 break;
            end
            iter_in = iter_in+1;
        end
    end
    x_near=[];
    %Step 2: 遍历树，从树中找到最近邻近点x_near 
    %提示：x_near已经在树T里
    x_near = [T.v(1).x T.v(1).y];
    min_dist = cal_dist(x_near,x_rand);
    x_near_idx = 1;
    for i = 1:size(T.v,2)
        x_near_temp = [T.v(i).x T.v(i).y];
        dist_temp = cal_dist(x_near_temp,x_rand);
        if dist_temp < min_dist
            x_near = x_near_temp;
            x_near_idx = i;
            min_dist = dist_temp;
        end
    end
    x_new=[];
    %Step 3: 扩展得到x_new节点
    %提示：注意使用扩展步长Delta
    theta = atan2(x_rand(2)-x_near(2),x_rand(1)-x_near(1));
    x_new = [Delta*cos(theta) Delta*sin(theta)] + x_near;
    %检查节点是否是collision-free
    if ~collisionChecking(x_near,x_new,Imp) 
       continue;
    end
    count=count+1;
    
    % Rewire for RRT*
    % cal path length from start to x_new
    c_path_dist = Delta + T.v(x_near_idx).dist;
    if FLAG_REWIRE == 1 || FLAG_INFORMED == 1
        % find neighbor nodes in rewire_radius
        for i = 1:size(T.v,2)
            if i == x_near_idx
                continue % skip x_near
            end
            x_near_temp = [T.v(i).x T.v(i).y];
            dist_temp = cal_dist(x_near_temp,x_new);
            if dist_temp < rewire_radius
                % rewire if path from x_near_i is shorter than its current
                % path from x_near
                n_path_dist = dist_temp + T.v(i).dist;
                if c_path_dist > n_path_dist
                    x_near_idx = i;
                    c_path_dist = n_path_dist;
                    x_near = [T.v(i).x T.v(i).y];
                end
            end
        end 
    end
    
    %Step 4: 将x_new插入树T 
    %提示：新节点x_new的父节点是x_near
    current_size = size(T.v,2);
    new_idx = current_size+1;
    T.v(new_idx).x = x_new(1);         % T是我们要做的树，v是节点
    T.v(new_idx).y = x_new(2); 
    T.v(new_idx).xPrev = x_near(1);     
    T.v(new_idx).yPrev = x_near(2);
    T.v(new_idx).dist = c_path_dist; %从父节点到该节点的距离，这里可取欧氏距离
    T.v(new_idx).indPrev = x_near_idx;     %
    if FLAG_INFORMED == 1 && flag_path_found == 1
        index_list = [goal_index];
        pathIndex = goal_index;
        while 1
            pathIndex = T.v(pathIndex).indPrev;
            index_list = [index_list;pathIndex];
            if pathIndex == 1
                break
            end
        end
        for i = length(index_list)-1:-1:1
            idx = index_list(i);
            idx_p = index_list(i+1);
            dist_temp = cal_dist([T.v(idx).x T.v(idx).y],[T.v(idx_p).x T.v(idx_p).y]);
            T.v(idx).dist = T.v(idx_p).dist + dist_temp;
        end
    end
    %Step 5:检查是否到达目标点附近 
    %提示：注意使用目标点阈值Thr，若当前节点和终点的欧式距离小于Thr，则跳出当前for循环
    fprintf(1,'path_length: %.3f\n',record_path_length);
    goal_dist = cal_dist(x_new, [x_G y_G]);
    if goal_dist < Thr 
        if FLAG_INFORMED == 0
            break
        else
            flag_path_found = 1;
            goal_index = new_idx;
        end
    end
   %Step 6:将x_near和x_new之间的路径画出来
   %提示 1：使用plot绘制，因为要多次在同一张图上绘制线段，所以每次使用plot后需要接上hold on命令
   %提示 2：在判断终点条件弹出for循环前，记得把x_near和x_new之间的路径画出来
   line_seg = [x_near;x_new];
   plot(line_seg(:,1),line_seg(:,2),'r.-')
   hold on
   pause(0.005); %暂停0.1s，使得RRT扩展过程容易观察
end
%% 路径已经找到，反向查询
if iter < 2000 || flag_path_found == 1
    path.pos(1).x = x_G; path.pos(1).y = y_G;
    path.pos(2).x = T.v(end).x; path.pos(2).y = T.v(end).y;
    pathIndex = T.v(end).indPrev; % 终点加入路径
    j=0;
    while 1
        path.pos(j+3).x = T.v(pathIndex).x;
        path.pos(j+3).y = T.v(pathIndex).y;
        pathIndex = T.v(pathIndex).indPrev;
        if pathIndex == 1
            break
        end
        j=j+1;
    end  % 沿终点回溯到起点
    path.pos(end+1).x = x_I; path.pos(end).y = y_I; % 起点加入路径
    for j = 2:length(path.pos)
        plot([path.pos(j).x; path.pos(j-1).x;], [path.pos(j).y; path.pos(j-1).y], 'b', 'Linewidth', 3);
    end
else
    disp('Error, no path found!');
end


