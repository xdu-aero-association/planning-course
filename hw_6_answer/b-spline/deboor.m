function p = deboor(degree, c_pts, knot, t, i)
% deboor algorithm
% degree: degree of b-spline
% c_pts: control points
% knot: array of knot positions
% t: position
% i: index of knot interval that contains t

% pia = zeros(20,2);
n = size(c_pts,1)-1;
pia = zeros(n+2,2); % auxiliary array

for j = i-degree:i
    pia(j,:) = c_pts(j,:);
end

% fprintf(1,'n: %d, i:, %d, t:%.1f\n',n,i,t)

t1_set = [];
t2_set = [];
for k = 1:degree
%     for j = i+1:-1:i-degree+1 % working & applied in scu
    for j = i:-1:i-degree+k % working
%     for j = i-degree+1:i+1
%         i,j
        if knot(j+degree-k+1) ~= knot(j)
            t1 = (knot(j+degree-k+1)-t)/(knot(j+degree-k+1)-knot(j));
        else
            t1 = 0;
        end
        t2 = 1-t1;
        pia(j,:) = t1*pia(j-1,:) + t2*pia(j,:);
        t1_set = [t1_set;t1];
        t2_set = [t2_set;t2];
    end
end
p = pia(i,:);

assignin('base','t1_set',t1_set);
assignin('base','t2_set',t2_set);

% hold on
% plot(p(:,1),p(:,2),'.');
% axis equal
% pause();

end

