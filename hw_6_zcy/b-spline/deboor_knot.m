function [knot] = deboor_knot(degree,n,type)
% compute deboor knots
knot = zeros(n+degree+2,1);
if type == 0 % clamped at start
    for i = 0:n+degree+1
        if i <= degree
            knot(i+1) = 0;
        else 
            knot(i+1) = i-degree;
        end
    end
elseif type == 1 % no clamp
    for i = 0:n+degree+1
        knot(i+1) = i;
    end
else % clamp both ends
    for i = 0:n+degree+1
        if i <= degree
            knot(i+1) = 0;
        elseif i>degree && i <= n
            knot(i+1) = i-degree;
        else
            knot(i+1) = n-degree+1;
        end
    end
end
% for i = 0:n+degree+1
%     knot(i+1) = i;
% end
% knot = zeros(n+1,1);
% for i = 0:n
%     if i <= degree
%         knot(i+1) = 0;
%     elseif i>degree && i <= n
%         knot(i+1) = i-degree;
%     end
% end

% knot
% size(knot)


end

