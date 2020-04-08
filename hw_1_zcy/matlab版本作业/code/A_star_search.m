function path = A_star_search(map,MAX_X,MAX_Y)
%%
%This part is about map/obstacle/and other settings
    %pre-process the grid map, add offset
    size_map = size(map,1);
    Y_offset = 0;
    X_offset = 0;
    
    %Define the 2D grid map array.
    %Obstacle=-1, Target = 0, Start=1
    MAP=2*(ones(MAX_X,MAX_Y));
    
    %Initialize MAP with location of the target
    xval=floor(map(size_map, 1)) + X_offset;
    yval=floor(map(size_map, 2)) + Y_offset;
    xTarget=xval;
    yTarget=yval;
    MAP(xval,yval)=0;
    
    %Initialize MAP with location of the obstacle
    for i = 2: size_map-1
        xval=floor(map(i, 1)) + X_offset;
        yval=floor(map(i, 2)) + Y_offset;
        MAP(xval,yval)=-1;
    end 
    
    %Initialize MAP with location of the start point
    xval=floor(map(1, 1)) + X_offset;
    yval=floor(map(1, 2)) + Y_offset;
    xStart=xval;
    yStart=yval;
    MAP(xval,yval)=1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %LISTS USED FOR ALGORITHM
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %OPEN LIST STRUCTURE
    %--------------------------------------------------------------------------
    %IS ON LIST 1/0 |X val |Y val |Parent X val |Parent Y val |h(n) |g(n)|f(n)|
    %--------------------------------------------------------------------------
    OPEN=[];
    %CLOSED LIST STRUCTURE
    %--------------
    %X val | Y val |
    %--------------
    % CLOSED=zeros(MAX_VAL,2);
    CLOSED=[];

    %Put all obstacles on the Closed list
    k=1;%Dummy counter
    for i=1:MAX_X
        for j=1:MAX_Y
            if(MAP(i,j) == -1)
                CLOSED(k,1)=i;
                CLOSED(k,2)=j;
                k=k+1;
            end
        end
    end
    CLOSED_COUNT=size(CLOSED,1);
    %set the starting node as the first node
    xNode=xStart;
    yNode=yStart;
    OPEN_COUNT=1;
    goal_distance=distance(xNode,yNode,xTarget,yTarget);
    path_cost=0;
    OPEN(OPEN_COUNT,:)=insert_open(xNode,yNode,xNode,yNode,goal_distance,path_cost,goal_distance);
%     OPEN(OPEN_COUNT,1)=0;
%     CLOSED_COUNT=CLOSED_COUNT+1;
%     CLOSED(CLOSED_COUNT,1)=xNode;
%     CLOSED(CLOSED_COUNT,2)=yNode;
    NoPath=1;
    current_x = xStart;
    current_y = yStart;
%%
%This part is your homework
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START ALGORITHM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    while(1) %you have to dicide the Conditions for while loop exit 
        % check if open set is empty or not
        current_idx = min_fn(OPEN,OPEN_COUNT,xTarget,yTarget);
        if current_idx == -1
            fprintf(1,'No Path Found !!!\n')
            break
        end
        current_x = OPEN(current_idx,2);
        current_y = OPEN(current_idx,3);
        OPEN(current_idx,1) = 0;
        % check if we have reached goal
        if xTarget == current_x && yTarget == current_y
            fprintf(1,'Goal Reached !!!\n')
            NoPath = 0;
            break;
        end
        % find node with lowest fn & set as current node & pop it out and
        % add it to closed set
        CLOSED_COUNT=CLOSED_COUNT+1;
        CLOSED(CLOSED_COUNT,1)=current_x;
        CLOSED(CLOSED_COUNT,2)=current_y;
        % find all valid neighbors
        current_gn = OPEN(current_idx,7);
        exp_array = expand_array(current_x,current_y,current_gn,xTarget,yTarget,CLOSED,MAX_X,MAX_Y);
        for i = 1:size(exp_array,1)
            neighbor_x = exp_array(i,1);
            neighbor_y = exp_array(i,2);
            neighbor_idx = node_index(OPEN,neighbor_x,neighbor_y);
            % if neighbor in closed set & new gn is smaller
            if neighbor_idx <= size(OPEN,1) % && OPEN(neighbor_idx,1) == 0
                if OPEN(neighbor_idx,7) >= exp_array(i,4)
                    OPEN(neighbor_idx,6) = exp_array(i,3);
                    OPEN(neighbor_idx,7) = exp_array(i,4);
                    OPEN(neighbor_idx,8) = exp_array(i,5);
                    OPEN(neighbor_idx,4) = current_x;
                    OPEN(neighbor_idx,5) = current_y;
                end
%             elseif neighbor_idx <= size(OPEN,1) && OPEN(neighbor_idx,1) == 1
%                 if OPEN(neighbor_idx,7) >= exp_array(i,4)
%                     OPEN(neighbor_idx,6) = exp_array(i,3);
%                     OPEN(neighbor_idx,7) = exp_array(i,4);
%                     OPEN(neighbor_idx,8) = exp_array(i,5);
%                     OPEN(neighbor_idx,4) = current_x;
%                     OPEN(neighbor_idx,5) = current_y;
%                 end
            % if it's a new node
            else
                OPEN_COUNT = OPEN_COUNT+1;
                OPEN(OPEN_COUNT,:)=insert_open(neighbor_x,neighbor_y,current_x,current_y,exp_array(i,3),exp_array(i,4),exp_array(i,5));
            end
        end
    end %End of While Loop
    
    %Once algorithm has run The optimal path is generated by starting of at the
    %last node(if it is the target node) and then identifying its parent node
    %until it reaches the start node.This is the optimal path
    
    %
    %How to get the optimal path after A_star search?
    %please finish it
    %
    assignin('base','OPEN',OPEN);
    path = [];
    if NoPath == 0
        path = [current_x current_y;path];
        while current_x ~= xStart || current_y ~= yStart
            node_idx = node_index(OPEN,current_x,current_y);
            current_x = OPEN(node_idx,4);
            current_y = OPEN(node_idx,5);
            path = [current_x current_y;path];
        end
        path = [current_x current_y;path];
    end
end
