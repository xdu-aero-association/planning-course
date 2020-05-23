function path = A_star_search(map,MAX_X,MAX_Y)
%%
%This part is about map/obstacle/and other settings
%pre-process the grid map, add offset
size_map = size(map,1);
Y_offset = 0;
X_offset = 0;

%Define the 2D grid map array.
%Obstacle=-1, Target = 0, Start=1
MAP=2*(ones(MAX_X,MAX_Y));%  2 means normal  not obs

%Initialize MAP with location of the target
xval=floor(map(size_map, 1)) + X_offset;
yval=floor(map(size_map, 2)) + Y_offset;
xTarget=xval;
yTarget=yval;
MAP(xval,yval)=0;% means target

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
% means start point

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%LISTS USED FOR ALGORITHM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%OPEN LIST STRUCTURE
%--------------------------------------------------------------------------
%| 1X val |2Y val |3Parent X val |4Parent Y val |5h(n) |6g(n)|7f(n)|
%--------------------------------------------------------------------------
OPEN=[];
%CLOSED LIST STRUCTURE
%--------------
%X val | Y val |  Parent X val |  Parent y val  |
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
            CLOSED(k,3)=0;
            CLOSED(k,4)=0;
            k=k+1;
        end
    end
end
CLOSED_COUNT=size(CLOSED,1);
%set the starting node as the first node
xNode=xval;
yNode=yval;
OPEN_COUNT=1;
goal_distance=distance(xNode,yNode,xTarget,yTarget);
path_cost=0;
OPEN(OPEN_COUNT,:)=[xNode,yNode,xNode,yNode,goal_distance,path_cost,goal_distance];
% OPEN(OPEN_COUNT,1)=0;
CLOSED_COUNT=CLOSED_COUNT+1;
CLOSED(CLOSED_COUNT,1)=xNode;
CLOSED(CLOSED_COUNT,2)=yNode;
CLOSED(CLOSED_COUNT,3)=0;
CLOSED(CLOSED_COUNT,4)=0;
NoPath=1;

%%
%This part is your homework
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START ALGORITHM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

not_finish = 1;

a_star_fail = 0;

target_index_inopen = 0;
while(not_finish) %you have to dicide the Conditions for while loop exit
    
    
    OPEN = sortrows(OPEN,7);
    now_point(1) = OPEN(1,1);
    now_point(2) = OPEN(1,2);
    
    for i_nei = now_point(1)-1:now_point(1)+1 % search the neighbor
        for j_nei = now_point(2)-1:now_point(2)+1
            
            if((i_nei== now_point(1))&&(j_nei== now_point(2))) % is not itself
            else
                if ((i_nei<1)||(i_nei>MAX_X)||(j_nei<1)||(j_nei>MAX_Y)) % is not outof the range
                else
                    
                    in_close =0;
                    for i_clo =1: size(CLOSED,1)
                        
                        if((i_nei== CLOSED(i_clo,1))&&(j_nei== CLOSED(i_clo,2)))%
                            in_close =1;
                            
                            break;
                            
                        end
                        
                    end
                    
                    if(~(in_close))% not in close list
                        
                        in_open =0;
                        index_in_open=0;
                        for i_ope =1: size(OPEN,1)
                            
                            if ((OPEN(i_ope,1)==i_nei)&&(OPEN(i_ope,2)==j_nei))%in open list
                                in_open =1;
                                index_in_open = i_ope;
                                break;
                                
                            end
                            
                        end
                        
                        if in_open
                            
                            g = OPEN(1,6)+ distance(i_nei,j_nei,now_point(1),now_point(2));
                            
                            if OPEN(index_in_open,6)>g %change its parent
                                OPEN(index_in_open,6) =g;
                                OPEN(index_in_open,3) =now_point(1);
                                OPEN(index_in_open,4) =now_point(2);
                                OPEN(index_in_open,7) = OPEN(index_in_open,5) + OPEN(index_in_open,6);
                            end
                            
                            
                        else
                            g = OPEN(1,6)+ distance(i_nei,j_nei,now_point(1),now_point(2));
                            h =  distance(i_nei,j_nei,xTarget,yTarget);
                            OPEN = [OPEN;[i_nei,j_nei, now_point(1), now_point(2),h,g,h+g]];
                            
                            
                        end
                        
                    end
                    
                    
                end
            end
            
        end
    end
    
    
    CLOSED=[CLOSED;[now_point(1),now_point(2),OPEN(1,3),OPEN(1,4)]];
    OPEN(1,:)=[];% delete this point and add it to close list
    
    
    % if finish ?
    
    if isempty(OPEN)
        
        not_finish = 0;
        a_star_fail = 1;
        
    end
    for i =1: size(OPEN,1)
        
        if ((OPEN(i,1)==xTarget)&&(OPEN(i,2)==yTarget))
            a_star_fail = 0;
            not_finish = 0;
            target_index_inopen = i;
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
path = [];
if (~a_star_fail)
    
    path(1,1) = xTarget;
    path(1,2) = yTarget;
    
    find_start = 0;
    parent_x  = OPEN(target_index_inopen,3);
    parent_y  = OPEN(target_index_inopen,4);
    path = [path; [parent_x,parent_y]];
    while(~find_start)
        
        for i_clo =1: size(CLOSED,1)
            
            if((parent_x== CLOSED(i_clo,1))&&(parent_y== CLOSED(i_clo,2)))%
                
                parent_x  = CLOSED(i_clo,3);
                parent_y  = CLOSED(i_clo,4);
                path = [path; [parent_x,parent_y]];
                break;
                
            end
            
        end
        
        if ((parent_x==xStart )&&(parent_y==yStart))
            find_start = 1;
            
%              path = [path; [xStart,yStart]];
        end
        
    end
    
end



end
