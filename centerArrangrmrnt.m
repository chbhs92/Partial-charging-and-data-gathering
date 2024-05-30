function [bs1,sensor_x,sensor_y,circle_nodes3,dist]=centerArrangrmrnt(bs,dmat,x_co,y_co,radius)

% calculate  nodes inside circle or hexagonal cell
count=-1;
for i=1:length(bs)
   for j=1:length(y_co) 
     if( dmat(i,j)<= radius)
           index(i,j)=j;
           count=count+1;
           no_of_point(i,1)=count;
     end
   end
   count=-1;
end
% sort the points based on neighbour 
[value,Index]= sort(no_of_point,'descend');         %sort(point);
[n,m]=size(index);

s=1;
 
%circle nodes 
for i=1:n
    for j=1:m
        if(index(i,j)~=0)
           circle_point(i,s)=j;
           s=s+1;
        end
       
        
    end
    s=1;
end
% remove duplicate nodes
 
 [a,c] = unique(circle_point,'first');
 out = zeros(size(circle_point));
 out(c) = a;
  [rt,ct]=size(out);
 for i=1:rt
    for j=1:ct
        
          if( out(i,j)~=0)
             circle_nodes1(i,s)=out(i,j);
             s=s+1;
          end
 
        
    end
    s=1;
 end
circle_nodes2=circle_nodes1;
 [r1,c1]=size(circle_nodes1);
 s=1;
for i=1:r1
     if(  circle_nodes1(i,1)~=0)
         for j=1:c1
            circle_nodes3(s,j)=circle_nodes1(i,j);
            index_circle(s)=i; % circle_nodes1(i,1);
            
        end
        s=s+1;
     end
 
 
end
 % coordinate of centers and nodes
  [r2,c2]=size(circle_nodes3);
  for i=1:r2
    
          center_coordinate(i,1)=bs(index_circle(i),1);
          center_coordinate(i,2)=bs(index_circle(i),2);
          for j=1:c2
              if(circle_nodes3(i,j)~=0)
                nodes_x_coordinate(i,j)=x_co(1,circle_nodes3(i,j));
                nodes_y_coordinate(i,j)=y_co(1,circle_nodes3(i,j)); 
              end
          end
  end
 
%  plot(center_coordinate(:,1),center_coordinate(:,2),'g*')
bs1=center_coordinate;
sensor_x=nodes_x_coordinate;
sensor_y=nodes_y_coordinate;
 
  % calculate distace between sensor nodes and cenetr of circle
[r3,c3]=size(nodes_y_coordinate);
for i=1:length(bs1)
    for j=1:c3
        if(sensor_x(i,j)~=0 && sensor_y(i,j)~=0)
        dist(i,j)=sqrt((bs1(i,1)-sensor_x(i,j))^2+(bs1(i,2)-sensor_y(i,j))^2);
        end
    end
end

 a=10;


end