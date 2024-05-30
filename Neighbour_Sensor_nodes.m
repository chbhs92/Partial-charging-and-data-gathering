function [Neighbour_node] =  Neighbour_Sensor_nodes(bs,sensor_x,sensor_y,circle_nodes,RC, CR)

for i=1:size(bs,1)
    for j = 1:size(sensor_x,1)
        d(i,j) = sqrt((bs(i,1) - sensor_x(j))^2 + (bs(i,2) - sensor_y(j))^2 );  % distance from base station to all sensor nodes 
    end
end
% condition for being neighbour
s=1;
for i =1 : size(d,1)
    for j = 1: size (d,2)
        if (d(i,j) > CR && d (i,j) <= RC)
            Neighbour_node(i,s) = j;
            s=s+1;
        end
    end
    s=1;
end


end