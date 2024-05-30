function [sensors] = DisplaySensor(fin, N, CR, L, W)
% This function display sensors with circles

sensors = zeros(N,6); % x1, y1, data-generation rate(DGR), energy consumption rate(ECR)

%figure;
rectangle('Position', [0, 0, L, W] );


i=0; % sensor index
while (i< N)
 
% if (round ==13 && i==10)
%     disp(round);
% end   
 
    i=i+1;
    [sensors(i,1:5), count] = fscanf(fin,'%f %f %f ',[1,5]); % Read sensors coordinates(xi,yi, data-generated, ECR))from file    
    pos = [sensors(i,1)-3 sensors(i,2)-3 6 6];  
    rectangle('Position', pos,'Curvature',[1 1], 'FaceColor',[0 .5 .5]); % Display sensor
    text(sensors(i,1)+5, sensors(i,2)-25, num2str(i));  % Display coordinates   
    DisplayCircle(sensors(i,1), sensors(i,2),CR);      % Display Communication region of the sub-sink 
end
       

    
end