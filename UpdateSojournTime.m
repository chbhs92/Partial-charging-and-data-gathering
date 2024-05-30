function  [sojourn_time] = UpdateSojournTime(T,DGR,Emax,U_sensor, circle_nodes)
% charging time and data collection time
G=500 ;% data collection rate of MD in kbps
[r1,c] = size(circle_nodes);
for i=1:r1
    for j=1:c
        if(circle_nodes(i,j)~=0)
%             charging_time(i,j)=(T*ECR(i,j))/U_sensor(i,j);
             charging_time(i,j)=Emax/U_sensor(i,j);
           % charging_time(i,j)=randi([ceil(Emin/U_sensor(i,j)),ceil(max(0,Emax/U_sensor(i,j)))],1,1); %(T*ECR(i,j))/U_sensor(i,j);
            %data_collection_time(i,j)=(T*DGR(i,j)/G);
            data_collection_time(i,j)=(T*DGR(i,j))/G;
            Data_generated(i,j) = T*DGR(i,j);  % data generated and transmitted to MV.
            % data_collection_time1(i,j)=(charging_time1(i,j)*DGR(i,j))/G;
            
        end
    end
end
% charging max and data collection max in each cell
for i=1:r1
    charge_max(i,1)=max(charging_time(i,:));
    data_collection_max(i,1)=max(data_collection_time(i,:));
end
% sojourn time in each cell

for i=1:r1
    
    sojourn_time(i,1)=max(charge_max(i,1),data_collection_max(i,1));
    
end