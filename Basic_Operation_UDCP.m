function [ sojourn_time,working_time,Tmax,DGR,U_sensor,ECR,Energy,charging_time,data_collection_time,Data_generated, RemStorage]=Basic_Operation_UDCP(bs,T,G,radius,Emax,Emin,circle_nodes,p,r,dist,U_full,SCsensor)
circle_nodes3=circle_nodes;
% data collection time
n=length(r);
[~,m]=size(circle_nodes);
% for i=1:n*m
%
%     DataCollection_time1(i)=T*r(i)/G;
%
% end
% coordinate of centers and nodes
[r1,c]=size(circle_nodes3);
for i=1:r1
    for j=1:c
        if(circle_nodes3(i,j)~=0)
            ECR(i,j)=p(1,circle_nodes3(i,j));
            DGR(i,j)=r(1,circle_nodes3(i,j));
            
        end
    end
end
% charging recieving rate
%%% equation for charging receiving x =0%(0 + (2.7).*rand(1,5));
% U_sensor=(-0.0985*dist^2 - 0.0377*dist +1)*5
for i=1:r1
    for j=1:c
        if(circle_nodes3(i,j)~=0)
            if(dist(i,j)>radius)
                U_sensor(i,j)=0;
            end
        if (dist(i,j)>0 && dist(i,j)< radius)
               U_sensor(i,j)=(-0.0985*dist(i,j)^2 - 0.0377*dist(i,j) + 1)*U_full;
               %U_sensor(i,j)=(4.32*10^-4/(dist(i,j) + 0.2316)^2)*U_full;
        end
            %             if(dist(i,j)>0 && dist(i,j)<=0.18*radius )  % energy efficiency 4%
%                 U_sensor(i,j)=4.7921 + (5.0-4.7921)*rand(1);
%             end
    
%             if(dist(i,j)>0.18*radius && dist(i,j)<=0.44*radius )   % energy efficiency 6.25%
%                 U_sensor(i,j)=4.081 + (4.7921-4.081)*rand(1);
%             end
%             if(dist(i,j)>0.44*radius && dist(i,j)<=0.74*radius )   % energy efficiency 14%
%                 U_sensor(i,j)= 2.6573 + (4.081- 2.6573)*rand(1);
%             end
%             if(dist(i,j)>0.74*radius && dist(i,j)<=0.85*radius )     % energy efficiency 20%
%                 U_sensor(i,j)= 1.9734 + (2.6573-1.9734)*rand(1);
%             end
%             if(dist(i,j)>0.85*radius && dist(i,j)<=radius )    % energy efficiency 24%
%                 U_sensor(i,j)=   0.9007 + (1.9734 - 0.9007)*rand(1);
%             end
            if (dist(i,j)==0)
                U_sensor(i,j) = U_full;
            end
        end
    end
end

% charging time and data collection time

for i=1:r1
    for j=1:c
        if(circle_nodes3(i,j)~=0)
           charging_time(i,j)=(T*ECR(i,j))/U_sensor(i,j);
             %charging_time(i,j)=Emax/U_sensor(i,j);
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
% total sojourn time
working_time=sum(sojourn_time);

% optimised time period
for i=1:r1
    for j=1:c
        if(circle_nodes3(i,j)~=0)
            Kp(i,j)=ECR(i,j)/U_sensor(i,j);
            Kr(i,j)=DGR(i,j)/G;
        end
    end
end
for i=1:r1
    kp_max(i,1)=max(Kp(i,:));
    kr_max(i,1)=max(Kr(i,:));
    Kj(i,1)=max(kp_max(i,1),kr_max(i,1));
end

for i=1:r1
    for j=1:c
        if(circle_nodes3(i,j)~=0)
            T_cell(i,1)=max((Emax-Emin)/((1-Kj(i,1))*ECR(i,j)));
        end
    end
end
for i=1:r1
    Kj(i,1)=Kj(i,1)*T;
end
working_time2=sum(Kj);
Tmax=min(T_cell);
%fprintf("optimised time period =%f\n",Tmax);
for i=1:r1
    for j=1:c
        if(circle_nodes3(i,j)~=0)
            Energy(i,j)=Emax;
            RemStorage(i,j) = SCsensor;
        end
    end
end
end

