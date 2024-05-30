function [] = GenerateSensors()
% Generate N random sensors withn L*W area and store them in fin

N=200*50;

L=100;
W=100;

Max_Packets=10; % Data packets generation rate in a sensor randomly in between [1, 10] \cite{Ma:2018,Lyu:2020} ; [1,50] \cite{Xu:2020} 
Packet_size=1; % Packet size in kilo bits 

Min_ecr=0.05; % Joule/Sec
Max_ecr=0.5; % Joule/Sec \cite{Liu:2020} 
fin=fopen('input_sensors_100m_100m_area.txt','a');
% fin=fopen
i=1;
   while (i<= N)
%     dgr = randi([0 Max_Packets],1,1)*Packet_size;  % Data generation rate of a sensor
      ec = (0.1 + (1-0.1).*rand(1,1));     % Energy consumption rate of a sensor in between [Min_ecr to Max_ecr]  \cite{Z. Wei et al.}     
%     x =randi([0 L],1,1);  
%     y = randi([0 W],1,1);
      Rem_Energy = (1080+ (10800-1080).*rand(1,1));
      dgr=randi([1 10],1,1); %linspace(1,10,190);%data generation rate kbps
      x_coor=0+rand(1,1).*100 ;
      y_coor=0+rand(1,1).*100 ;
        % x,y, dgr, ecr
      i=i+1;
      fprintf(fin,'%7d  %7d   %7d   %7d  %7d\n', x_coor, y_coor, dgr, ec, Rem_Energy);
   end
fclose(fin);
 
end