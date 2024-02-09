clc
clear all;
% variable initialization
Emax=10800; % maximum energy of a node in joule
Emin=540 ; % minimum energy of node in joule
Etrav=20000; % travelling enrgy of MV 20KJ
Ech=20000; % travelling enrgy of MV 20KJ
BCsensor=Emax;
BCMV=200000; %  Battery capacity of MV in joule 400Kj
v=5; % velocity in m/sec
Pt=600 ; %J/m (Ma et al. 2018; Tomar et al. 2019; Liang et al.2016) ;% power of MD in travelling
G=500 ;% data collection rate of MD in kbps
% T=3600;
T=10800; % for 3 hours
U_full=10;  % maximum charging rate of MD in watt (joule/sec)[Z. Wei et al. 2020]
radius = 2.7;  % in  meter
SCsensor = 216*10^6; % in byte (216MB)
load input_sensors_100m_100m_area.txt
fin=fopen('input_sensors_1000m_1000m_area.txt','r');
CR=2.7;  % charging range ( in m )
RC = 40; % communication range (in m)
L=100;
W=100;
delta=Emax/5;    % mJ Amount of energy transfered is in the multiples of delta
ett1 = delta;
%ett=delta/err;
round=1;
BSX = 50;
BSY = 50;
ObjectiveFileIGWO = fopen('Objective_valueIGWO_125.txt','a');
ObjectiveFileDFA = fopen('Objective_value_DFA_125.txt','a');
NMV = fopen('Number_of MV_300.txt','a');

%%%%%%%%%%%%% IGWO %%%%%%%%%%%%%%%%%%%%%%%%%%% Average
TimeTaken = fopen('Average_Total_Time_taken.txt','a');
EnergyGain = fopen('Average_Energy_Gained.txt','a');
ChargeTime = fopen('Average_Charging_Time.txt','a');
DeadPeriod = fopen('Average_Dead_Period.txt','a');
TotalDistance = fopen('Average_Distance.txt','a');
EUeff = fopen('Average_Energy_Utilization.txt','a');
EnergyConsumed = fopen('Average_Energy_consumed.txt','a');
deadSensor = fopen('Average_Dead_sensor.txt','a');
TDeadPeriod = fopen('Average_Total_Dead_period.txt','a');


%%%%%%%%%%%%% IGWO %%%%%%%%%%%%%%%%%%%%%%%%%%% total
TimeTaken1 = fopen('Total_Total_Time_taken_100.txt','a');
EnergyGain1 = fopen('Total_Energy_Gained_100.txt','a');
ChargeTime1 = fopen('Total_Charging_Time_100.txt','a');
TotalDistance1 = fopen('Total_Distance_100.txt','a');
EUeff1 = fopen('Total_Energy_Utilization_100.txt','a');
EnergyConsumed1 = fopen('Total_Energy_consumed_100.txt','a');
deadSensor1 = fopen('Total_Dead_Sensor_100.txt','a');
DeadPeriod1 = fopen('Total_Dead_Period_100.txt','a');
% AnchorPoint = fopen('Number of anchor point_N300.txt','a');
% AverageAnchorPoint = fopen('Average_number of anchor point_N300.txt','a');
N=300; % no. of sensor nodes
%anchor = zeros(20,3);
while( round <=20)
    %% Read sensors input from file
    x_co1= input_sensors_100m_100m_area(((round-1)*N+1:round*N),1);
    y_co1= input_sensors_100m_100m_area(((round-1)*N+1:round*N),2) ;
    dgr=input_sensors_100m_100m_area(((round-1)*N+1:round*N),3);
    ecr=input_sensors_100m_100m_area(((round-1)*N+1:round*N),4);
    x_co=x_co1'; % x co-ordinate
    y_co=y_co1';  % y co-ordinate
    r=dgr';   % data generation rate
    p=ecr';    % energy consumption rate
    %fprintf(fout,'%2d ', round)
    %%%% covering the sensor nodes using minimum no. of circles Proposed Approaach
    [bs1,sensor_x1,sensor_y1,circle_nodes1,dist1] = UDCP_greedy_technique( x_co,y_co,radius);
    [bs,sensor_x,sensor_y,circle_nodes,dist] = IntersectionCircle(bs1,radius,circle_nodes1,x_co,y_co);
    [Neighbour_nodes] =  Neighbour_Sensor_nodes(bs,sensor_x1,sensor_y1,circle_nodes,RC, CR);  %%%% neighbour nodes
    [sojourn_time,working_time,Tmax,DGR,CHR,ECR, RemainingEnergy,charging_time, data_collection_time,Data_generated, RemStorage] = Basic_Operation_UDCP(bs,T,G,radius,Emax,Emin,circle_nodes,p,r,dist,U_full,SCsensor);
   ett = ett1./charging_time;
    s=1;
    for i=1:size(CHR,1)
        for j=1:size(CHR,2)
            if (CHR(i,j) ~= 0)
                ProposedEeff(s) = 1-((CHR(i,j)*U_full)/100); %% energy loss in proposed scheme
                s = s+1;
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%% hexagonal structure proposed by Lyu:2020  %%%%%%%%%%%%%%%%%% %%%%%%%%%%%
    [hex_bs,hex_sensor_x,hex_sensor_y,hex_circle_nodes,hex_dist] = cell_layout(radius,x_co,y_co);  % 100*100 m^2 area hexagonal cell
    [hex_sojourn_time,hex_working_time,hex_Tmax,hex_DGR,hex_CHR,hex_ECR,hex_Energy,hex_charging_time1,hex_data_collection_time,hex_Data_generated, hex_RemStorage] = Basic_Operation_UDCP(hex_bs,T,G,radius,Emax,Emin,hex_circle_nodes,p,r,hex_dist, U_full,SCsensor);
     hex_charging_time =  hex_charging_time1; %  charging time
    s1=1;
    for i=1:size(hex_CHR,1)
        for j=1:size(hex_CHR,2)
            if (hex_CHR(i,j) ~= 0)
                hexEeff(s1) = 1-((hex_CHR(i,j)*U_full)/100); %%% energy loss in hexagonal structure
                s1 = s1+1;
            end
        end
    end
    
    %    fprintf(AnchorPoint,' %d %d %d \n', length(bs), length(bs1), length(hex_bs));
    %    anchor(round,1) = length(bs); anchor(round,2) = length(bs1); anchor(round,3) = length(hex_bs);
    %     round = round+1
    % end
    %  fprintf(AverageAnchorPoint,' %d %d %d \n', ceil(mean(anchor(:,1))), ceil(mean(anchor(:,2))), ceil(mean(anchor(:,3))));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Anchor points proposed by Boukerche 2021 %%%%%%%%%%%%%%%%
    [bs11,sensor_x11,sensor_y11,circle_nodes11,dist11] = UDCP_greedy_technique( x_co,y_co,radius);
    [bs21,sensor_x21,sensor_y21,circle_nodes21,dist21] = WeightedAnchorPoints(bs11,radius,circle_nodes11,x_co,y_co);
    [sojourn_time21,working_time21,Tmax21,DGR21,CHR21,ECR21, RemainingEnergy21,charging_time212, data_collection_time21,Data_generated21, RemStorage21] = Basic_Operation_UDCP(bs11,T,G,radius,Emax,Emin,circle_nodes11,p,r,dist11,U_full,SCsensor);
     charging_time21 = 0.9 * charging_time212; % 90% charging time
    s2=1;
    for i=1:size(CHR21,1)
        for j=1:size(CHR21,2)
            if (CHR21(i,j) ~= 0)
                MPSSEeff(s2) = 1-((CHR21(i,j)*U_full)/100);
                s2 = s2+1;
            end
        end
    end
    
    
    %     [RoadSide_nodes] = RoadSide_sensor_nodes()
    %     DisplayCircle(sensor_x,sensor_y,bs);
    %     [sensors] = DisplaySensor(fin, N, CR, L, W);
    %     DisplayCircle1(bs,radius,sensor_x,sensor_y);
    
    
    %%%%%%%%%%%%% parameters for network partitioning  %%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%
    MaxNumberMV=ceil(0.3*length(bs));
    LowerTh=ceil(size(bs,1)/(MaxNumberMV)); % minimum no. of anchor points in a region
    UpperTH=ceil((8*size(bs,1))/(MaxNumberMV)); % maximum no. of anchor points in a region
    NumberOfMV=2;
    FitnessValueIGWO =[];
    while(NumberOfMV <= MaxNumberMV)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% apply I-GWO for partitioning the network  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        round1 = round
        [Fbest,PartitionAnchorPointIGWO,Convergence_curve] = I_GWO(NumberOfMV,bs,LowerTh,UpperTH, MaxNumberMV);
        FitnessValueIGWO(NumberOfMV-1) = Fbest;
        ChargingSequenceIGWO{round,NumberOfMV-1} = PartitionAnchorPointIGWO;
         
        %          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% apply PE-DFA for partitioning the network  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        %         [fitnessValue, PartitionAnchorPointDFA] = FAMAIN(NumberOfMV,bs,LowerTh,UpperTH, MaxNumberMV);
        %         FitnessValueDFA(round,NumberOfMV-1) = fitnessValue;
        %         ChargingSequenceDFA{round,NumberOfMV-1} = PartitionAnchorPointDFA;
        
        %         fprintf(ObjectiveFileIGWO,' %8.2f ', Fbest);
        %         fprintf(ObjectiveFileDFA,' %8.2f  ',fitnessValue);
 
        NumberOfMV=NumberOfMV+1
    end
%     [nmv, index] = sort(FitnessValueIGWO);
%     fprintf(NMV, '%d \n',index(1)+1);
  
    %     fprintf(ObjectiveFileIGWO,' \n ');
    %     fprintf(ObjectiveFileDFA,' \n ' );
    %%%%%%%%%%%%%%%%%%%%%%%%  For IGWO
    [Sort_FitnessValueIGWO, Index_FitnessValueIGWO] = sort(FitnessValueIGWO);
    ChargingSequence_IGWO1 = ChargingSequenceIGWO (round, Index_FitnessValueIGWO(1));
%     nvm = size(ChargingSequence_IGWO1{1}(:,:),1);
%      fprintf(NMV, '%d \n',nvm);
%       round = round+1;
% end
%     % %%%%%%%%%%%%%%%%%%%%%%  for DFA
%     [Sort_FitnessValueDFA, Index_FitnessValueDFA] = sort(FitnessValueDFA(round,:));
%     ChargingSequence_DFA1 = ChargingSequenceDFA (round, Index_FitnessValueDFA(1));
    %%%%%%%%%%%%%%%%%%%%%%partitioning data into region using K-means algorithm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [partnSPointMPSS]=Kmeans_clusteringMPSS(bs11,N);
    [partnSPointPMCDC]=Kmeans_clusteringPMCDC(hex_bs,N);
    [partnSPointHAN]=Kmeans_clusteringHAN(bs,N);
    [partnSPointMOAC]=Kmeans_clusteringMOAC(bs,N);

    %%%%%%%%%%%%%%%%%   Comparative analysis %%%%%%%%%%%%%%%%%%%%%%%%%
    [avgTimeTaken11,avgEnergyGained12, avgChargingTime13, avgDeadPeriod14, avgToatlDistance15, avgEUeff16, avgEnergyConsumed17, avgDeadNode18] = ProposedAlgorithmModified(ChargingSequence_IGWO1 ,BSX , BSY , v, bs, sojourn_time,T,DGR,CHR,ECR,RemainingEnergy,charging_time, data_collection_time, circle_nodes, Emax,U_full ,Pt,Data_generated,SCsensor,RemStorage,ett);    
    [avgTimeTaken21,avgEnergyGained22, avgChargingTime23, avgDeadPeriod24, avgToatlDistance25, avgEUeff26, avgEnergyConsumed27, avgDeadNode28] = MPSS(partnSPointMPSS  ,BSX , BSY , v, bs1, sojourn_time21,T,DGR21,CHR21,ECR21,RemainingEnergy21,charging_time21, data_collection_time21, circle_nodes11, Emax,U_full ,Pt,Data_generated,SCsensor,RemStorage21); % Boukerche 2021
    [avgTimeTaken31,avgEnergyGained32, avgChargingTime33, avgDeadPeriod34, avgToatlDistance35, avgEUeff36, avgEnergyConsumed37, avgDeadNode38] = PMCDC(partnSPointPMCDC , BSX , BSY , v, hex_bs, hex_sojourn_time,T,hex_DGR,hex_CHR,hex_ECR,hex_Energy,hex_charging_time, hex_data_collection_time, hex_circle_nodes, Emax,U_full ,Pt,hex_Data_generated,SCsensor,hex_RemStorage);% Lyu:2020
    [avgTimeTaken41,avgEnergyGained42, avgChargingTime43, avgDeadPeriod44, avgToatlDistance45, avgEUeff46, avgEnergyConsumed47, avgDeadNode48] = MOAC(partnSPointMOAC , BSX , BSY , v, bs, sojourn_time,T,DGR,CHR,ECR,RemainingEnergy,charging_time, data_collection_time, circle_nodes, Emax,U_full ,Pt,Data_generated,SCsensor,RemStorage); % Wei:2020
    [avgTimeTaken51,avgEnergyGained52, avgChargingTime53, avgDeadPeriod54, avgToatlDistance55, avgEUeff56, avgEnergyConsumed57, avgDeadNode58] = HAN2018(partnSPointHAN  ,BSX , BSY , v, bs1, sojourn_time21,T,DGR21,CHR21,ECR21,RemainingEnergy21,charging_time21, data_collection_time21, circle_nodes11, Emax,U_full ,Pt,Data_generated21,SCsensor,RemStorage21); % HAN:2018
    
    %%%%%%%%%%%%%%%% write data into file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    fprintf(TimeTaken1,'%8.4f %8.4f %8.4f %8.4f %8.4f  \n ',avgTimeTaken11,avgTimeTaken21,avgTimeTaken31,avgTimeTaken41,avgTimeTaken51);
    fprintf(EnergyGain1,'%8.4f %8.4f %8.4f %8.4f %8.4f  \n',avgEnergyGained12,avgEnergyGained22,avgEnergyGained32, avgEnergyGained42, avgEnergyGained52);
    fprintf(DeadPeriod1,'%8.4f %8.4f %8.4f %8.4f %8.4f  \n ',avgDeadPeriod14/avgDeadNode18,avgDeadPeriod24/avgDeadNode28, avgDeadPeriod34/avgDeadNode38, avgDeadPeriod44/avgDeadNode48, avgDeadPeriod54/avgDeadNode58);
    fprintf(TotalDistance1,'%8.4f %8.4f %8.4f %8.4f %8.4f  \n',avgToatlDistance15,avgToatlDistance25, avgToatlDistance35, avgToatlDistance45, avgToatlDistance55);
    fprintf(EUeff1,'%8.4f %8.4f %8.4f %8.4f %8.4f  \n ',avgEUeff16,avgEUeff26, avgEUeff36, avgEUeff46, avgEUeff56);
    fprintf(ChargeTime1,'%8.4f %8.4f %8.4f %8.4f %8.4f  \n ',avgChargingTime13, avgChargingTime23, avgChargingTime33, avgChargingTime43, avgChargingTime53);
    fprintf(EnergyConsumed1,'%8.4f %8.4f %8.4f %8.4f %8.4f  \n ',avgEnergyConsumed17, avgEnergyConsumed27, avgEnergyConsumed37, avgEnergyConsumed47, avgEnergyConsumed57);
    fprintf(deadSensor1,'%d %d %d %d %d  \n ',floor(avgDeadNode18), floor(avgDeadNode28), floor(avgDeadNode38), floor(avgDeadNode48), floor(avgDeadNode58));
   
    
    avgTimeTaken1(round,:) = [avgTimeTaken11/3600, avgTimeTaken21/3600, avgTimeTaken31/3600, avgTimeTaken41/3600, avgTimeTaken51/3600];
    avgEnergyGained1(round,:) = [avgEnergyGained12/1000, avgEnergyGained22/1000, avgEnergyGained32/1000, avgEnergyGained42/1000, avgEnergyGained52/1000];
    avgChargingTime1(round,:) = [avgChargingTime13/3600, avgChargingTime23/3600, avgChargingTime33/3600, avgChargingTime43/3600, avgChargingTime53/3600];
    avgDeadPeriod1(round,:) = [avgDeadPeriod14/avgDeadNode18, avgDeadPeriod24/avgDeadNode28, avgDeadPeriod34/avgDeadNode38, avgDeadPeriod44/avgDeadNode48, avgDeadPeriod54/avgDeadNode58];
    avgTDeadPeriod1(round,:) = [avgDeadPeriod14/3600, avgDeadPeriod24/3600, avgDeadPeriod34/3600, avgDeadPeriod44/3600, avgDeadPeriod54/3600];
    avgToatlDistance1(round,:) = [avgToatlDistance15, avgToatlDistance25, avgToatlDistance35, avgToatlDistance45, avgToatlDistance55];
    avgEUeff1(round,:) = [avgEUeff16, avgEUeff26, avgEUeff36, avgEUeff46, avgEUeff56];
    avgEnergyConsumed1(round,:) = [avgEnergyConsumed17/1000, avgEnergyConsumed27/1000, avgEnergyConsumed37/1000, avgEnergyConsumed47/1000, avgEnergyConsumed57/1000];
    avgDeadSensor1(round,:) = [avgDeadNode18, avgDeadNode28, avgDeadNode38, avgDeadNode48, avgDeadNode58];
   
    round = round+1;
end
%%%%%%%%%%%%%%% average result write in to the file %%%%%%%%%%%%%%%%%%%%%%%%
 
fprintf(TimeTaken,'%d %8.4f %8.4f %8.4f %8.4f %8.4f  \n',N,mean(avgTimeTaken1(:,1)),mean(avgTimeTaken1(:,2)), mean(avgTimeTaken1(:,3)), mean(avgTimeTaken1(:,4)), mean(avgTimeTaken1(:,5)));
fprintf(EnergyGain,'%d %8.4f %8.4f %8.4f %8.4f %8.4f  \n ',N, mean(avgEnergyGained1(:,1)),mean(avgEnergyGained1(:,2)), mean(avgEnergyGained1(:,3)), mean(avgEnergyGained1(:,4)), mean(avgEnergyGained1(:,5)));
fprintf(ChargeTime,'%d %8.4f %8.4f %8.4f %8.4f %8.4f  \n ',N,mean(avgChargingTime1(:,1)),mean(avgChargingTime1(:,2)), mean(avgChargingTime1(:,3)), mean(avgChargingTime1(:,4)), mean(avgChargingTime1(:,5)));
fprintf(DeadPeriod,'%d %8.4f %8.4f %8.4f %8.4f %8.4f  \n ',N,mean(avgDeadPeriod1(:,1)),mean(avgDeadPeriod1(:,2)), mean(avgDeadPeriod1(:,3)), mean(avgDeadPeriod1(:,4)), mean(avgDeadPeriod1(:,5)));
fprintf(TDeadPeriod,'%d %8.4f %8.4f %8.4f %8.4f %8.4f  \n ',N,mean(avgTDeadPeriod1(:,1)),mean(avgTDeadPeriod1(:,2)), mean(avgTDeadPeriod1(:,3)), mean(avgTDeadPeriod1(:,4)), mean(avgTDeadPeriod1(:,5)));
fprintf(TotalDistance,'%d %8.4f %8.4f %8.4f %8.4f %8.4f  \n ',N, mean(avgToatlDistance1(:,1)),mean(avgToatlDistance1(:,2)), mean(avgToatlDistance1(:,3)), mean(avgToatlDistance1(:,4)), mean(avgToatlDistance1(:,5)));
fprintf(EUeff,'%d %8.4f %8.4f %8.4f %8.4f %8.4f  \n ', N,mean(avgEUeff1(:,1)),mean(avgEUeff1(:,2)), mean(avgEUeff1(:,3)), mean(avgEUeff1(:,4)), mean(avgEUeff1(:,5)));
fprintf(EnergyConsumed,'%d %8.4f %8.4f %8.4f %8.4f %8.4f  \n ',N,mean(avgEnergyConsumed1(:,1)),mean(avgEnergyConsumed1(:,2)), mean(avgEnergyConsumed1(:,3)), mean(avgEnergyConsumed1(:,4)), mean(avgEnergyConsumed1(:,5)));
fprintf(deadSensor,'%d %d %d %d %d %d  \n ',N,floor(mean(avgDeadSensor1(:,1))),floor(mean(avgDeadSensor1(:,2))), floor(mean(avgDeadSensor1(:,3))), floor(mean(avgDeadSensor1(:,4))), floor(mean(avgDeadSensor1(:,5))));

%%%%%%%%%%%%%% closing file %%%%%%%%%%%%%%%%%%%%%%%%%%

fclose(TimeTaken);
fclose(EnergyGain);
fclose(DeadPeriod);
fclose(TotalDistance);
fclose(EUeff);
fclose(ChargeTime);
fclose(EnergyConsumed);
fclose(TimeTaken1);
fclose(EnergyGain1);
fclose(DeadPeriod1);
fclose(TotalDistance1);
fclose(EUeff1);
fclose(ChargeTime1);
fclose(EnergyConsumed1);

