function [AvailableEnergyIGWO1, TotalDeadPeriodIGWO, deadSensor] = InitialCycle( ChargingSequence_IGWO1,BSX, BSY, v, bs, sojourn_time,DGR,CRR,ECR,Energy,charging_time, data_collection_time, circle_nodes, MVX, MVY,Emax,U_full,Pt,Data_generated, SCsensor,RemStorage)
% size of charging sequence
GWOR = size(ChargingSequence_IGWO1,1);

%%%%%%% count the number of nodes in each circle
CurPosX = MVX;
CurPosY = MVY;
count = 0;
for i=1: size(circle_nodes,1)
    for j = 1: size(circle_nodes,2)
        if (circle_nodes(i,j) ~= 0)
            count = count +1;
            NSensor(i,1) = count;
        end
    end
    count = 0;
end

AvailableEnergyIGWO1 = zeros(size(circle_nodes,1), size(circle_nodes,2));
 %%% for IGWO
 EnergyIGWO = Energy;
 t1=0;
for i = 1: GWOR
    %%%%% from base station to first anchor points
    TotalTimeTaken = 0;
    TotalTime = 0;
    TravellTimeIGWO = [];
    flag = 0;
   % LifeTimeIGWO1 = zeros(GWOR,size(circle_nodes,2));
    ChargingSequence_IGWO2 =[];
    ChargingSequence_IGWO2 = cell2mat(ChargingSequence_IGWO1{1}(i,:));
    %%%%%%%%% sort them according to distance and remaining energy %%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [ChargingSequence_IGWO] = HeuristicApproach(ChargingSequence_IGWO2, Energy, BSX, BSY, bs, Emax, circle_nodes, SCsensor, RemStorage);
    DeadPeriodIGWO = [];
    %%%%%%%%%%%%%%%%% 
    for j = 1: length( ChargingSequence_IGWO)
        MVPositionIGWOX(j) = bs(ChargingSequence_IGWO(j),1); %%% current position of MV
        MVPositionIGWOY(j) = bs(ChargingSequence_IGWO(j),2);
        distIGWO(j) = sqrt((CurPosX - MVPositionIGWOX(j))^2 + (CurPosY - MVPositionIGWOY(j))^2);
        TravellTimeIGWO(j) = distIGWO(j)/v;
        ChargingTimeIGWO(j) = (sojourn_time(ChargingSequence_IGWO(j),1));
        TotalTimeTaken = TotalTimeTaken + TravellTimeIGWO(j);

       p = 0;
        for k =1: size(circle_nodes,2)  %%% charge and collect data from multiple sensor nodes at jth sojourn location. 
            if (circle_nodes(ChargingSequence_IGWO(j),k) ~= 0)
                p = p + 1;
                RemainingEnergyIGWO1(j,p) = (EnergyIGWO( ChargingSequence_IGWO(j),k)) - TotalTimeTaken * (ECR(ChargingSequence_IGWO(j),k));
                EnergyGainedIGWO1(j,k) = min(Emax, (sojourn_time(ChargingSequence_IGWO(j),1) * CRR(ChargingSequence_IGWO(j),k)) ) ;
                LifeTimeIGWO1(j,p) =  RemainingEnergyIGWO1(j,k)/(ECR(ChargingSequence_IGWO(j),1));
                DataGenerated(j,k) = (TotalTimeTaken) *(DGR(ChargingSequence_IGWO(j),k));
                RemStorageSensor(j,k) = (RemStorage( ChargingSequence_IGWO(j),k)) - DataGenerated(j,k);
                RemLifeStorage1(j,p) =  RemStorageSensor(j,k)/DGR(ChargingSequence_IGWO(j),k);
                AvailableEnergyIGWO1(ChargingSequence_IGWO(j),k) = min(Emax,(RemainingEnergyIGWO1(j,k) + EnergyGainedIGWO1(j,k)));
                
            end
        end
        RemainingEnergyIGWO(j) = min(RemainingEnergyIGWO1(j,p));
        EnergyGainedIGWO(j) = sum(EnergyGainedIGWO1(j,:));
        RemLifeStorage(j) = min(RemLifeStorage1(j,p));
        RemLifeEnergy(j) = min(LifeTimeIGWO1(j,p));
        LifeTimeIGWO(j) = min(RemLifeEnergy(j),RemLifeStorage(j));
        ArrivalTimeIGWO(j) =  (TotalTimeTaken) ; % TotalTimeTaken;
      
       
        %   for j = 2: length( ChargingSequence_IGWO)
        if (ArrivalTimeIGWO(j) <= LifeTimeIGWO(j))
            TotalTimeTaken = TotalTimeTaken +   ChargingTimeIGWO(j);
            CurPosX = MVPositionIGWOX(j);
            CurPosY = MVPositionIGWOY(j);
            diff =0;
        else
            success = 0;
            diff =0;
            flag = 1;
%             MVPositionIGWOX(j) = bs(ChargingSequence_IGWO(j),1);
%             MVPositionIGWOY(j) = bs(ChargingSequence_IGWO(j),2);
%             distIGWO(j) = sqrt((CurPosX - MVPositionIGWOX(j))^2 + (CurPosY - MVPositionIGWOY(j))^2);
           % TravellTimeIGWO(j) = distIGWO(j)/v;
           % ChargingTimeIGWO(j) = (sojourn_time(ChargingSequence_IGWO(j),1));
%             TotalTimeTaken = TotalTimeTaken +  TravellTimeIGWO(j) +  ChargingTimeIGWO(j);
            p = 0;
            for k =1: size(circle_nodes,2)
                if (circle_nodes(ChargingSequence_IGWO(j),k) ~= 0)
                    if (ArrivalTimeIGWO(j) > LifeTimeIGWO1(j,k))
                        %%%% % Do Partial charging of the sensors and check
                        p = p +1;
                        lambda =0.9;
                        A2 = charging_time(ChargingSequence_IGWO(1:j),: );
                        [max1,index1] = max( charging_time(ChargingSequence_IGWO(1:j),: ));   %% select anchor           % Select the sensor with longest charging time
                        [~,index2] = max(max1);
                        index =index1(index2);
                        q = 0.2*charging_time(ChargingSequence_IGWO(index),index2);
                        while( charging_time(ChargingSequence_IGWO(index),index2) > 0.2*charging_time(ChargingSequence_IGWO(index),index2) && lambda >= 0.2)                            % remaining energy is more than 2*delta then reduce the charging time by delta
                            
                            charging_time(ChargingSequence_IGWO(index),index2) = max(0, charging_time(ChargingSequence_IGWO(index),index2) -  lambda *charging_time(ChargingSequence_IGWO(index),index2));       % Reduce the charging time by ett
                            data_collection_time(ChargingSequence_IGWO(index),index2) = data_collection_time(ChargingSequence_IGWO(index),index2) ;       % Reduce the charging time by ett
                            sojourn_time(ChargingSequence_IGWO(index),1) = max( charging_time(ChargingSequence_IGWO(index),index2), data_collection_time(ChargingSequence_IGWO(index),index2));
                            TotalTimeTaken = TotalTimeTaken  - lambda * charging_time(ChargingSequence_IGWO(index),index2); % remove delta time into time period
                            RemainingEnergyIGWO1(index,index2) = (EnergyIGWO( ChargingSequence_IGWO(index),index2)) - (TotalTimeTaken * (ECR(ChargingSequence_IGWO(index),index2)));
                            LifeTimeIGWO1(index,index2) =  RemainingEnergyIGWO1(index,index2)/(ECR(ChargingSequence_IGWO(index),1));
                            DataGenerated(index,index2) = TotalTimeTaken * DGR(ChargingSequence_IGWO(index),index2);
                            RemStorageSensor(index,index2) = (RemStorage( ChargingSequence_IGWO(j),k)) - DataGenerated(index,index2);
                            RemLifeStorage1(index,index2) =  RemStorageSensor(index,index2)/DGR(ChargingSequence_IGWO(index),index2);
                            AvailableEnergyIGWO1(ChargingSequence_IGWO(index),index2) = min(Emax,(RemainingEnergyIGWO1(index,index2) + EnergyGainedIGWO1(index,index2)));
                          % TotalTimeTaken;
                            RemainingEnergyIGWO(index) = min(RemainingEnergyIGWO1(index,:));
                            EnergyGainedIGWO(index) = sum(EnergyGainedIGWO1(index,:));
                            RemLifeStorage(index) = min(RemLifeStorage1(index,:));
                            RemLifeEnergy(index) = min(LifeTimeIGWO1(index,:));
                            LifeTimeIGWO(index) = min(RemLifeEnergy(index),RemLifeStorage(index));
                            ArrivalTimeIGWO(index) =  (TotalTimeTaken); 
                            
                            for t = index:j   % Update the arival time
                                RemainingEnergyIGWO(t) = min(RemainingEnergyIGWO1(t,p));
                                EnergyGainedIGWO(t) = sum(EnergyGainedIGWO1(t,:));
                                RemLifeStorage(t) = min(RemLifeStorage1(t,p));
                                RemLifeEnergy(t) = min(LifeTimeIGWO1(t,p));
                                LifeTimeIGWO(t) = min(RemLifeEnergy(t),RemLifeStorage(t));
                                ArrivalTimeIGWO(t) =  (TotalTimeTaken );
                                
                            end
                            
                            if (ArrivalTimeIGWO(j) < LifeTimeIGWO1(j,k))
                                success=1;
                                CurPosX = MVPositionIGWOX(j);
                                CurPosY = MVPositionIGWOY(j);
                                TotalTimeTaken = TotalTimeTaken +   TravellTimeIGWO(j) +  ChargingTimeIGWO(j);
                                break;
                            end
                            [max1,index1] = max( charging_time(ChargingSequence_IGWO(1:j),: ));   %% select anchor           % Select the sensor with longest charging time
                            [~,index2] = max(max1);
                            index =index1(index2);
                            lambda = lambda - 0.1;
                        end
                    end
                    if success == 0
                        diff1 = abs(ArrivalTimeIGWO(j) - LifeTimeIGWO1(j,k));
                        diff = diff + diff1;
                        t1 = t1+1;
                        deadSensor (t1) = circle_nodes(ChargingSequence_IGWO(j),k);
                        CurPosX = MVPositionIGWOX(j);
                        CurPosY = MVPositionIGWOY(j);
                        TotalTimeTaken = TotalTimeTaken +   TravellTimeIGWO(j) +  ChargingTimeIGWO(j);
                    end
                    RemainingEnergyIGWO1(j,k) = (EnergyIGWO( ChargingSequence_IGWO(j),k)) - (TotalTimeTaken)*(ECR(ChargingSequence_IGWO(j),k));
                    EnergyGainedIGWO1(j,k) = min(Emax, (sojourn_time(ChargingSequence_IGWO(j),1)) * CRR(ChargingSequence_IGWO(j),k)) ;
                    LifeTimeIGWO1(j,k) =  RemainingEnergyIGWO1(j,k)/(ECR(ChargingSequence_IGWO(j),k));
                    DataGenerated(j,k) = (TotalTimeTaken) *(DGR(ChargingSequence_IGWO(j),k));
                    RemStorageSensor(j,k) = SCsensor - DataGenerated(j,k);
                    RemLifeStorage1(j,k) =  RemStorageSensor(j,k)/DGR(ChargingSequence_IGWO(j),k);
                    AvailableEnergyIGWO1(ChargingSequence_IGWO(j),k) = min(Emax,(RemainingEnergyIGWO1(j,k) + EnergyGainedIGWO1(j,k)));
                end
            end
            RemainingEnergyIGWO(j) = min(RemainingEnergyIGWO1(j,p));
            EnergyGainedIGWO(j) = sum(EnergyGainedIGWO1(j,:));
            RemLifeStorage(j) = min(RemLifeStorage1(j,p));
            RemLifeEnergy(j) = min(LifeTimeIGWO1(j,p));
            LifeTimeIGWO(j) = min(RemLifeEnergy(j),RemLifeStorage(j));
            ArrivalTimeIGWO(j) = TotalTimeTaken;
            DeadPeriodIGWO(j) = diff;
        end
         if (j == length( ChargingSequence_IGWO))  %%% from last position to base station 
%             MVPositionIGWOX(j+1) = bs(ChargingSequence_IGWO(1),1);
%             MVPositionIGWOY(j+1) = bs(ChargingSequence_IGWO(1),2);
            distIGWO(j+1) = sqrt((MVPositionIGWOX(j) - MVX)^2 + (MVPositionIGWOY(j) - MVY)^2);
            TravellTimeIGWO(j+1) = distIGWO(j+1)/v;
            TotalTime =  TotalTimeTaken + TravellTimeIGWO(j+1) + ChargingTimeIGWO(j);
         end      
    end
    TotalTimeTakenIGWO(i,1) = TotalTime;
    TotalEnergyGainedIGWO(i,1) = sum(EnergyGainedIGWO);
    if (flag ==1)
         TotalDeadPeriodIGWO(i,1) = sum(DeadPeriodIGWO);
    else
         TotalDeadPeriodIGWO(i,1) = 0;
    end
    TotalDeadPeriodIGWO(i,1) = TotalDeadPeriodIGWO(i,1);
    TotalChargingTimeIGWO(i,1) = sum(ChargingTimeIGWO);
    TotalDistanceIGWO(i,1) = sum(distIGWO);
    EnergyUtilIGWO(i,1) = TotalEnergyGainedIGWO(i,1)/(TotalChargingTimeIGWO(i,1) * U_full + (TotalDistanceIGWO(i,1)/v)*Pt);
    TotalEnergyConsumedIGWO(i,1) = (TotalChargingTimeIGWO(i,1) * U_full + (TotalDistanceIGWO(i,1)/v)*Pt);
end
EnergyIGWO = AvailableEnergyIGWO1;
% 
%   % update energy
% %%%%%%%%%%%%%%% IGWO %%%%%%%%%%%%%%%%
% avgTimeTakenIGWO1(cycle) = mean(TotalTimeTakenIGWO);
% avgEnergyGainedIGWO1(cycle) = mean(TotalEnergyGainedIGWO);
% avgChargingTimeIGWO1(cycle) = mean(TotalChargingTimeIGWO);
% avgDeadPeriodIGWO1(cycle) = mean(TotalDeadPeriodIGWO);
% avgToatlDistanceIGWO1(cycle) = mean(TotalDistanceIGWO);
% avgEUIGWO1(cycle) = mean(EnergyUtilIGWO);
% avgEnergyConsumedIGWO1(cycle) = mean(TotalEnergyConsumedIGWO);
% 
%  cycle = cycle +1;
% end
% avgTimeTaken = mean(avgTimeTakenIGWO1);
% avgEnergyGained = mean(avgEnergyGainedIGWO1);
% avgChargingTime = mean(avgChargingTimeIGWO1);
% avgDeadPeriod = mean(avgDeadPeriodIGWO1);
% avgToatlDistance = mean(avgToatlDistanceIGWO1);
% avgEUeff = mean(avgEUIGWO1);
% avgEnergyConsumed = mean(avgEnergyConsumedIGWO1);

end
 
 
 
 
 
 
 
 
%  %%%%%%%%%%%%%%%%%%
%  t = 0;
% for i = 1: GWOR
%     TotalTimeTaken = 0;
%     TotalTime = 0;
%     flag = 0;
%     ChargingSequence_IGWO2 =[];
%     ChargingSequence_IGWO2 = cell2mat(ChargingSequence_IGWO1{1}(i,:));
%     [ChargingSequence_IGWO] = HeuristicApproach(ChargingSequence_IGWO2, EnergyIGWO,BSX,BSY,bs, Emax,circle_nodes,SCsensor,RemStorage);
%     MVPositionIGWOX(1) = bs(ChargingSequence_IGWO(1),1); %%% current position of MV
%     MVPositionIGWOY(1) = bs(ChargingSequence_IGWO(1),2);
%     distIGWO(1) = sqrt((CurPos(1,1) - MVPositionIGWOX(1))^2 + (CurPos(1,2) - MVPositionIGWOY(1))^2);
%     TravellTimeIGWO(1) = distIGWO(1)/v;
%     TotalTimeTaken = TotalTimeTaken + TravellTimeIGWO(1);
%     ArrivalTimeIGWO(1) =   TravellTimeIGWO(1);
%     ChargingTimeIGWO(1) = (sojourn_time(ChargingSequence_IGWO(1),1));
%    
%     for k =1: size(circle_nodes,2)
%         if (circle_nodes(ChargingSequence_IGWO(1),k) ~= 0)
%             RemainingEnergyIGWO1(1,k) = max(0,(Energy( ChargingSequence_IGWO(1),k)) - TotalTimeTaken *(ECR(ChargingSequence_IGWO(1),k)));
%             EnergyGainedIGWO1(1,k) = min(Emax, (sojourn_time(ChargingSequence_IGWO(1),1)) * CRR(ChargingSequence_IGWO(1),1) ) ;
%             DataGenerated(1,k) = TotalTimeTaken*DGR(ChargingSequence_IGWO(1),k);
%             RemainingStorage(1,k) = SCsensor -  DataGenerated(1,k); 
%             LifeTimeIGWO1(1,k) = min( RemainingStorage(1,k)/DGR(ChargingSequence_IGWO(1),k) , RemainingEnergyIGWO1(1,k)/(ECR(ChargingSequence_IGWO(1),k)));
%             AvailableEnergyIGWO1( ChargingSequence_IGWO(1),k) = min(Emax,(RemainingEnergyIGWO1(1,k) + EnergyGainedIGWO1(1,k)));
%         end
%     end
%     RemainingEnergyIGWO(1) = min(RemainingEnergyIGWO1(1,:));
%     EnergyGainedIGWO(1) = sum(EnergyGainedIGWO1(1,:));
%     LifeTimeIGWO(1) = min(LifeTimeIGWO1(1,:));  
%     for j = 2: length( ChargingSequence_IGWO)
%         if (ArrivalTimeIGWO(j-1) <= LifeTimeIGWO(j-1))
%             MVPositionIGWOX(j) = bs(ChargingSequence_IGWO(j),1);
%             MVPositionIGWOY(j) = bs(ChargingSequence_IGWO(j),2);
%             distIGWO(j) = sqrt((MVPositionIGWOX(j) - MVPositionIGWOX(j-1))^2 + (MVPositionIGWOY(j) - MVPositionIGWOY(j-1))^2);
%             TravellTimeIGWO(j) = distIGWO(j)/v;
%             ChargingTimeIGWO(j) = (sojourn_time(ChargingSequence_IGWO(j),1));
%             TotalTimeTaken = TotalTimeTaken +   TravellTimeIGWO(j) +  ChargingTimeIGWO(j-1);
%             for k =1: size(circle_nodes,2)
%                 if (circle_nodes(ChargingSequence_IGWO(j),k) ~= 0)
%                     RemainingEnergyIGWO1(j,k) = max(0,(Energy( ChargingSequence_IGWO(j),k)) - TotalTimeTaken * (ECR(ChargingSequence_IGWO(j),k)));
%                     EnergyGainedIGWO1(j,k) = min(Emax, (sojourn_time(ChargingSequence_IGWO(j),1) * CRR(ChargingSequence_IGWO(j),k)) ) ;
%                     DataGenerated(j,k) = TotalTimeTaken*DGR(ChargingSequence_IGWO(j),k);
%                     RemainingStorage(j,k) = SCsensor -  DataGenerated(j,k);
%                     LifeTimeIGWO1(j,k) = min( RemainingStorage(j,k)/DGR(ChargingSequence_IGWO(j),k) , RemainingEnergyIGWO1(j,k)/(ECR(ChargingSequence_IGWO(j),k)));
%                     AvailableEnergyIGWO1(ChargingSequence_IGWO(j),k) = min(Emax,(RemainingEnergyIGWO1(j,k) + EnergyGainedIGWO1(j,k)));
%                 end
%             end
%             RemainingEnergyIGWO(j) = min(RemainingEnergyIGWO1(j,:));
%             EnergyGainedIGWO(j) = sum(EnergyGainedIGWO1(j,:));
%             LifeTimeIGWO(j) = min(LifeTimeIGWO1(j,:));
%             ArrivalTimeIGWO(j) =   TotalTimeTaken;
%              diff =0;
%         else
%            diff =0;
%            flag = 1;
%             MVPositionIGWOX(j) = bs(ChargingSequence_IGWO(j),1);
%             MVPositionIGWOY(j) = bs(ChargingSequence_IGWO(j),2);
%             distIGWO(j) = sqrt((MVPositionIGWOX(j) - MVPositionIGWOX(j-1))^2 + (MVPositionIGWOY(j) - MVPositionIGWOY(j-1))^2);
%             TravellTimeIGWO(j) = distIGWO(j)/v;
%             ChargingTimeIGWO(j) = (sojourn_time(ChargingSequence_IGWO(j),1));
%             TotalTimeTaken = TotalTimeTaken +  TravellTimeIGWO(j) +  ChargingTimeIGWO(j-1);
%             for k =1: size(circle_nodes,2)
%                 if (circle_nodes(ChargingSequence_IGWO(j),k) ~= 0)
%                     if (ArrivalTimeIGWO(j-1) > LifeTimeIGWO1(j-1))
%                         diff1 = abs(ArrivalTimeIGWO(j-1) - LifeTimeIGWO1(j-1));
%                         diff = diff + diff1;
%                         t = t+1;
%                         deadSensor (t) = circle_nodes(j,k);
%                       
%                         
%                     end
%                     RemainingEnergyIGWO1(j,k) = max(0,(Energy( ChargingSequence_IGWO(j),k)) - TotalTimeTaken*(ECR(ChargingSequence_IGWO(j),k)));
%                     EnergyGainedIGWO1(j,k) = min(Emax, (sojourn_time(ChargingSequence_IGWO(j),1)) * CRR(ChargingSequence_IGWO(j),k)) ;
%                     DataGenerated(j,k) = TotalTimeTaken*DGR(ChargingSequence_IGWO(j),k);
%                     RemainingStorage(j,k) = SCsensor -  DataGenerated(j,k);
%                     LifeTimeIGWO1(j,k) = min( RemainingStorage(j,k)/DGR(ChargingSequence_IGWO(j),k) , RemainingEnergyIGWO1(j,k)/(ECR(ChargingSequence_IGWO(j),k)));
%                     AvailableEnergyIGWO1(ChargingSequence_IGWO(j),k) = min(Emax,(RemainingEnergyIGWO1(j,k) + EnergyGainedIGWO1(j,k)));
%                 end
%             end
%             RemainingEnergyIGWO(j) = min(RemainingEnergyIGWO1(j,:));
%             EnergyGainedIGWO(j) = sum(EnergyGainedIGWO1(j,:));
%             LifeTimeIGWO(j) = min(LifeTimeIGWO1(j,:));
%             ArrivalTimeIGWO(j) = TotalTimeTaken;
%             DeadPeriodIGWO(j) = diff;
%         end
%          if (j == length( ChargingSequence_IGWO))
%             MVPositionIGWOX(j+1) = bs(ChargingSequence_IGWO(1),1);
%             MVPositionIGWOY(j+1) = bs(ChargingSequence_IGWO(1),2);
%             distIGWO(j+1) = sqrt((MVPositionIGWOX(j) - MVPositionIGWOX(j+1))^2 + (MVPositionIGWOY(j) - MVPositionIGWOY(j+1))^2);
%             TravellTimeIGWO(j+1) = distIGWO(j+1)/v;
%             TotalTime =  TotalTimeTaken + TravellTimeIGWO(j+1) + ChargingTimeIGWO(j);
%          end      
%     end
%     TotalTimeTakenIGWO(i,1) = TotalTime;
%     TotalEnergyGainedIGWO(i,1) = sum(EnergyGainedIGWO);
%     if (flag ==1)
%          TotalDeadPeriodIGWO(i,1) = sum(DeadPeriodIGWO);
%     else
%          TotalDeadPeriodIGWO(i,1) = 0;
%     end
%     TotalChargingTimeIGWO(i,1) = sum(ChargingTimeIGWO);
%     TotalDistanceIGWO(i,1) = sum(distIGWO);
%     EnergyUtilIGWO(i,1) = TotalEnergyGainedIGWO(i,1)/(TotalChargingTimeIGWO(i,1) * U_full + (TotalDistanceIGWO(i,1)/v)*Pt);
%     TotalEnergyConsumedIGWO(i,1) = (TotalChargingTimeIGWO(i,1) * U_full + (TotalDistanceIGWO(i,1)/v)*Pt);
% end         
% T = max(TotalChargingTimeIGWO);
% end