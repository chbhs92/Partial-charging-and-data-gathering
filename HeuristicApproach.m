function  [ChargingSeq,dist_bs] = HeuristicApproach(ChargingSequence, RemEnergy,BSX,BSY,bs, Emax1,circle_nodes,SCsensor1,RemStorage)
%%%%%%%%%%%% calculate distance from base station %%%%%%%%%%%%%%%%%%%%%%%%
Emax = Emax1;
SCsensor = SCsensor1;
 theta = 0.7; % selection percentage between distance and Remaining 
 phi = 0.6; % 

for i = 1: length(ChargingSequence)
    dist_bs(i) = sqrt((BSX-bs(ChargingSequence(i),1))^2 + (BSY-bs(ChargingSequence(i),2))^2);
end
maxDistance = max(dist_bs);
%%%%%%%%%%%%%%%%%%% selection weitghted term %%%%%%%%%%%%%%%%%%%%%%%%%%%%
s=1;
for i = 1: length(ChargingSequence)
 for k =1: size(circle_nodes,2)
        if (circle_nodes(ChargingSequence(i),k) ~= 0)
            Energy(i,s) = RemEnergy(ChargingSequence(i),k);  % Remaining energy
            DataStore(i,s) = RemStorage(ChargingSequence(i),k); % Remaining storage
            s = s+1;
        end
 end
 s=1;
 minEnergy(i) = min(Energy(i,:));
 minStorage(i) = min(DataStore(i,:));
end
for i = 1: length(ChargingSequence)
    wightedTerm(i) =(1-theta)* dist_bs(i)/maxDistance + theta*(phi*minEnergy(i)/Emax +(1-phi)*minStorage(i)/SCsensor);
end
[element, imdex] = sort(wightedTerm);
for i = 1: length(wightedTerm)
    ChargingSeq(i) = ChargingSequence(imdex(i));
end
end