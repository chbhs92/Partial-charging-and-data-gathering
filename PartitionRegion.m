function [partitionNode4, PartitionAngle1] = PartitionRegion(partitionNode1, PartitionAngle, Thresold, theta, bs, Pt,v,EnergyCharge)
BSX = 100;
BSY = 100;
Dmax = BSX * 1.414;
lambda = 0.9;
 for i=1:size(partitionNode1,1)
    [ TotalEnergy(i)] = GetSchedule(bs, partitionNode1(i,:),Pt,v, EnergyCharge, theta);
 end
[NumberAnchorPoint] = CalculateNumberOfAnchorPoints(partitionNode1);
 t=1;
for i=1:size(partitionNode1,1)
    if (TotalEnergy(i) > 2*Thresold)
              for j=1:size(partitionNode1,2)
                  if (partitionNode1(i,j) ~= 0)
                      ElementAngle(t,j) = theta(partitionNode1(i,j));  % find region who has greater energy than upper thresold 
                      partitionIndex(t) = i; % partitioning region
                      DistanceToAnchor(t,j) = sqrt((BSX - bs(partitionNode1(i,j),1))^2 + (BSY - bs(partitionNode1(i,j),2))^2);
                      WeightedTerm(t,j) = lambda * ElementAngle(t,j)/360 + (1-lambda) * DistanceToAnchor(t,j)/Dmax;
                  end
              end
               t=t+1;
              
    end    
end

[NumberElement] = CalculateNumberOfAnchorPoints(ElementAngle);  %%%% find number of anchor points
midElement=ceil(NumberElement/2);  % middle element
partitionNode2=partitionNode1;    %%% replica 
s=1;
for i=1:size(ElementAngle,1)
    [partitionElement(i,:), elemntIndex(i,:)] = sort(WeightedTerm(i,:));   % sort the anchor point ( big region ) with angle
     for j = 1:size(ElementAngle,2) 
         if ( partitionElement(i,j) ~= 0)   %%% removing zero items from sorted elements
             SortedElement(i,s) = partitionElement(i,j);   
             SortedIndex(i,s) = elemntIndex(i,j);
             s=s+1;
         end
         
     end
     s=1;  
end
% %%%%%%%%%  enter partition anchor point into another region
% [PartitionAnchorPoint, PartitionLine] = WeightedPartition(SortedElement, SortedIndex, bs, PartitionAngle, Thresold, partitionNode1);

len=size(partitionNode1,2);
r=size(partitionNode1,1);
p=1;
w=1;
flag = 0;
partitionNode3=zeros(r,len);
TotalEnergyConsume =0;
for i=1:size(partitionNode1,1)
    TotalEnergyConsume =0;
    if (find(i == partitionIndex))
        k=find(i == partitionIndex);
        MVX = BSX;
        MVY = BSY;
        for j=1:size(partitionNode1,2)
            if (partitionNode1(i,j) ~= 0)
                distance(i,j) = sqrt((MVX - bs(partitionNode1(i,SortedIndex( k,j)),1))^2 + (MVY - bs(partitionNode1(i,SortedIndex( k,j)),2))^2);
                travEnergy(i,j) = (distance(i,j)/v)*Pt;
                ChargeEnergy(i,j) = EnergyCharge(partitionNode1(i,SortedIndex( k,j)));
                TotalEnergyConsume = TotalEnergyConsume + travEnergy(i,j) + ChargeEnergy(i,j);
                if ( TotalEnergyConsume < 2*Thresold )
                    partitionNode3(p,j) = partitionNode1(i,SortedIndex( k,j));
                    MVX =  bs(partitionNode1(i,SortedIndex( k,j)),1);
                    MVY = bs(partitionNode1(i,SortedIndex( k,j)),2);
                else
                   % p=p+1;
                   flag = 1;
                    partitionNode3(p+1,w) = partitionNode1(i,SortedIndex(k ,j));
                    MVX =  bs(partitionNode1(i,SortedIndex( k,j)),1);
                    MVY = bs(partitionNode1(i,SortedIndex( k,j)),2);
                    w=w+1;
                end
            end
        end
        w=1;
          if flag == 1
           p = p+2;
          else 
            p = p+1;
          end
    else
        partitionNode3(p,:) = partitionNode1(i,:);
        p = p+1;
    end
end
g=1;
%%%%%% remove the zero from each region
for i = 1: size(partitionNode3,1)
    for j = 1: size(partitionNode3,2)
        if (partitionNode3(i,j) ~= 0)
            partitionNode4(i,g) = partitionNode3(i,j);
            g = g+1;
        end
    end
    g=1;
end

%%%%%%%%%%%%%%%%%% add  partition angle to original angle set

len1 = length(PartitionAngle);
PartitionAngle1=zeros(1,len1); %% replica
p=1;
for i=1:len1
    if (find(i == partitionIndex))
          k = find( i == partitionIndex );
          PartitionAngle1(p)= PartitionAngle(i);
          PartitionAngle1(p+1) = SortedElement(k, midElement(k));
          p=p+2;
    else
         PartitionAngle1(p)= PartitionAngle(i);
         p=p+1;
    end
   len1 = length(PartitionAngle1); 
end

end