function [partitionNode2,partitionDegree] = mergeAngle(partitionNode1,theta,Thresold,Degree, PartitionAngle,bs,Pt,v,EnergyCharge)
s=1;
flag=0;
mergeIndex=0;
MergeWith=0;
partitionNode = partitionNode1;
len = size(partitionNode,1);
for i=1:len
    %%%% count number of anchor points in eaach partition
  [NumberAnchorPoint] = CalculateNumberOfAnchorPoints(partitionNode);
  for k=1:size(partitionNode,1)
    [ TravelEnergy(k)] = GetSchedule(bs, partitionNode(k,:),Pt,v,EnergyCharge, theta);
  end
    %%% Case 1 Assign a MV in each partition
    if (TravelEnergy(i) >= Thresold && TravelEnergy(i) <= 2 * Thresold)
        partitionNode(i,:)= partitionNode(i,:);
    end
    %%%% merge the partition q(1,1:6)=[ q(1,1:3) w(1,1:3)]
    totalMergeElement =0;
    if (i<size(partitionNode,1) && i>1)  %% Except 1st and last partition region
        if (TravelEnergy(i) < Thresold && TravelEnergy(i) ~=0)
            if (TravelEnergy(i-1) < TravelEnergy(i+1))
                k=1;
                for j=1:size(partitionNode,2)
                    if (partitionNode(i-1,j) == 0 && j <= size(partitionNode,2))
                        partitionNode(i-1,j) = partitionNode(i,k);
                        partitionNode(i,k) = 0;
                        mergeIndex(s) = i;
                        MergeWith(s) = i-1;
                        totalMergeElement = k;
                        k=k+1;
                    elseif (j == size(partitionNode,2) && (NumberAnchorPoint(i) - totalMergeElement) >0)
                        q=1;
                        for t=1:(NumberAnchorPoint(i) - totalMergeElement)
                            partitionNode(i-1,j+q) = partitionNode(i,k);
                            partitionNode(i,k) = 0;
                            mergeIndex(s) = i;
                            MergeWith(s) = i-1;
                            k=k+1;
                            q=q+1;
                        end
                    end
                end
                k=1;
                s=s+1;
            else   %if (NumberAnchorPoint(i-1,1) > NumberAnchorPoint(i+1,1))
                  k=1;
                for j=1:size(partitionNode,2)
                     if (partitionNode(i+1,j) == 0 && j <= size(partitionNode,2))
                       partitionNode(i+1,j) = partitionNode(i,k);
                       partitionNode(i,k) = 0;
                       mergeIndex(s) = i;
                       MergeWith(s) = i+1;
                       totalMergeElement = k;
                       k=k+1;
                     elseif (j == size(partitionNode,2) && (NumberAnchorPoint(i) - totalMergeElement) >0)
                          q=1;
                          for t=1:(NumberAnchorPoint(i) - totalMergeElement)
                              partitionNode(i+1,j+q) = partitionNode(i,k);
                              partitionNode(i,k) = 0;
                              mergeIndex(s) = i;
                              MergeWith(s) = i+1;
                              k=k+1;
                              q=q+1;
                          end
                     end
                 end
                
                k=1;
               s=s+1;
            end
           % s=s+1;
            % RandomAngle=randi([Degree(1,i),Degree(1,i+1)],1,1);
        end
    elseif (i==length(TravelEnergy) && TravelEnergy(i) ~=0)
        if (TravelEnergy(i) < Thresold && TravelEnergy(i) ~=0)
            if (TravelEnergy(1) < TravelEnergy(i-1) && flag == 0)
                k=1;
                for j=1:size(partitionNode,2)
                    if (partitionNode(1,j)==0 && j <= size(partitionNode,2))
                        partitionNode(1,j)=partitionNode(i,k);
                        partitionNode(i,k)=0;
                        mergeIndex(s)=i;
                        MergeWith(s) = 1;
                        totalMergeElement = k;
                        k=k+1;
                    elseif (j == size(partitionNode,2)&& (NumberAnchorPoint(i) - totalMergeElement) >0)
                          q=1;
                          for t=1:(NumberAnchorPoint(i) - totalMergeElement)
                              partitionNode(1,j+q)=partitionNode(i,k);
                              partitionNode(i,k)=0;
                              mergeIndex(s)=i;
                              MergeWith(s) = 1;
                              k=k+1;
                              q=q+1;
                          end
                    end
                end
                k=1;
                s=s+1;
            else
                  k=1;
                for j=1:size(partitionNode,2)
                    if (partitionNode(i-1,j)==0 && j <= size(partitionNode,2))
                        partitionNode(i-1,j)=partitionNode(i,k);
                        partitionNode(i,k)=0;
                        mergeIndex(s)=i;
                        MergeWith(s) = i-1;
                        totalMergeElement = k;
                        k=k+1;
                    elseif (j == size(partitionNode,2)&& (NumberAnchorPoint(i) - totalMergeElement) >0)
                        q=1;
                        for t=1:(NumberAnchorPoint(i) - totalMergeElement)
                            partitionNode(i-1,j+q) = partitionNode(i,k);
                            partitionNode(i,k)=0;
                            mergeIndex(s)=i;
                            MergeWith(s) = 1;
                            k=k+1;
                            q=q+1;
                        end
                    end
                end
                k=1;
                s=s+1;
                
            end
            % RandomAngle=randi([Degree(1,i),Degree(1,i+1)],1,1);
        end
    elseif (i==1 && TravelEnergy(i) ~=0 ) %%% for 1st partition region
        if (TravelEnergy(1) < Thresold && TravelEnergy(i) ~=0)
               flag=1;
            if (TravelEnergy(length(NumberAnchorPoint)) < TravelEnergy(i+1))
                k=1;
                for j=1:size(partitionNode,2) 
                    if (partitionNode(length(NumberAnchorPoint),j)==0 && j <= size(partitionNode,2))
                        partitionNode(length(NumberAnchorPoint),j) = partitionNode(1,k);
                        partitionNode(1,k)=0;
                        mergeIndex(s)=i;
                        MergeWith(s) = length(NumberAnchorPoint);
                        %partitionNode(1,j)=0
                        totalMergeElement = k;
                        k=k+1;
                    elseif ( j == size(partitionNode,2) && (NumberAnchorPoint(i) - totalMergeElement) >0)
                        q=1;
                        for t=1:(NumberAnchorPoint(i) - totalMergeElement)
                            partitionNode(length(NumberAnchorPoint),j+q) = partitionNode(1,k);
                            partitionNode(1,k)=0;
                            mergeIndex(s)=i;
                            MergeWith(s) = length(NumberAnchorPoint);
                            k=k+1;
                            q=q+1;
                        end
                    end
                end
                k=1;
                s=s+1;
            else
                  k=1;
                for j=1:size(partitionNode,2)
                    if (partitionNode(i+1,j) == 0 && j <= size(partitionNode,2))
                        partitionNode(i+1,j) = partitionNode(i,k);
                        partitionNode(i,k)=0;
                        mergeIndex(s)=i;
                        MergeWith(s) = i+1;
                        totalMergeElement = k;
                        k=k+1;
                        
                    elseif ( j == size(partitionNode,2) && (NumberAnchorPoint(i) - totalMergeElement) >0)
                        q=1;
                        for t=1:(NumberAnchorPoint(i) - totalMergeElement)
                            partitionNode(i+1,j+t) = partitionNode(i,k);
                            partitionNode(i,k)=0;
                            mergeIndex(s)=i;
                            MergeWith(s) = i+1;
                            k=k+1;
                        end
                    end
                end
                k=1;
                s=s+1;
            end
            % RandomAngle=randi([Degree(1,i),Degree(1,i+1)],1,1);
        end
    end
    h=1;
    partitionNode3 =[];
    for p=1:size(partitionNode,1)
        % for j=1:size(partitionNode,2)
        if ( partitionNode(p,1) ~= 0)
            partitionNode3(h,:)= partitionNode(p,:);
            h=h+1;
        end
        
    end
    partitionNode = partitionNode3;
    len = size(partitionNode,1);
    if (len <= i)
        break;
    end
end
%%%%%% delete the merge angle
if(length(mergeIndex) >=1 && mergeIndex(1) ~= 0)
    for i = 1:length(mergeIndex)
        if ( mergeIndex(i) > MergeWith(i))
            PartitionAngle( mergeIndex(i))=0;
        else
            PartitionAngle(MergeWith(i))=0;
        end
    end
end
k=1;
%partitionNode1=partitionNode; 
for i=1:size(partitionNode,1)
   % for j=1:size(partitionNode,2)
        if ( partitionNode(i,1) ~= 0)
            partitionNode2(k,:)= partitionNode(i,:);
            k=k+1;
        end
  
   % end
    
end
p=1;
for i=1:length(PartitionAngle)
    if (PartitionAngle(i) ~= 0)
        partitionDegree(p) = PartitionAngle(i);
        p=p+1;
    end
end



end