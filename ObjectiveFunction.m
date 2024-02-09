function [Objective, ChargingSequence ]=ObjectiveFunction(Positions,bs,dim,N,LowerTh,UpperTH, MaxNumberMV)

% degree of partitioning 
for i=1:size(Positions,1)
    for j=1:size(Positions,2)
        Degree(i,j)=360*Positions(i,j); % angle of partitions
    end
end

%% find anchor point in each partition
% calculate theta of each anchor points
for i=1:length(bs)
    
    theta(i,1)= atan2(bs(i,2)-500,bs(i,1)-500)*(180/pi); % find angle of points
    
end
%%%%%%%%%%%%%%%%%%% convert negative angle into positive
for i=1:length(bs)
    if (theta(i,1)<0)
        theta(i,1)=theta(i,1)+360;
    end
end

[PartitionAngle, PartitionIndex] = sort(Degree); %%%%% sort the partition angle
k=1;
for i=1:length(Degree)
    if (i~=length(Degree))
        for j=1:length(bs)
            if(theta(j,1) > PartitionAngle(1,i) && theta(j,1) < PartitionAngle(1,i+1) )   %%%%%%%%%%%%%%%%%%%%%
                partitionNode(i,k)=j;                    %%%% find anchor points in each region
                k=k+1;
            end
        end
        k=1;
    end
    if (i==length(Degree))
        for j=1:length(bs)
            if(theta(j,1) > PartitionAngle(1,i) || theta(j,1) < PartitionAngle(1,1) )
                partitionNode(i,k)=j;
                k=k+1;
            end
        end
    end
end
%%%%%%%%%%%%%%% if some partition got zero anchor points  %%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
[partitionNode] = partitionGotZero(partitionNode);
 
%%%%%%%%%%%%%%%%%% case 1 condition for merging the partition  %%%%%%%%%%%%%%

[NumberAnchorPoint] = CalculateNumberOfAnchorPoints(partitionNode);
c=1;
while (1)
    if (find(  NumberAnchorPoint < LowerTh))
        [partitionNode1,PartitionAngle] = mergeAngle (partitionNode,theta,LowerTh,UpperTH,Degree, PartitionAngle);
        [NumberAnchorPoint] = CalculateNumberOfAnchorPoints(partitionNode1);
        partitionNode = partitionNode1;
        c=c+1;
    else
      break;
      
    end
end
if (c>1)
   partitionNode1 = partitionNode1;
else
    partitionNode1=partitionNode;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% case 2 partition the big area (partition) into two partition  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[NumberAnchorPoint1] = CalculateNumberOfAnchorPoints(partitionNode1);

while (1)
    if (find( NumberAnchorPoint1 > UpperTH ))
      [partitionNode3, PartitionAngle1] = PartitionRegion(partitionNode1, PartitionAngle, UpperTH, theta);
      [NumberAnchorPoint1] = CalculateNumberOfAnchorPoints(partitionNode3);
      partitionNode1 = partitionNode3;
    else
         break;
         
    end   
end
  
 
%%%%%%%%%%%%%%%% sub tour length %%%%%%%%%%%%%%%%%%%%%%%% 

for i=1:size(partitionNode1,1)
    
        [SubTourLength(1,i), Sequence] =   tsp_greedy ( bs, partitionNode1(i,:) );%TSP_GA1(bs, partitionNode1(i,:));
        ChargingSequence{i,:} = Sequence;
end
len2 = size(partitionNode1,1); 
MaxTour = max(SubTourLength);
MinTour = min(SubTourLength);
DiffTour = (MaxTour -  MinTour)/2;
Mean=mean(SubTourLength);
Std=std(SubTourLength);
Objective=Std/DiffTour + len2/MaxNumberMV; % + (max(SubTourLength)/Mean);  %% objective value
 NumberMV =  size(partitionNode1,1); 
partitionAnchor = partitionNode1;



