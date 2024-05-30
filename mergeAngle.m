function [partitionNode2,partitionDegree] = mergeAngle(partitionNode,theta,LowerTh,UpperTH,Degree, PartitionAngle)
s=1;
flag=0;
mergeIndex=0;
MergeWith=0;
partitionNode1 = partitionNode;
for i=1:size(partitionNode,1)
    %%%% count number of anchor points in eaach partition
  [NumberAnchorPoint] = CalculateNumberOfAnchorPoints(partitionNode);
    %%% Case 1 Assign a MV in each partition
    if (NumberAnchorPoint(i,1) >= LowerTh && NumberAnchorPoint(i,1) <= UpperTH )
        partitionNode(i,1)= partitionNode(i,1);
    end
    %%%% merge the partition
    if (i<length(NumberAnchorPoint) && i>1)
        if (NumberAnchorPoint(i,1) < LowerTh)
            if (NumberAnchorPoint(i-1,1) < NumberAnchorPoint(i+1,1))
                k=1;
                for j=1:size(partitionNode,2)
                    if (partitionNode(i-1,j) == 0)
                        partitionNode(i-1,j) = partitionNode(i,k);
                        partitionNode(i,k) = 0;
                        mergeIndex(s) = i;
                        MergeWith(s) = i-1;
                        k=k+1;
                    elseif (partitionNode(i-1,j) ~= 0 && j == size(partitionNode,2))
                        partitionNode(i-1,j) = partitionNode(i,k);
                        partitionNode(i,k) = 0;
                        mergeIndex(s) = i;
                        MergeWith(s) = i-1;
                        k=k+1;
                    end
                end
                k=1;
                s=s+1;
            else   %if (NumberAnchorPoint(i-1,1) > NumberAnchorPoint(i+1,1))
                  k=1;
                for j=1:size(partitionNode,2)
                     if (partitionNode(i+1,j) == 0)
                       partitionNode(i+1,j) = partitionNode(i,k);
                       partitionNode(i,k) = 0;
                       mergeIndex(s) = i;
                       MergeWith(s) = i+1;
                       k=k+1;
                     elseif (partitionNode(i+1,j) ~= 0 && j == size(partitionNode,2))
                       partitionNode(i+1,j) = partitionNode(i,k);
                       partitionNode(i,k) = 0;
                       mergeIndex(s) = i;
                       MergeWith(s) = i+1;
                       k=k+1;
                     end
                 end
                
                k=1;
               s=s+1;
            end
           % s=s+1;
            % RandomAngle=randi([Degree(1,i),Degree(1,i+1)],1,1);
        end
    elseif (i==length(NumberAnchorPoint))
        if (NumberAnchorPoint(i,1) < LowerTh)
            if (NumberAnchorPoint(1,1) < NumberAnchorPoint(i-1,1) && flag == 0)
                k=1;
                for j=1:size(partitionNode,2)
                    if (partitionNode(1,j)==0)
                        partitionNode(1,j)=partitionNode(i,k);
                        partitionNode(1,k)=0;
                        mergeIndex(s)=i;
                        MergeWith(s) = 1;
                        k=k+1;
                    elseif (partitionNode(1,j) ~=0 && j == size(partitionNode,2) )
                         partitionNode(1,j)=partitionNode(i,k);
                        partitionNode(1,k)=0;
                        mergeIndex(s)=i;
                        MergeWith(s) = 1;
                        k=k+1;
                    end
                end
                k=1;
                s=s+1;
            else
                  k=1;
                for j=1:size(partitionNode,2)
                    if (partitionNode(i-1,j)==0)
                        partitionNode(i-1,j)=partitionNode(i,k);
                        partitionNode(i,k)=0;
                        mergeIndex(s)=i;
                        MergeWith(s) = i-1;
                        k=k+1;
                    elseif (partitionNode(1,j) ~= 0 && j == size(partitionNode,2))
                        partitionNode(1,j+1) = partitionNode(i,k);
                        partitionNode(i,k)=0;
                        mergeIndex(s)=i;
                        MergeWith(s) = 1;
                        k=k+1;
                        
                    end
                end
                k=1;
                s=s+1;
                
            end
            % RandomAngle=randi([Degree(1,i),Degree(1,i+1)],1,1);
        end
    elseif (i==1 )
        if (NumberAnchorPoint(1,1) < LowerTh)
               flag=1;
            if (NumberAnchorPoint(length(NumberAnchorPoint),1) < NumberAnchorPoint(i+1,1))
                k=1;
                for j=1:size(partitionNode,2)
                    if (partitionNode(length(NumberAnchorPoint),j)==0)
                        partitionNode(length(NumberAnchorPoint),j)=partitionNode(1,k);
                        partitionNode(1,k)=0;
                        mergeIndex(s)=i;
                        MergeWith(s) = length(NumberAnchorPoint);
                        %partitionNode(1,j)=0
                        k=k+1;
                    elseif (partitionNode(length(NumberAnchorPoint),j) ~=0 && j == size(partitionNode,2))
                        partitionNode(length(NumberAnchorPoint),j)=partitionNode(1,k);
                        partitionNode(1,k)=0;
                        mergeIndex(s)=i;
                        MergeWith(s) = length(NumberAnchorPoint);
                        k=k+1;
                    end
                end
                k=1;
                s=s+1;
            else
                  k=1;
                for j=1:size(partitionNode,2)
                    if (partitionNode(i+1,j) == 0)
                        partitionNode(i+1,j) = partitionNode(i,k);
                        partitionNode(i,k)=0;
                        mergeIndex(s)=i;
                        MergeWith(s) = i+1;
                        k=k+1;
                        
                    elseif (partitionNode(i+1,j) ~= 0 && j == size(partitionNode,2))
                        partitionNode(i+1,j+1) = partitionNode(i,k);
                        partitionNode(i,k)=0;
                        mergeIndex(s)=i;
                        MergeWith(s) = i+1;
                        k=k+1;
                        
                    end
                end
                k=1;
                s=s+1;
            end
            % RandomAngle=randi([Degree(1,i),Degree(1,i+1)],1,1);
        end
    end
end
%%%%%% delate the merge angle
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