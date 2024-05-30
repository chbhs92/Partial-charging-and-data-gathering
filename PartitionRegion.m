function [partitionNode4, PartitionAngle1] = PartitionRegion(partitionNode1, PartitionAngle, UpperTH, theta)

[NumberAnchorPoint] = CalculateNumberOfAnchorPoints(partitionNode1);
 k=1;
for i=1:size(partitionNode1,1)
    if (NumberAnchorPoint(i,1) > UpperTH)
              for j=1:size(partitionNode1,2)
                  if (partitionNode1(i,j) ~= 0)
                      ElementAngle(k,j) = theta(partitionNode1(i,j));  % find region who has greater element than upper thresold
                      partitionIndex(k) = i;
                  end
              end
               k=k+1;
    end    
end
[NumberElement] = CalculateNumberOfAnchorPoints(ElementAngle);  %%%% find number of anchor points
midElement=ceil(NumberElement/2);  % middle element
partitionNode2=partitionNode1;    %%% replica 
s=1;
for i=1:size(ElementAngle,1)
    [partitionElement(i,:), elemntIndex(i,:)]=sort(ElementAngle(i,:));   % sort the anchor point ( big region ) with angle
     for j = 1:size(ElementAngle,2) 
         if ( partitionElement(i,j) ~= 0)   %%% removing zero items
             SortedElement(i,s) = partitionElement(i,j);   
             SortedIndex(i,s) = elemntIndex(i,j);
             s=s+1;
         end
         
     end
     s=1;  
end
%%%%%%%%%  enter partition anchor point into another region
len=size(partitionNode1,2);
r=size(partitionNode1,1);
p=1;
w=1;
partitionNode3=zeros(r,len);
for i=1:size(partitionNode1,1)
    if (find(i == partitionIndex))
        k=find(i == partitionIndex);
        for j=1:midElement(k)    %%% first half
            if (partitionNode1(i,j) ~= 0)
                partitionNode3(p,j) = partitionNode1(i,SortedIndex( k,j));
                
            end
        end
        for x = midElement(k)+1: len      %%% second half
            if (partitionNode1(i,x) ~= 0)
                partitionNode3(p+1,w) = partitionNode1(i,SortedIndex(k ,x));
                w=w+1;
            end
        end
        
        w=1;
        p = p+2;
    else
        partitionNode3(p,:) = partitionNode1(i,:);
        p = p+1;
    end
end
g=1;
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