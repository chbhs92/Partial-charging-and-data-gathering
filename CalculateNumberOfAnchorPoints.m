function [NumberAnchorPoint] = CalculateNumberOfAnchorPoints ( partitionNode )
count=0;
for i=1:size(partitionNode,1)
    for j=1:size(partitionNode,2)
        if (partitionNode(i,j) ~=0 )
            count=count+1;
            NumberAnchorPoint(i,1) = count;
        end
    end
    count=0;
end