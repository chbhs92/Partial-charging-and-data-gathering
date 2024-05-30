function [partitionNode1] = partitionGotZero(partitionNode)
s=1;
for i= 1 : size(partitionNode,1)
    %for j = 1:size(partitionNode,2)
        if (partitionNode(i,1) ~= 0)
            partitionNode1(s,:) = partitionNode(i,:);
            s=s+1;
        end
   % end
    
end

end