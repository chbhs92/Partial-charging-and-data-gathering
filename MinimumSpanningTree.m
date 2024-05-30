function [total_distance,tour_path]= MinimumSpanningTree ( bs1, partitionNode1 )
 %%%%%%%%%%%%%%%%% distance matrix for each region
 bs=bs1;
 BSX = 50;
BSY = 50;
partitionNode = partitionNode1;
l=0;
for i=1:length(partitionNode)
    if (partitionNode(i) ~= 0)
        l =l+1;
    end
end
 partitionNode(l+1) = length(bs)+1;
 bs(length(bs)+1,1) = BSX;
 bs(length(bs)+1,2) = BSY;
 for i=1:l+1
      for j=1:l+1
          distance(i,j) = sqrt((bs(partitionNode(i),1)-bs(partitionNode(j),1))^2 + (bs(partitionNode(j),1)-bs(partitionNode(j),1))^2);
      end
 end
 
 
%
 
  [ m, n ] = size ( distance );
 
 
  cost_best = Inf;
  start = n ; % from base station
  %for start = 1 : n

    p = path_greedy ( n, distance, start );
    cost = path_cost ( n, distance, p );
    total_distance=cost;
    p(1) = [];   %% remove base station to charging schedule
    tour_path = partitionNode(p(:,1));

 
%  end
end
   

function cost = path_cost ( n, distance, p )

  cost = 0.0;
  i1 = n;
  for i2 = 1 : n
    cost = cost + distance ( p(i1), p(i2) );
    i1 = i2;
  end

  return
end

function p = path_greedy ( n, distance, start )
 
%
  p = zeros ( n, 1 );
  p(1) = start;

  d = distance(1:n,1:n);
  d(:,start) = Inf;

  for i = 1 : n
    d(i,i) = Inf;
  end

  from = start;
  for j = 2 : n
    [ ~, to ] = min ( d(from,:) );
    p(j) = to;
    d(:,to) = Inf;
  end
  
 
end
 

 
