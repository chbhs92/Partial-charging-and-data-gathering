function [total_distance,tour_path]= tsp_greedyDiscrete( dmat, partitionNode,  dmatBStoSensor )
 %%%%%%%%%%%%%%%%% distance matrix for each region
  [r,c] = size(dmat);
  distance = dmat;
  distance(r+1,:) = dmatBStoSensor(1,:);
  t = dmatBStoSensor(1,:)';
  t(r+1,1) = 0;
  distance (:,c+1) = t(:,1);
 
 
%
 
  [ m, n ] = size ( distance );
 
 
  cost_best = Inf;
  start = n ; % from base station
  %for start = 1 : n

    p = path_greedy ( n, distance, start );
    cost = path_cost ( n, distance, p );
    total_distance=cost;
    p(1) = [];   %% remove base station from charging schedule
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
 

 
