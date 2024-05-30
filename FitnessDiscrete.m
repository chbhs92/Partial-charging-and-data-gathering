function [totalDist] = FitnessDiscrete( pop,dmat, dmatBStoSensor)

%%%%%%%%%%%%%  Calcualte distance from base station to the anchor points
 
 % [st
  
     d= 0;
     f = size(pop,2);
d1=dmatBStoSensor(1,pop(1,1))+dmatBStoSensor(1,pop(1,f)); % from base station to first sensor nodes and last sensor nodes to base station
for i=2:size(pop,2)
    d=d+dmat(pop(1,i),pop(1,i-1));
end
totalDist=d+d1;

end