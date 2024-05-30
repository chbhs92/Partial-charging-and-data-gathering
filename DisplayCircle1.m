function h = DisplayCircle1(bs,r,sensor_x,sensor_y)
figure(2)
hold on
th = 0:pi/50:2*pi;
for i=1:length(bs)
xunit = r * cos(th) + bs(i,1);
yunit = r * sin(th) + bs(i,2);
h = plot(xunit, yunit);
end
plot(sensor_x,sensor_y,'r.')
plot(bs(:,1),bs(:,2),'g+')
hold off