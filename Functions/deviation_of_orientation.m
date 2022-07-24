function [deviation] = deviation_of_orientation(x,y,vx,vy,SaliencyMap,DrawArea,Dis_step)
%%% Area to Analysis
xmin = DrawArea(1,1);
xmax = DrawArea(1,2);
ymin = DrawArea(2,1);
ymax = DrawArea(2,2);

for i = 1:size(x,1)
    for j = 1:size(x,2)
        tempV = [cos(SaliencyMap(round((x(i,j) - xmin)./Dis_step+1),round((y(i,j) - ymin)./Dis_step+1))),sin(SaliencyMap(round((x(i,j) - xmin)./Dis_step+1),round((y(i,j) - ymin)./Dis_step+1)))];
        deviation(i,j) = wrapToPi(acos(dot([vx(i,j) vy(i,j)],tempV)/(norm([vx(i,j) vy(i,j)])*norm(tempV))));
    end
end
end

