function [deviation] = deviation_of_position(x,y,SaliencyMap,DrawArea,Dis_step)
%%% Area to Analysis
xmin = DrawArea(1,1);
xmax = DrawArea(1,2);
ymin = DrawArea(2,1);
ymax = DrawArea(2,2);
[peak_x,peak_y,~] = find(imregionalmax( SaliencyMap ));
Tran_peak_x = (peak_x-1)*Dis_step+xmin;    Tran_peak_y = (peak_y-1)*Dis_step+ymin;

for i = 1:size(x,1)
    for j = 1:size(x,2)
        deviation(i,j) = min(sqrt((x(i,j) - Tran_peak_x).^2+ (y(i,j) - Tran_peak_y).^2));
    end
end
end

