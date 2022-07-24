function [Shuffled_dis,Shuffled_Angle,dis_CI,Angle_CI] = Shuffled_level(Utility,pos,ori,DrawArea,Dis_step,n)
xmin = DrawArea(1,1);
ymin = DrawArea(2,1);

if nargin < 8
    n = 1000;
end

for tryi = 1:n
    rand_u_i = randi([1 length(Utility)]);
    U = Utility{rand_u_i};
    
    rand_p_i = randi([1 size(pos,1)]);
    P = pos(rand_p_i,:);
    
    O = ori(rand_p_i);
    
    [peak_x,peak_y,~] = find(imregionalmax( U{1} ));
    Tran_peak_x = (peak_x-1)*Dis_step+xmin;    Tran_peak_y = (peak_y-1)*Dis_step+ymin;
    Shuffled_dis(tryi) = min(sqrt((P(1) - Tran_peak_x).^2+ (P(2) - Tran_peak_y).^2));
    
    tempV = [cos(U{5}(round((P(1) - xmin)./Dis_step+1),round((P(2) - ymin)./Dis_step+1))),sin(U{5}(round((P(1) - xmin)./Dis_step+1),round((P(2) - ymin)./Dis_step+1)))];
    Shuffled_Angle(tryi) = rad2deg(wrapToPi(acos(dot([cos(O) sin(O)],tempV)/(norm([cos(O) sin(O)])*norm(tempV)))));
end
    [~, dis_CI]=MD_CI(Shuffled_dis ,0);
    [~, Angle_CI]=MD_CI(Shuffled_Angle ,0);
end
