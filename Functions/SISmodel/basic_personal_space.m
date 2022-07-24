function f_p=basic_personal_space(DrawArea,Dis_step,VH_pos)
%%% Area or Point to Analysis
if size(DrawArea,1)>1
    xmin = DrawArea(1,1);
    xmax = DrawArea(1,2);
    ymin = DrawArea(2,1);
    ymax = DrawArea(2,2);
else
    xmin = DrawArea(1);
    xmax = DrawArea(1);
    ymin = DrawArea(2);
    ymax = DrawArea(2);
    Dis_step = 1;
end

sig_0_px = 0.45;
sig_0_py = 0.45;
Ap =1;

%%%
for x = xmin:Dis_step:xmax
    for y = ymin:Dis_step:ymax
        x_index = round((x-xmin)/Dis_step+1);y_index = round((y-ymin)/Dis_step+1);
        d = sqrt((x - VH_pos(:,1)).^2+(y - VH_pos(:,2)).^2);
        theta = atan2(y - VH_pos(:,2), x - VH_pos(:,1));
        f_p(x_index,y_index) = max(Ap*exp(-((d.*cos(theta-VH_pos(:,3))./(sqrt(2)*sig_0_px)).^2 +(d.*sin(theta-VH_pos(:,3))./(sqrt(2)*sig_0_py)).^2)));
    end
end

