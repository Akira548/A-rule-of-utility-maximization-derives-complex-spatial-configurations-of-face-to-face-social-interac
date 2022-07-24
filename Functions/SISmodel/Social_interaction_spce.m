function f_g=Social_interaction_spce(DrawArea,Dis_step,VH_pos)

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

[r,cx,cy]=SIS_circle_fit(VH_pos);

sig_gx = r/2;
sig_gy = r/2;
Ap =1;

%%%
for x = xmin:Dis_step:xmax
    for y = ymin:Dis_step:ymax
        x_index = round((x-xmin)/Dis_step+1);y_index = round((y-ymin)/Dis_step+1);
        d = sqrt((x - cx)^2+(y - cy)^2);
        theta = atan2(y - cy, x - cx);
        f_g(x_index,y_index) = max(Ap*exp(-((d*cos(theta)./(sqrt(2)*sig_gx)).^2 +(d*sin(theta)./(sqrt(2)*sig_gy)).^2)));
    end
end
