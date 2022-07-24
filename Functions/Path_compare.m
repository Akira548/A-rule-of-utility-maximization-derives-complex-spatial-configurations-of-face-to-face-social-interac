function [Rsquared_of_Path,p,F,frechet_dis] = Path_compare(RealPath,ModelPath,StartPoint,draw)

if nargin<4
    draw = 0;
end
    
    
poly1fun=fit(RealPath(:,1),RealPath(:,2),'poly1');
flip_angle = 45-atand(poly1fun.p1);

% Rotates all data to 45 degree for a better fitting

% Rotates around a fixed point
% x1=cos(angle)*x-sin(angle)*y;
% y1=cos(angle)*y+sin(angle)*x;
fliped_mean_Trace = [((RealPath(:,1)-StartPoint(1))*cosd(flip_angle)-(RealPath(:,2)-StartPoint(2))*sind(flip_angle))+StartPoint(1),...
    (RealPath(:,2)-StartPoint(2))*cosd(flip_angle)+(RealPath(:,1)-StartPoint(1))*sind(flip_angle)+StartPoint(2)];
fliped_ModelPath = [(ModelPath(:,1)-StartPoint(1))*cosd(flip_angle)-(ModelPath(:,2)-StartPoint(2))*sind(flip_angle)+StartPoint(1),...
    (ModelPath(:,2)-StartPoint(2))*cosd(flip_angle)+(ModelPath(:,1)-StartPoint(1))*sind(flip_angle)+StartPoint(2)];
f = fit(fliped_mean_Trace(:,1),fliped_mean_Trace(:,2),'poly2');
% f = fit(fliped_mean_Trace(:,1),fliped_mean_Trace(:,2),'linearinterp');
% f = fit(fliped_mean_Trace(:,1),fliped_mean_Trace(:,2),'smoothingspline');

[~,~,~,~,stats]=regress(fliped_ModelPath(:,2),[ones(size(fliped_ModelPath(:,2))) f(fliped_ModelPath(:,1))]);
Rsquared_of_Path= double(stats(1));
p=double(stats(3));
F=double(stats(2));
frechet_dis =  frechet(RealPath(:,1),RealPath(:,2),ModelPath(:,1),ModelPath(:,2));

if draw
figure
subplot(1,2,1)
plot(f,fliped_mean_Trace(:,1),fliped_mean_Trace(:,2),'.');hold on
plot(fliped_ModelPath(:,1),fliped_ModelPath(:,2),'.');
% plot(f((min(fliped_mean_Trace)-.3):0.1:(max(fliped_mean_Trace(:,2)+.3))),(min(mean_Trace{i_ori}(:,2))-.3):0.1:(max(mean_Trace{i_ori}(:,2)+.3)),'-','color',[1,0,0]);hold on
title(['Sequences: ' num2str(double(Rsquared_of_Path))]);
legend('off');
end

end

