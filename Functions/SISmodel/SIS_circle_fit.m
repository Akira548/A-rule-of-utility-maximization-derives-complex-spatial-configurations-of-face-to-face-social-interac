%%
function [r,xc,yc]=SIS_circle_fit(v)

if size(v,1) == 2

        %%% =2
        % close all;
%         theta = pi/3;
%         v=[-1,0,-theta;1 ,0,theta-pi];
        
        uv = [v(2,1)-v(1,1), v(2,2)-v(1,2)];
        [Vector_Orientation1,interpersonal_dis] = cart2pol(uv(:,1),uv(:,2));
        Theta1 = wrapToPi(v(1,3) - Vector_Orientation1);
        [Vector_Orientation2,~] = cart2pol(-uv(:,1),-uv(:,2));
        Theta2 = -wrapToPi(v(2,3) - Vector_Orientation2);
        
        r=(interpersonal_dis/2)./cos((Theta1+Theta2)/5);
        cir_move_ang = Vector_Orientation1 + (Theta1+Theta2)/5;
        xc = v(1,1)+r*cos(cir_move_ang);
        yc = v(1,2)+r*sin(cir_move_ang);
        
%         figure
%         quiver(v(:,1),v(:,2),cos(v(:,3)),sin(v(:,3)),.2);
%         rectangle('Curvature',[1 1],'Position',[xc-r yc-r 2*r 2*r]);hold on; % Plot circle fit
%         plot(v(:,1),v(:,2),'-k');hold on;
%         plot(xc,yc,'.r');hold on;
%         axis equal
else        
    if size(v,1) > 2
        %%% >=3
        % close all;
%         theta = pi/3;
%         v=[-1,0,theta;1 ,0,-theta-pi;2 ,1,-theta-pi];
        
        [r,xc,yc] = circfit(v(:,1),v(:,2));
%         
%         figure
%         quiver(v(:,1),v(:,2),cos(v(:,3)),sin(v(:,3)),.2);
%         rectangle('Curvature',[1 1],'Position',[xc-r yc-r 2*r 2*r]);hold on; % Plot circle fit
%         plot(v(:,1),v(:,2),'-k');hold on;
%         plot(xc,yc,'.r');hold on;
%         axis equal
    end
end
end