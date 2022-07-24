function [SaliencyMap,Best_Ori,Theta1_Intensity] = CalculatePotentialField(DrawArea,Dis_step,Ori_step,VH_pos,draw,mode)
%%% Formula
para = [1.523, 0.904, 4.173;
        0.042, 1.548, 0.247];
fun = @(beta,x) 1-exp(-( beta(3)./  ( (exp(beta(1).*abs(x(:,1)))).*x(:,2) )).^abs(beta(2)));
 
%%% Area to Analysis
xmin = DrawArea(1,1);
xmax = DrawArea(1,2);
ymin = DrawArea(2,1);
ymax = DrawArea(2,2);
[mesh_x,mesh_y]=meshgrid(xmin:Dis_step:xmax,ymin:Dis_step:ymax);

%%%
if nargin < 6
    mode = 0;
end
Theta1_Intensity = [];

%%%
for ix = xmin:Dis_step:xmax
    for iy = ymin:Dis_step:ymax
        x_index = round((ix-xmin)/Dis_step+1);y_index = round((iy-ymin)/Dis_step+1);
        %%% theta2
        Orientation2 = wrapToPi(VH_pos(:,3));
        v = [(repmat(ix,size(VH_pos,1),1)-VH_pos(:,1)),(repmat(iy,size(VH_pos,1),1)-VH_pos(:,2))];
        [Vector_Orientation,Rho] = cart2pol(v(:,1),v(:,2));
        Theta2 = Orientation2 - Vector_Orientation;
%         P1(x_index,y_index) = prod(fun(para(1,:),[wrapToPi(Theta2),Rho]).*(1-fun(para(2,:),[wrapToPi(Theta2),Rho])));
        
        %%% theta1
        Orientation1_test = -pi:Ori_step:pi;
        for pi_i = 1:length(Orientation1_test)
            Theta1_test = wrapToPi(((Vector_Orientation-pi) - Orientation1_test(pi_i)));
%             P2_all_angles(pi_i) = prod(fun(para(1,:),[Theta1_test,Rho]).*(1-fun(para(2,:),[Theta1_test,Rho]))).*P1(x_index,y_index);
            P2_all_angles(pi_i) = geomean((((fun(para(1,:),[wrapToPi(Theta2),Rho]).*(1-fun(para(2,:),[wrapToPi(Theta2),Rho]))).*(fun(para(1,:),[Theta1_test,Rho]).*(1-fun(para(2,:),[Theta1_test,Rho])))).^.5));
            if mode ~= 0
                Theta1_Intensity(x_index,y_index,pi_i) = P2_all_angles(pi_i);
            end
        end
        max_Ori1_index = find(P2_all_angles  == max(P2_all_angles));
        if ~isempty(max_Ori1_index)
            random_max_Ori1_index = max_Ori1_index(randi([1 length(max_Ori1_index)]));%random_max_Ori1_index = max_Ori1_index(1);
            Theta1(x_index,y_index,:) = wrapToPi(((Vector_Orientation-pi) - Orientation1_test(random_max_Ori1_index)));
            Best_Ori(x_index,y_index) = Orientation1_test(random_max_Ori1_index);
            P2(x_index,y_index) = P2_all_angles(random_max_Ori1_index);
        else
            P2(x_index,y_index) = nan;
        end
%         SaliencyMap(x_index,y_index) = P1(x_index,y_index).*P2(x_index,y_index);
        SaliencyMap(x_index,y_index) = P2(x_index,y_index);
    end
end

if draw
    figure
    contourf(mesh_x,mesh_y,(SaliencyMap'./max(SaliencyMap,[],'all')));hold on;
    [dx,dy]=gradient((SaliencyMap'./max(SaliencyMap,[],'all')));
    h=streamslice(mesh_x,mesh_y,[],dx,dy,[],0,150,0,5,'nearest');hold on;axis equal;
    for i = 1:length(h)
        h(i).Color = [0.1,0.1,0.1];
    end
    title({'Merged Potential'});xlabel('Distance(m)');ylabel('Distance(m)');
    axis equal;colorbar;
    shading interp
end

end

