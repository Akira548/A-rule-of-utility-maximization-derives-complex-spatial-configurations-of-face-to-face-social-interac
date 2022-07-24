function [Contrast,Best_ori] = CalculatePotentialFieldContrast(DrawArea,Dis_step,Ori_step,VH_pos,draw,mode)
%%% Formula
para = [-1.523, 0.904, 4.173;
        -0.042, 1.548, 0.247];
fun = @(beta,x) 1-exp(-( beta(3).*(exp(beta(1).*abs(x(:,1))))./x(:,2)).^beta(2));

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
Best_ori = [];

%%%
for ix = xmin:Dis_step:xmax
    for iy = ymin:Dis_step:ymax
        x_index = round((ix-xmin)/Dis_step+1);y_index = round((iy-ymin)/Dis_step+1);
        Orientation2 = wrapToPi(VH_pos(:,3));
        v = [(repmat(ix,size(VH_pos,1),1)-VH_pos(:,1)),(repmat(iy,size(VH_pos,1),1)-VH_pos(:,2))];
        [Vector_Orientation,Rho] = cart2pol(v(:,1),v(:,2));
        Theta2 = wrapToPi(Orientation2 - Vector_Orientation);
%%  HVFF modified
        % r = 0.6m. sigma x,y^2 = 4m
        r = 0.6; sigma1 = .2;sigma2 = 1.5;
        o_center = [VH_pos(:,1)+r*cos(VH_pos(:,3)), VH_pos(:,2)+r*sin(VH_pos(:,3))];
        Contrast{2}(x_index,y_index) = geomean(exp(-((ix-o_center(:,1)).^2+(iy-o_center(:,2)).^2)./(2*sigma2^2))-exp(-((ix-o_center(:,1)).^2+(iy-o_center(:,2)).^2)./(2*sigma1^2)));


%         Contrast{2}(x_index,y_index) = sum(-20./Rho.^3 + 15./Rho.^0.8);%SPF
%%  CPM concept proximec model  
%         Contrast{3}(x_index,y_index) = prod(1.214*exp(-0.312*Rho)-0.096);%PM
%         temp =exp(-0.9228*Rho);
        temp = exp(-(Rho-1.2).^2./2*(1)^2);
%         temp =1.214*exp(-0.312*Rho)-0.096;
%         temp(temp<0)=0;
        Contrast{3}(x_index,y_index) = geomean(temp);%PM
%%
%         temp = [];
%         for Ri=1:length(Rho)
%             if (Theta2 > pi/3 | Theta2 < -pi/3)
%                 temp(Ri) = 0;%VAM
%             else
%                 temp(Ri) = (exp(-Rho(Ri)/0.5));%VAM
%             end
%         end
%         Contrast{4}(x_index,y_index) = geomean(temp);
% %%        
%         Aij = [];
%         for Ri=1:length(Rho)
%                 Aij(Ri) = exp(-Rho(Ri)^2./0.5-(Theta2).^2./0.5);%VAM
%         end
%         Contrast{7}(x_index,y_index) = geomean(temp);
%         r = 0.6; sigma1 = .6;sigma2 = 1.5;
%         o_center = [VH_pos(:,1)+r*cos(VH_pos(:,3)), VH_pos(:,2)+r*sin(VH_pos(:,3))];
%         Contrast{4}(x_index,y_index) = geomean(exp(-((ix-o_center(:,1)).^2+(iy-o_center(:,2)).^2)./(2*sigma2^2))-exp(-((ix-o_center(:,1)).^2+(iy-o_center(:,2)).^2)./(2*sigma1^2)));
        r = 0.6; sigma = .2;
        o_center = [VH_pos(:,1)+r*cos(VH_pos(:,3)), VH_pos(:,2)+r*sin(VH_pos(:,3))];
        Contrast{4}(x_index,y_index) = geomean(exp(-((ix-o_center(:,1)).^2+(iy-o_center(:,2)).^2)./(2*sigma^2)));

%%         
%         P1(x_index,y_index) = prod(fun(para(1,:),[wrapToPi(Theta2),Rho]).*(1-fun(para(2,:),[wrapToPi(Theta2),Rho])));
        %%% theta2 for our model
        Orientation1_test = -pi:Ori_step:pi;
        for pi_i = 1:length(Orientation1_test)
            Theta1_test = wrapToPi(((Vector_Orientation-pi) - Orientation1_test(pi_i)));
%             P2_all_angles(pi_i) = prod(fun(para(1,:),[Theta1_test,Rho]).*(1-fun(para(2,:),[Theta1_test,Rho])));
            P2_all_angles(pi_i) = geomean((((fun(para(1,:),[wrapToPi(Theta2),Rho]).*(1-fun(para(2,:),[wrapToPi(Theta2),Rho]))).*(fun(para(1,:),[Theta1_test,Rho]).*(1-fun(para(2,:),[Theta1_test,Rho])))).^.5));
            if mode ~= 0
                Theta1_Intensity(x_index,y_index,pi_i) = P2_all_angles(pi_i);
            end
        end
        max_Ori1_index = find(P2_all_angles  == max(P2_all_angles));
        if ~isempty(max_Ori1_index)
            random_max_Ori1_index = max_Ori1_index(randi([1 length(max_Ori1_index)]));%random_max_Ori1_index = max_Ori1_index(1);
            Theta1(x_index,y_index,:) = wrapToPi(((Vector_Orientation-pi) - Orientation1_test(random_max_Ori1_index)));
            Contrast{5}(x_index,y_index) = Orientation1_test(random_max_Ori1_index);
            P2(x_index,y_index) = P2_all_angles(random_max_Ori1_index);
        else
            P2(x_index,y_index) = nan;
            Contrast{5}(x_index,y_index)=nan;
        end
        if mode ~= 0
            Best_ori(x_index,y_index) = Orientation1_test(random_max_Ori1_index);
        end
%         Contrast{1}(x_index,y_index) = P1(x_index,y_index).*P2(x_index,y_index);%SPM
        Contrast{1}(x_index,y_index) = P2(x_index,y_index);%SPM
        %%% theta for VAM
        Orientation1_test = -pi:Ori_step:pi;
        for pi_i = 1:length(Orientation1_test)
%             Theta1_test = wrapToPi(((Vector_Orientation-pi) - Orientation1_test(pi_i)));
            o_center2 = [ix+r*cos(Orientation1_test(pi_i)), iy+r*sin(Orientation1_test(pi_i))];
            %             for Ri=1:length(Rho)
%                 if (Theta1_test > pi*1/3 | Theta1_test < -pi*1/3)
%                     temp(Ri) = 0;%P
%                 else
%                     temp(Ri) = (exp(-Rho(Ri)/0.5));
%                 end
%             end
            PM_test2(pi_i) = geomean(exp(-((VH_pos(:,1)-o_center2(:,1)).^2+(VH_pos(:,2)-o_center2(:,2)).^2)./(2*sigma2^2))-exp(-((VH_pos(:,1)-o_center2(:,1)).^2+(VH_pos(:,2)-o_center2(:,2)).^2)./(2*sigma1^2)));
        end
        
        max_Ori1_index2 = find(PM_test2  == max(PM_test2));
        if ~isempty(max_Ori1_index2)
            random_max_Ori1_index2 = max_Ori1_index2(randi([1 length(max_Ori1_index2)]));%             random_max_Ori1_index = max_Ori1_index(1);
            Contrast{6}(x_index,y_index) = Orientation1_test(random_max_Ori1_index2);
        else
            Contrast{6}(x_index,y_index) = nan;
        end
        %%%%%%%%%%%%%%%%%%%%%%  
        
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

