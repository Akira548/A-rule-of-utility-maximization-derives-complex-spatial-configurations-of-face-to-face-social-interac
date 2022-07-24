function [Contrast,Best_ori] = CalculatePotentialContrast(DrawArea,Dis_step,Ori_step,VH_pos,draw,mode)
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
        r = 0.6; sigma1 = .2;sigma2 = 1.5;
        o_center = [VH_pos(:,1)+r*cos(VH_pos(:,3)), VH_pos(:,2)+r*sin(VH_pos(:,3))];
        Contrast{2}(x_index,y_index) = geomean(exp(-((ix-o_center(:,1)).^2+(iy-o_center(:,2)).^2)./(2*sigma2^2))-exp(-((ix-o_center(:,1)).^2+(iy-o_center(:,2)).^2)./(2*sigma1^2)));

%%  CPM concept proximec model  
        temp = exp(-(Rho-1.2).^2./2*(1)^2);
        Contrast{3}(x_index,y_index) = geomean(temp);%PM
%%
        r = 0.6; sigma = .2;
        o_center = [VH_pos(:,1)+r*cos(VH_pos(:,3)), VH_pos(:,2)+r*sin(VH_pos(:,3))];
        Contrast{4}(x_index,y_index) = geomean(exp(-((ix-o_center(:,1)).^2+(iy-o_center(:,2)).^2)./(2*sigma^2)));

%%         
        Orientation1_test = -pi:Ori_step:pi;
        for pi_i = 1:length(Orientation1_test)
            Theta1_test = wrapToPi(((Vector_Orientation-pi) - Orientation1_test(pi_i)));
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
        Contrast{1}(x_index,y_index) = P2(x_index,y_index);%SPM
        Orientation1_test = -pi:Ori_step:pi;
        for pi_i = 1:length(Orientation1_test)
            o_center2 = [ix+r*cos(Orientation1_test(pi_i)), iy+r*sin(Orientation1_test(pi_i))];
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

end

