function [] = fitting_test(mode,fun,para,DrawArea,Dis_step,Ori_step,VH_pos,draw)

%%% Area to Analysis
xmin = DrawArea(1,1);
xmax = DrawArea(1,2);
ymin = DrawArea(2,1);
ymax = DrawArea(2,2);
[mesh_x,mesh_y]=meshgrid(xmin:Dis_step:xmax,ymin:Dis_step:ymax);
switch mode
    %%%
    case {'Att', 'Rep'}
        for ix = xmin:Dis_step:xmax
            for iy = ymin:Dis_step:ymax
                x_index = round((ix-xmin)/Dis_step+1);y_index = round((iy-ymin)/Dis_step+1);
                Orientation2 = wrapToPi(VH_pos(:,3));
                v = [(repmat(ix,size(VH_pos,1),1)-VH_pos(:,1)),(repmat(iy,size(VH_pos,1),1)-VH_pos(:,2))];
                [Vector_Orientation,Rho] = cart2pol(v(:,1),v(:,2));
                Theta2 = Orientation2 - Vector_Orientation;
                P = fun(para(1,:),[wrapToPi(Theta2),Rho]);
                SaliencyMap(x_index,y_index) = prod(P);
            end
        end
        
    case {'Field'}
        
        for ix = xmin:Dis_step:xmax
            for iy = ymin:Dis_step:ymax
                x_index = round((ix-xmin)/Dis_step+1);y_index = round((iy-ymin)/Dis_step+1);
                Orientation2 = wrapToPi(VH_pos(:,3));
                v = [(repmat(ix,size(VH_pos,1),1)-VH_pos(:,1)),(repmat(iy,size(VH_pos,1),1)-VH_pos(:,2))];
                [Vector_Orientation,Rho] = cart2pol(v(:,1),v(:,2));
                Theta2 = Orientation2 - Vector_Orientation;
                P(x_index,y_index) = fun(para(1,:),[wrapToPi(Theta2),Rho]).*(1-fun(para(2,:),[wrapToPi(Theta2),Rho]));
                
                %%% theta2
                Orientation1_test = -pi:Ori_step:pi;
                for pi_i = 1:length(Orientation1_test)
                    Theta1_test = wrapToPi(((Vector_Orientation-pi) - Orientation1_test(pi_i)));
                    P2_all(pi_i) = prod(fun(para(1,:),[Theta1_test,Rho]).*(1-fun(para(2,:),[Theta1_test,Rho])));
                end
                
                max_Ori1_index = find(P2_all  == max(P2_all));
                if ~isempty(max_Ori1_index)
                    random_max_Ori1_index = max_Ori1_index(randi([1 length(max_Ori1_index)]));%             random_max_Ori1_index = max_Ori1_index(1);
                    Theta1(x_index,y_index,:) = wrapToPi(((Vector_Orientation-pi) - Orientation1_test(random_max_Ori1_index)));
                    Orientation1(x_index,y_index) = Orientation1_test(random_max_Ori1_index);
                    P2(x_index,y_index) = P2_all(random_max_Ori1_index);
                else
                    P2(x_index,y_index) = nan;
                end
                
                SaliencyMap(x_index,y_index) = P(x_index,y_index).*P2(x_index,y_index);
            end
        end
        
    case {'V'}
        
        for ix = xmin:Dis_step:xmax
            for iy = ymin:Dis_step:ymax
                x_index = round((ix-xmin)/Dis_step+1);y_index = round((iy-ymin)/Dis_step+1);
                Orientation2 = wrapToPi(VH_pos(:,3));
                v = [(repmat(ix,size(VH_pos,1),1)-VH_pos(:,1)),(repmat(iy,size(VH_pos,1),1)-VH_pos(:,2))];
                [Vector_Orientation,Rho] = cart2pol(v(:,1),v(:,2));
                Theta2 = Orientation2 - Vector_Orientation;
                P(x_index,y_index) = fun(para(1,:),[wrapToPi(Theta2),Rho]).*(1-fun(para(2,:),[wrapToPi(Theta2),Rho]));
                SaliencyMap(x_index,y_index) = P(x_index,y_index);
            end
        end
        
    case {'AR'}
        for ix = xmin:Dis_step:xmax
            for iy = ymin:Dis_step:ymax
                x_index = round((ix-xmin)/Dis_step+1);y_index = round((iy-ymin)/Dis_step+1);
                Orientation2 = wrapToPi(VH_pos(:,3));
                v = [(repmat(ix,size(VH_pos,1),1)-VH_pos(:,1)),(repmat(iy,size(VH_pos,1),1)-VH_pos(:,2))];
                [Vector_Orientation,Rho] = cart2pol(v(:,1),v(:,2));
                Theta2 = Orientation2 - Vector_Orientation;
                P(x_index,y_index) = fun(para(1,:),[wrapToPi(Theta2),Rho]).*(1-fun(para(2,:),[wrapToPi(Theta2),Rho]));
                [dx,dy]=gradient(SaliencyMap'./max(SaliencyMap,[],'all'));

                SaliencyMap(x_index,y_index) = P(x_index,y_index);
            end
        end
end

if draw
    figure
%     contourf(mesh_x,mesh_y,(SaliencyMap'./max(SaliencyMap,[],'all')));hold on;
    surf(mesh_x,mesh_y,SaliencyMap');shading interp;hold on;
    view(0,90);
%     [dx,dy]=gradient((SaliencyMap'./max(SaliencyMap,[],'all')));
%     h=streamslice(mesh_x,mesh_y,[],dx,dy,[],0,150,0,1,'nearest');hold on;axis equal;
    %[max_x_index,max_y_index] = find(SaliencyMap == max(SaliencyMap,[],1:2));
%     [max_x_index,max_y_index,~] = find(imregionalmax((SaliencyMap./max(SaliencyMap,[],'all'))));
%     max_x = (max_x_index-1)*Dis_step+xmin;
%     max_y = (max_y_index-1)*Dis_step+ymin;
%     scatter(VH_pos(:,1),VH_pos(:,2),100,'.r');hold on;
%     scatter(max_x,max_y,150,'xr');hold on;
%     quiver(VH_pos(:,1) ,VH_pos(:,2),cos(VH_pos(:,3)),sin(VH_pos(:,3)) ,.3,'.r');hold on;
%     axis equal;ylim([ymin ymax]);xlim([xmin xmax]);box off;grid off;
%     for i = 1:length(h)
%         h(i).Color = [0.1,0.1,0.1];
%     end
    title({'Merged Potential'});xlabel('Distance(m)');ylabel('Distance(m)');
%     axis equal;
    
end

end

