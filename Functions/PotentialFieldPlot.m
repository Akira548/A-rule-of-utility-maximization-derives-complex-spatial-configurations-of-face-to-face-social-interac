function SaliencyMap = PotentialFieldPlot(DrawArea,Dis_step,Ori_step,VH_pos,mode,draw)
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
if nargin < 5
    mode = 'Psoc';draw=0;
end
%%%
if nargin < 6
    draw=0;
end

for ix = xmin:Dis_step:xmax
    for iy = ymin:Dis_step:ymax
        x_index = round((ix-xmin)/Dis_step+1);y_index = round((iy-ymin)/Dis_step+1);
        %%% theta2
        Orientation2 = wrapToPi(VH_pos(:,3));
        v = [(repmat(ix,size(VH_pos,1),1)-VH_pos(:,1)),(repmat(iy,size(VH_pos,1),1)-VH_pos(:,2))];
        [Vector_Orientation,Rho] = cart2pol(v(:,1),v(:,2));
        Theta2 = Orientation2 - Vector_Orientation;
        switch mode
            case 'Patt'
                SaliencyMap(x_index,y_index) = fun(para(1,:),[wrapToPi(Theta2),Rho]);
            case 'Prep'
                SaliencyMap(x_index,y_index) = fun(para(2,:),[wrapToPi(Theta2),Rho]);
            case 'Psoc'
                SaliencyMap(x_index,y_index) = fun(para(1,:),[wrapToPi(Theta2),Rho]).*(1-fun(para(2,:),[wrapToPi(Theta2),Rho]));
            case 'Pdyadic'
                SaliencyMap(x_index,y_index) = geomean(fun(para(1,:),[wrapToPi(Theta2),Rho].*(1-fun(para(2,:),[wrapToPi(Theta2),Rho]))));
        end                
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

