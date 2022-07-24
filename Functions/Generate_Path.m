function [ModelPath] = Generate_Path(PotentialMap,DrawArea,Dis_step,StartPoint,Step,Mode)

if nargin < 6
    Mode = 'normal';
    ratio = [];
end
% if nargin < 7
%     ratio = 2.2;
% end

%%% Area to Analysis
xmin = DrawArea(1,1);
xmax = DrawArea(1,2);
ymin = DrawArea(2,1);
ymax = DrawArea(2,2);

[dx,dy] = gradient(PotentialMap'./max(PotentialMap',[],'all'));
norm_dx= dx./sqrt(dx.^2+dy.^2);
norm_dy = dy./sqrt(dx.^2+dy.^2);
% norm_dx= dx;
% norm_dy = dy;

[peak_x,peak_y] = find(imregionalmax( PotentialMap ));
regionalmax_pos = [((peak_x-1)*Dis_step+xmin), ((peak_y-1)*Dis_step+ymin)];

Agent_pos = StartPoint;
ModelPath = StartPoint;

% if (strcmp(Mode, 'Combination'))
%     [max_x, max_y] = find(max( PotentialMap ));
%     max_pos = [((max_x-1)*Dis_step+xmin), ((max_y-1)*Dis_step+ymin)];
% end

% contour(xmin:Dis_step:xmax,ymin:Dis_step:ymax,PotentialMap');hold on
% scatter(regionalmax_pos(:,1),regionalmax_pos(:,2),'.k');hold on
% h=streamslice(xmin:Dis_step:xmax,ymin:Dis_step:ymax,[],dx,dy,[],0,150,0,1,'nearest');hold on;

for turn = 1:100
    %     scatter(Agent_pos(1),Agent_pos(2),'.r');hold on
    %     xlim([xmin xmax]);ylim([ymin ymax]);
    if any(sqrt((regionalmax_pos(:,1)-Agent_pos(1)).^2+(regionalmax_pos(:,2)-Agent_pos(2)).^2) < Step)
        break;
    end
    
    Agent_pos_index = [round((Agent_pos(2)-ymin)/Dis_step+1),round((Agent_pos(1)-xmin)/Dis_step+1)];
    %%%
    if (strcmp(Mode, 'Combination'))
        Relative_Distance = (regionalmax_pos-Agent_pos);
        [~,nearest]=min(sqrt(Relative_Distance(:,1).^2 + Relative_Distance(:,2).^2));
        PotentialMap_ = PotentialMap'./max(PotentialMap',[],'all');
        %         norm_target = 0.1*PotentialMap_(Agent_pos_index(1),Agent_pos_index(2)).*(Relative_Distance(nearest,:)/norm(Relative_Distance(nearest,:)));
        norm_target = PotentialMap_(round((StartPoint(2)-ymin)/Dis_step+1),round((StartPoint(1)-xmin)/Dis_step+1))...
            .*(Relative_Distance(nearest,:)/norm(Relative_Distance(nearest,:)));
        Associate_Vector = norm_target*ratio + [dx(Agent_pos_index(1),Agent_pos_index(2)),dy(Agent_pos_index(1),Agent_pos_index(2))]./norm([dx(Agent_pos_index(1),Agent_pos_index(2)),dy(Agent_pos_index(1),Agent_pos_index(2))]);
        Agent_pos = Agent_pos + Step.*(Associate_Vector./norm(Associate_Vector));
    else
        Agent_pos = Agent_pos + Step.*[dx(Agent_pos_index(1),Agent_pos_index(2)),dy(Agent_pos_index(1),Agent_pos_index(2))]./norm([dx(Agent_pos_index(1),Agent_pos_index(2)),dy(Agent_pos_index(1),Agent_pos_index(2))]);
    end
    
    %%%
    
    ModelPath = [ModelPath;Agent_pos];
end


end

