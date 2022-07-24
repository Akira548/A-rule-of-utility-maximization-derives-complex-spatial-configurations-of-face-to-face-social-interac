function result = ExpContrast_o_Distance(SaliencyMap,pos_I,DrawArea,Dis_step,VH_pos,isdraw)
%%% Area to Analysis
xmin = DrawArea(1,1);
xmax = DrawArea(1,2);
ymin = DrawArea(2,1);
ymax = DrawArea(2,2);
if nargin < 6
    isdraw = 0;
end
c2 = parula(10);
for modi = 1:size(SaliencyMap,2)
    tempdindex = pos_I;
    if modi == 1
        [peak_x,peak_y,~] = find(imregionalmax( SaliencyMap{modi} ));
        Tran_peak_x = (peak_x-1)*Dis_step+xmin;    Tran_peak_y = (peak_y-1)*Dis_step+ymin;
        O_dis(modi) = min(sqrt((tempdindex(1) - Tran_peak_x).^2+ (tempdindex(2) - Tran_peak_y).^2));
        
    else
        if modi == 2
%             if size(VH_pos,1) == 1
%                 r = 0.6; sigma1 = .2;sigma2 = 1.5;
%                 o_center = [VH_pos(:,1)+r*cos(VH_pos(:,3)), VH_pos(:,2)+r*sin(VH_pos(:,3))];
%                 testr = 0:0.01:3;
%                 testu= exp(-((testr-r).^2)./(2*sigma2^2))-exp(-((testr-r).^2)./(2*sigma1^2));
%                 maxu_r=(max(testr(testu==max(testu)))-min(testr(testu==max(testu))))/2;
%                 O_dis(modi) = abs(sqrt((tempdindex(1)-o_center(1))^2 + (tempdindex(2)-o_center(2))^2)-maxu_r);
%             else
%                 O_dis(modi) = min(sqrt((tempdindex(1) - Tran_peak_x).^2 + (tempdindex(2) - Tran_peak_y).^2));
%             end
            q=Estimate_approaching_goal_pose(VH_pos);
            if ~isempty(q)
                O_dis(modi) = min(sqrt((tempdindex(1) - q(:,1)).^2 + (tempdindex(2) - q(:,2)).^2));
            else
                O_dis(modi) = nan;
            end
        else
            [peak_x,peak_y] = find(SaliencyMap{modi} == max(SaliencyMap{modi},[],'all'));
            Tran_peak_x = (peak_x-1)*Dis_step+xmin;    Tran_peak_y = (peak_y-1)*Dis_step+ymin;
            O_dis(modi) = min(sqrt((tempdindex(1) - Tran_peak_x).^2 + (tempdindex(2) - Tran_peak_y).^2));
        end
    end
end
if isdraw
[~,~,CI(1,:)]=ttest(O_dis(:,1));
[~,~,CI(2,:)]=ttest(O_dis(:,3));
[~,~,CI(3,:)]=ttest(O_dis(:,2));
figure
bar(1:size(O_dis(:,[1 3 2]),2),nanmean(O_dis(:,[1 3 2]),1),'FaceColor',[.5 .5 .5]);hold on
her = errorbar(1:size(O_dis(:,[1 3 2]),2),nanmean(O_dis(:,[1 3 2]),1),abs(CI(:,1)'-nanmean(O_dis(:,[1 3 2]),1))) ;
her.LineStyle = 'none';her.Color = 'k';
xticklabels({'OurModel','CPM','VAM'});gca.FontSize=12;
ylabel('Deviation Distance(m)','fontsize',12);grid off;box off;
set(gcf,'position',[0,0,700,400]);
clc
disp('done!');

MD_CI(O_dis(:,1));
MD_CI(O_dis(:,3));
MD_CI(O_dis(:,2));

end

result = O_dis(:,[1 3 2]);
end

