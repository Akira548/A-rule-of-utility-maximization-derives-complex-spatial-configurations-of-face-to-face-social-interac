function result = Contrast_o_Distance(SaliencyMap,pos_I,DrawArea,Dis_step,VH_pos)
%%% Area to Analysis
xmin = DrawArea(1,1);
xmax = DrawArea(1,2);
ymin = DrawArea(2,1);
ymax = DrawArea(2,2);

c2 = parula(10);
for i = 1:size(SaliencyMap,2)
    for modi = 1:size(SaliencyMap{i},2)
            tempdindex = pos_I(:,:,i);
        if modi == 1
            [peak_x,peak_y,~] = find(imregionalmax( SaliencyMap{i}{modi} ));
            Tran_peak_x = (peak_x-1)*Dis_step+xmin;    Tran_peak_y = (peak_y-1)*Dis_step+ymin;
            O_dis(i,modi) = min(sqrt((tempdindex(1) - Tran_peak_x).^2+ (tempdindex(2) - Tran_peak_y).^2));
            
        else
            if modi == 2
                q=Estimate_approaching_goal_pose(VH_pos{i});
                if ~isempty(q)
                    O_dis(i,modi) = min(sqrt((tempdindex(1) - q(:,1)).^2 + (tempdindex(2) - q(:,2)).^2));
                else
                    O_dis(i,modi) = nan;
                end
            else
            
            [peak_x,peak_y] = find(SaliencyMap{i}{modi} == max(SaliencyMap{i}{modi},[],'all'));
            Tran_peak_x = (peak_x-1)*Dis_step+xmin;    Tran_peak_y = (peak_y-1)*Dis_step+ymin;
            O_dis(i,modi) = min(sqrt((tempdindex(1) - Tran_peak_x).^2 + (tempdindex(2) - Tran_peak_y).^2));
            end
        end
    end
end

% for tryi = 1:1000
%     randpermi = randperm(size(SaliencyMap,2));
%     while sum(randpermi == 1:size(SaliencyMap,2)) ~= 0 
%         randpermi = randperm(size(SaliencyMap,2));
%     end
%     for i = 1:size(SaliencyMap,2)
%         tempdindex = pos_I(:,:,i);
%         [peak_x,peak_y,~] = find(imregionalmax( SaliencyMap{randpermi(i)}{1} ));
%         Tran_peak_x = (peak_x-1)*Dis_step+xmin;    Tran_peak_y = (peak_y-1)*Dis_step+ymin;
%         O_dis_stength_I_randperm(i,tryi) = min(sqrt((tempdindex(1) - Tran_peak_x).^2+ (tempdindex(2) - Tran_peak_y).^2));
%     end
%     clc
%     disp(['processing:' num2str(round(tryi/1000*100)) '%']);
% end

% figure
% bar(1:size(SaliencyMap{i},2),nanmean(O_dis,1),'FaceColor',[.5 .5 .5]);hold on
% her = errorbar(1:size(SaliencyMap{i},2),nanmean(O_dis,1),std(O_dis)) ;
% her.LineStyle = 'none';her.Color = 'k';
% xticklabels({'Dataset','AR','AR2','SPM','PM','Visual'});gca.FontSize=12;
% ylabel('The Distance from Peak(m)','fontsize',12);grid off;box off;



% figure
% bar(1:size(O_dis(:,[1 2 3 4]),2),nanmean(O_dis(:,[1 2 3 4]),1),'FaceColor',[.5 .5 .5]);hold on
% % her = errorbar(1:size(O_dis(:,[1 2 3 4]),2),nanmean(O_dis(:,[1 2 3 4]),1),std(O_dis(:,[1 2 3 4]))./sqrt(size(O_dis(:,[1 2 3 4]),1))) ;
% her = errorbar(1:size(O_dis(:,[1 2 3 4]),2),nanmean(O_dis(:,[1 2 3 4]),1),abs(CI(:,1)'-nanmean(O_dis(:,[1 2 3 4]),1))) ;
% her.LineStyle = 'none';her.Color = 'k';
% xticklabels({'OurModel','SPM','CPM','VAM'});gca.FontSize=12;
% ylabel('Deviation Distance(m)','fontsize',12);grid off;box off;

[~,~,CI(1,:)]=ttest(O_dis(:,1));
[~,~,CI(2,:)]=ttest(O_dis(:,3));
[~,~,CI(3,:)]=ttest(O_dis(:,2));
figure
bar(1:size(O_dis(:,[1 3 2]),2),nanmean(O_dis(:,[1 3 2]),1),'FaceColor',[.5 .5 .5]);hold on
% plot([0 size(O_dis(:,[1 3 4]),2)+1],[nanmean(O_dis_stength_I_randperm,'all') nanmean(O_dis_stength_I_randperm,'all')],'k');hold on;
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
result = O_dis(:,[1 3 2]);

% MD_CI(O_dis_stength_I_randperm);
end

