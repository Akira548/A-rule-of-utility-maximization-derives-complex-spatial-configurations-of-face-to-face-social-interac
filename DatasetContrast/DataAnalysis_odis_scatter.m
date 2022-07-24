function [] = DataAnalysis_odis_scatter(SaliencyMap,pos_I,DrawArea,Dis_step,y_limits)
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
            [peak_x,peak_y] = find(SaliencyMap{i}{modi} ==  max(SaliencyMap{i}{modi},[],'all'));
            Tran_peak_x = (peak_x-1)*Dis_step+xmin;    Tran_peak_y = (peak_y-1)*Dis_step+ymin;
            if modi == 2
                O_dis(i,modi) = (min(sqrt((tempdindex(1) - Tran_peak_x).^2+ (tempdindex(2) - Tran_peak_y).^2)));
            else
                O_dis(i,modi) = min(sqrt((tempdindex(1) - Tran_peak_x).^2+ (tempdindex(2) - Tran_peak_y).^2));
            end
        end
    end
end

figure
scatter(1:size(SaliencyMap,2),(O_dis(:,1)),10,c2(6,:),'linewidth',1);hold on;
ylim([-.1 max(O_dis(:,[1 3 2]),[],'all')+.1]);xlim([-5,size(SaliencyMap,2)+5]);
xlabel({'Interaction Scenes from Dataset'});gca.FontSize=12;
ylabel('Relative Angle(бу)','fontsize',12);grid off;box off;
% legend({'Model Prediction for One Scene','Model Prediction Average','Scrambled Average'},'location','northeast','Orientation','horizontal','NumColumns',2,'fontsize',12);
set(gcf,'position',[0,200,1000,450]);

figure
scatter(1:size(SaliencyMap,2),(O_dis(:,3)),10,c2(6,:),'linewidth',1);hold on;
ylim([-.1 max((O_dis(:,[1 3 2])),[],'all')+.1]);xlim([-5,size(SaliencyMap,2)+5]);
xlabel({'Interaction Scenes from Dataset'});gca.FontSize=12;
ylabel('Relative Angle(бу)','fontsize',12);grid off;box off;
% legend({'Model Prediction for One Scene','Model Prediction Average','Scrambled Average'},'location','northeast','Orientation','horizontal','NumColumns',2,'fontsize',12);
set(gcf,'position',[0,200,1000,450]);

figure
scatter(1:size(SaliencyMap,2),(O_dis(:,2)),10,c2(6,:),'linewidth',1);hold on;
ylim([-.1 max((O_dis(:,[1 3 2])),[],'all')+.1]);xlim([-5,size(SaliencyMap,2)+5]);
xlabel({'Interaction Scenes from Dataset'});gca.FontSize=12;
ylabel('Relative Angle(бу)','fontsize',12);grid off;box off;
% legend({'Model Prediction for One Scene','Model Prediction Average','Scrambled Average'},'location','northeast','Orientation','horizontal','NumColumns',2,'fontsize',12);

set(gcf,'position',[0,200,1000,450]);

end