function [] = DataAnalysis_Angle_scatter(SaliencyMap,pos_I,ori_I,DrawArea,Dis_step)
%%% Area to Analysis
xmin = DrawArea(1,1);
xmax = DrawArea(1,2);
ymin = DrawArea(2,1);
ymax = DrawArea(2,2);

c2 = parula(10);
for i = 1:size(SaliencyMap,2)
    tempdindex = pos_I(:,:,i);
    tempV = [cos(SaliencyMap{i}{5}(round((tempdindex(1) - xmin)./Dis_step+1),round((tempdindex(2) - ymin)./Dis_step+1))),sin(SaliencyMap{i}{5}(round((tempdindex(1) - xmin)./Dis_step+1),round((tempdindex(2) - ymin)./Dis_step+1)))];
    relativeangle_A(i,1) = wrapToPi(acos(dot([tempdindex(1)*cos(ori_I(i)) tempdindex(2)*sin(ori_I(i))],tempV)/(norm([tempdindex(1)*cos(ori_I(i)) tempdindex(2)*sin(ori_I(i))])*norm(tempV))));
    tempV2 = [cos(SaliencyMap{i}{6}(round((tempdindex(1) - xmin)./Dis_step+1),round((tempdindex(2) - ymin)./Dis_step+1))),sin(SaliencyMap{i}{6}(round((tempdindex(1) - xmin)./Dis_step+1),round((tempdindex(2) - ymin)./Dis_step+1)))];
    relativeangle_A(i,4) = wrapToPi(acos(dot([tempdindex(1)*cos(ori_I(i)) tempdindex(2)*sin(ori_I(i))],tempV2)/(norm([tempdindex(1)*cos(ori_I(i)) tempdindex(2)*sin(ori_I(i))])*norm(tempV2))));

end

for i = 1:size(SaliencyMap,2)
    SaliencyMap{i}{2}(round((pos_I(:,1,i) - xmin)./Dis_step+1),round((pos_I(:,2,i) - ymin)./Dis_step+1));
    [max_index_y,max_index_x] =find(SaliencyMap{i}{2}== max(SaliencyMap{i}{2},[],'all'));
    tempV3 =[((max_index_x-1)*Dis_step+xmin)-pos_I(:,1,i), ((max_index_y-1)*Dis_step+ymin)-pos_I(:,2,i)];
    relativeangle_A(i,2) = wrapToPi(acos(dot([tempdindex(1)*cos(ori_I(i)) tempdindex(2)*sin(ori_I(i))],tempV3)/(norm([tempdindex(1)*cos(ori_I(i)) tempdindex(2)*sin(ori_I(i))])*norm(tempV3))));
end

for i = 1:size(SaliencyMap,2)
    SaliencyMap{i}{3}(round((pos_I(:,1,i) - xmin)./Dis_step+1),round((pos_I(:,2,i) - ymin)./Dis_step+1));
    [max_index_y,max_index_x] =find(SaliencyMap{i}{3}== max(SaliencyMap{i}{3},[],'all'));
    tempV3 =[((max_index_x-1)*Dis_step+xmin)-pos_I(:,1,i), ((max_index_y-1)*Dis_step+ymin)-pos_I(:,2,i)];
    relativeangle_A(i,3) = wrapToPi(acos(dot([tempdindex(1)*cos(ori_I(i)) tempdindex(2)*sin(ori_I(i))],tempV3)/(norm([tempdindex(1)*cos(ori_I(i)) tempdindex(2)*sin(ori_I(i))])*norm(tempV3))));
end


figure
scatter(1:size(SaliencyMap,2),rad2deg(relativeangle_A(:,1)),10,c2(6,:),'linewidth',1);hold on;
% plot([1,size(SaliencyMap{i}{1},2)],[rad2deg(nanmean(relativeangle_A(i,1),1:2)),rad2deg(nanmean(relativeangle_A(i,1),1:2))],'color',[c2(6,:),0.5],'LineWidth',3);
% plot([1,size(SaliencyMap{i}{1},2)],[rad2deg(nanmean(relativeangle_A_randperm,1:3)),rad2deg(nanmean(relativeangle_A_randperm,1:3))],'color',[c2(9,:),0.5],'LineWidth',3);
ylim([-5 max(rad2deg(relativeangle_A),[],'all')+5]);xlim([-5,size(SaliencyMap,2)+5]);
xlabel({'Interaction Scenes from Dataset'});gca.FontSize=12;
ylabel('Relative Angle(бу)','fontsize',12);grid off;box off;
% legend({'Model Prediction for One Scene','Model Prediction Average','Scrambled Average'},'location','northeast','Orientation','horizontal','NumColumns',2,'fontsize',12);
set(gcf,'position',[0,200,1000,450]);

figure
scatter(1:size(SaliencyMap,2),rad2deg(relativeangle_A(:,3)),10,c2(6,:),'linewidth',1);hold on;
ylim([-5 max(rad2deg(relativeangle_A),[],'all')+5]);xlim([-5,size(SaliencyMap,2)+5]);
xlabel({'Interaction Scenes from Dataset'});gca.FontSize=12;
ylabel('Relative Angle(бу)','fontsize',12);grid off;box off;
% legend({'Model Prediction for One Scene','Model Prediction Average','Scrambled Average'},'location','northeast','Orientation','horizontal','NumColumns',2,'fontsize',12);
set(gcf,'position',[0,200,1000,450]);

figure
scatter(1:size(SaliencyMap,2),rad2deg(relativeangle_A(:,2)),10,c2(6,:),'linewidth',1);hold on;
ylim([-5 max(rad2deg(relativeangle_A),[],'all')+5]);xlim([-5,size(SaliencyMap,2)+5]);
xlabel({'Interaction Scenes from Dataset'});gca.FontSize=12;
ylabel('Relative Angle(бу)','fontsize',12);grid off;box off;
% legend({'Model Prediction for One Scene','Model Prediction Average','Scrambled Average'},'location','northeast','Orientation','horizontal','NumColumns',2,'fontsize',12);

set(gcf,'position',[0,200,1000,450]);

end