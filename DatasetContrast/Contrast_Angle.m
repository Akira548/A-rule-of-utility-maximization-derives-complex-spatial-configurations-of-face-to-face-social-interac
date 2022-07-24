function result = Contrast_Angle(SaliencyMap,pos_I,ori_I,DrawArea,Dis_step,VH_pos)
%%% Area to Analysis
xmin = DrawArea(1,1);
xmax = DrawArea(1,2);
ymin = DrawArea(2,1);
ymax = DrawArea(2,2);

c2 = parula(10);
for i = 1:size(SaliencyMap,2)
    tempdindex = pos_I(:,:,i);
    tempV = [cos(SaliencyMap{i}{5}(round((pos_I(:,1,i) - xmin)./Dis_step+1),round((pos_I(:,2,i) - ymin)./Dis_step+1))),sin(SaliencyMap{i}{5}(round((pos_I(:,1,i) - xmin)./Dis_step+1),round((pos_I(:,2,i) - ymin)./Dis_step+1)))];
    relativeangle_A(i,1) = wrapToPi(acos(dot([cos(ori_I(i)) sin(ori_I(i))],tempV)/(norm([cos(ori_I(i)) sin(ori_I(i))])*norm(tempV))));

end

% for i = 1:size(SaliencyMap,2)
%     SaliencyMap{i}{2}(round((pos_I(:,1,i) - xmin)./Dis_step+1),round((pos_I(:,2,i) - ymin)./Dis_step+1));
%     [max_index_y,max_index_x] =find(SaliencyMap{i}{2}== max(SaliencyMap{i}{2},[],'all'));
%     tempV3 =[((max_index_x-1)*Dis_step+xmin)-pos_I(:,1,i), ((max_index_y-1)*Dis_step+ymin)-pos_I(:,2,i)];
%     relativeangle_A(i,2) = wrapToPi(acos(dot([tempdindex(1)*cos(ori_I(i)) tempdindex(2)*sin(ori_I(i))],tempV3)/(norm([tempdindex(1)*cos(ori_I(i)) tempdindex(2)*sin(ori_I(i))])*norm(tempV3))));
% end

% for i = 1:size(SaliencyMap,2)
% %     SaliencyMap{i}{4}(round((pos_I(:,1,i) - xmin)./Dis_step+1),round((pos_I(:,2,i) - ymin)./Dis_step+1));
%     [max_index_x,max_index_y] =find(SaliencyMap{i}{4}== max(SaliencyMap{i}{4},[],'all'));
%     tempV2 =[((max_index_x-1)*Dis_step+xmin)-pos_I(:,1,i), ((max_index_y-1)*Dis_step+ymin)-pos_I(:,2,i)];
%     relativeangle_A(i,2) = wrapToPi(acos(dot([cos(ori_I(i)) sin(ori_I(i))],tempV2)/(norm([cos(ori_I(i)) sin(ori_I(i))])*norm(tempV2))));
% end

% %%%%%%%%%%%%%%%%%%%%%%%%%
% for i = 1:size(SaliencyMap,2)
%         tempdindex = pos_I(:,:,i);
%     %     tempV2 = [cos(SaliencyMap{i}{6}(round((pos_I(:,1,i) - xmin)./Dis_step+1),round((pos_I(:,2,i) - ymin)./Dis_step+1))),sin(SaliencyMap{i}{6}(round((pos_I(:,1,i) - xmin)./Dis_step+1),round((pos_I(:,2,i) - ymin)./Dis_step+1)))];
%     %     relativeangle_A(i,2) = wrapToPi(acos(dot([cos(ori_I(i)) sin(ori_I(i))],tempV2)/(norm([cos(ori_I(i)) sin(ori_I(i))])*norm(tempV2))));
% %     if size(VH_pos{i},1)==1
% %         theta = atan2(VH_pos{i}(2)-tempdindex(2),VH_pos{i}(1)-tempdindex(1));
% %         relativeangle_A(i,2) = wrapToPi(acos(dot([cos(ori_I(i)) sin(ori_I(i))],[cos(theta) sin(theta)])));
% %     else
%         [~,cx,cy]=SIS_circle_fit(VH_pos{i});
%         theta = atan2(cy-tempdindex(2), cx-tempdindex(1));
%         relativeangle_A(i,2) = wrapToPi(acos(dot([cos(ori_I(i)) sin(ori_I(i))],[cos(theta) sin(theta)])));
% %     end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:size(SaliencyMap,2)
        tempdindex = pos_I(:,:,i);
        q=Estimate_approaching_goal_pose(VH_pos{i});
        [~,k] = min(sqrt((tempdindex(1) - q(:,1)).^2 + (tempdindex(2) - q(:,2)).^2));
        [~,cx,cy]=SIS_circle_fit(VH_pos{i});
        theta = atan2(cy-q(k,2), cx-q(k,1));
        relativeangle_A(i,2) = wrapToPi(acos(dot([cos(ori_I(i)) sin(ori_I(i))],[cos(theta) sin(theta)])));
end

for i = 1:size(SaliencyMap,2)
%     SaliencyMap{i}{3}(round((pos_I(:,1,i) - xmin)./Dis_step+1),round((pos_I(:,2,i) - ymin)./Dis_step+1));
%     [max_index_x,max_index_y] =find(SaliencyMap{i}{3}== max(SaliencyMap{i}{3},[],'all'));
%     tempV3 =[((max_index_x-1)*Dis_step+xmin)-pos_I(:,1,i), ((max_index_y-1)*Dis_step+ymin)-pos_I(:,2,i)];
    rand_ang = deg2rad(randi([0 359]));
    tempV3 =[cos(rand_ang), sin(rand_ang)];
    relativeangle_A(i,3) = wrapToPi(acos(dot([cos(ori_I(i)) sin(ori_I(i))],tempV3)/(norm([cos(ori_I(i)) sin(ori_I(i))])*norm(tempV3))));
end

% for tryi = 1:100
%     randpermi = randperm(size(SaliencyMap,2));
%     while sum(randpermi == 1:size(SaliencyMap,2)) ~= 0
%         randpermi = randperm(size(SaliencyMap,2));
%     end
%     for i = 1:size(SaliencyMap,2)
%         tempdindex = pos_I(:,:,randpermi(i));
%         tempV = [cos(SaliencyMap{i}{5}(round((tempdindex(1) - xmin)./Dis_step+1),round((tempdindex(2) - ymin)./Dis_step+1))),sin(SaliencyMap{i}{5}(round((tempdindex(1) - xmin)./Dis_step+1),round((tempdindex(2) - ymin)./Dis_step+1)))];
%         relativeangle_A_randperm(i,tryi) = wrapToPi(acos(dot([tempdindex(1)*cos(ori_I(i)) tempdindex(2)*sin(ori_I(i))],tempV)/(norm([tempdindex(1)*cos(ori_I(i)) tempdindex(2)*sin(ori_I(i))])*norm(tempV)) ));
%     end
% end

% figure
% bar(1:size(relativeangle_A(:,[1 2 3 4]),2),rad2deg(nanmean(relativeangle_A(:,[1 2 3 4]),1)),'FaceColor',[.5 .5 .5]);hold on
% her = errorbar(1:size(relativeangle_A(:,[1 2 3 4]),2),rad2deg(nanmean(relativeangle_A(:,[1 2 3 4]),1)),nanstd(rad2deg(relativeangle_A(:,[1 2 3 4])))./sqrt(size(relativeangle_A(:,[1 2 3 4]),1))) ;
% her.LineStyle = 'none';her.Color = 'k';
% xticklabels({'OurModel','SPM','VAM','CPM'});gca.FontSize=12;
% ylabel('Deviation Angle(бу)','fontsize',12);grid off;box off;

[~,~,CI(1,:)]=ttest(relativeangle_A(:,1));
[~,~,CI(2,:)]=ttest(relativeangle_A(:,3));
[~,~,CI(3,:)]=ttest(relativeangle_A(:,2));

figure
bar(1:size(relativeangle_A(:,[1 3 2]),2),rad2deg(nanmean(relativeangle_A(:,[1 3 2]),1)),'FaceColor',[.5 .5 .5]);hold on
% plot([0 size(relativeangle_A(:,[1 3 2]),2)+1],[rad2deg(nanmean(relativeangle_A_randperm,'all')) rad2deg(nanmean(relativeangle_A_randperm,'all'))],'k');hold on;
her = errorbar(1:size(relativeangle_A(:,[1 3 2]),2),rad2deg(nanmean(relativeangle_A(:,[1 3 2]),1)),rad2deg(abs(CI(:,1)'-nanmean(relativeangle_A(:,[1 3 2]),1)))) ;
her.LineStyle = 'none';her.Color = 'k';
xticklabels({'OurModel','CPM','VAM'});gca.FontSize=12;
ylabel('Deviation Angle(бу)','fontsize',12);grid off;box off;
set(gcf,'position',[0,0,700,400]);
MD_CI(rad2deg(relativeangle_A(:,1)));
MD_CI(rad2deg(relativeangle_A(:,3)));
MD_CI(rad2deg(relativeangle_A(:,2)));
% MD_CI(rad2deg(relativeangle_A_randperm));
result = relativeangle_A(:,[1 3 2]);

end