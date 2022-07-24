function result = ExpContrast_Angle(SaliencyMap,pos_I,ori_I,DrawArea,Dis_step,VH_pos,isdraw)
%%% Area to Analysis
xmin = DrawArea(1,1);
xmax = DrawArea(1,2);
ymin = DrawArea(2,1);
ymax = DrawArea(2,2);
if nargin < 7
    isdraw = 0;
end

tempV = [cos(SaliencyMap{5}(round((pos_I(1) - xmin)./Dis_step+1),round((pos_I(2) - ymin)./Dis_step+1))),sin(SaliencyMap{5}(round((pos_I(1) - xmin)./Dis_step+1),round((pos_I(2) - ymin)./Dis_step+1)))];
relativeangle_A(1) = wrapToPi(acos(dot([cos(ori_I) sin(ori_I)],tempV)/(norm([cos(ori_I) sin(ori_I)])*norm(tempV))));


% [max_index_x,max_index_y] =find(SaliencyMap{4}== max(SaliencyMap{4},[],'all'));
% tempV2 =[((max_index_x-1)*Dis_step+xmin)-pos_I(1), ((max_index_y-1)*Dis_step+ymin)-pos_I(2)];
% if size(tempV2,1)>1
%     relativeangle_A(2) = wrapToPi(acos(dot([cos(ori_I) sin(ori_I)],tempV2(randi([1 length(size(tempV2,1))]),:))/(norm([cos(ori_I) sin(ori_I)])*norm(tempV2(randi([1 length(size(tempV2,1))]),:)))));
% else
%     relativeangle_A(2) = wrapToPi(acos(dot([cos(ori_I) sin(ori_I)],tempV2)/(norm([cos(ori_I) sin(ori_I)])*norm(tempV2))));
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tempV2 = [cos(SaliencyMap{6}(round((pos_I(1) - xmin)./Dis_step+1),round((pos_I(2) - ymin)./Dis_step+1))),sin(SaliencyMap{6}(round((pos_I(1) - xmin)./Dis_step+1),round((pos_I(2) - ymin)./Dis_step+1)))];
% relativeangle_A(2) = wrapToPi(acos(dot([cos(ori_I) sin(ori_I)],tempV2)/(norm([cos(ori_I) sin(ori_I)])*norm(tempV2))));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [max_index_x,max_index_y] =find(SaliencyMap{2}== min(SaliencyMap{2},[],'all'));
% tempV2 =[((max_index_x-1)*Dis_step+xmin)-pos_I(1), ((max_index_y-1)*Dis_step+ymin)-pos_I(2)];
% if size(tempV2,1)>1
%     relativeangle_A(2) = wrapToPi(acos(dot([cos(ori_I) sin(ori_I)],tempV2(randi([1 length(size(tempV2,1))]),:))/(norm([cos(ori_I) sin(ori_I)])*norm(tempV2(randi([1 length(size(tempV2,1))]),:)))));
% else
%     relativeangle_A(2) = wrapToPi(acos(dot([cos(ori_I) sin(ori_I)],tempV2)/(norm([cos(ori_I) sin(ori_I)])*norm(tempV2))));
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if size(VH_pos,1)==1
    theta = atan2(VH_pos(2)-pos_I(:,2),VH_pos(1)-pos_I(:,1));
    relativeangle_A(2) = wrapToPi(acos(dot([cos(ori_I) sin(ori_I)],[cos(theta) sin(theta)])));
else
    [~,cx,cy]=SIS_circle_fit(VH_pos);
    theta = atan2(cy-pos_I(:,2), cx-pos_I(:,1));
    relativeangle_A(2) = wrapToPi(acos(dot([cos(ori_I) sin(ori_I)],[cos(theta) sin(theta)])));
end


% if ~isempty(q)
%     [~,n] = min(sqrt((pos_I(1) - q(:,1)).^2 + (pos_I(2) - q(:,2)).^2));
% else
%     relativeangle_A(2) = nan;
% end
% 
% if ~isempty(q)
%     relativeangle_A(2) = wrapToPi(acos(dot([cos(ori_I) sin(ori_I)],[cos(q(n,3)) sin(q(n,3))])));
% else
%     relativeangle_A(2) = nan;
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rand_ang = deg2rad(randi([0 359]));
tempV3 =[cos(rand_ang), sin(rand_ang)];
relativeangle_A(3) = wrapToPi(acos(dot([cos(ori_I) sin(ori_I)],tempV3)/(norm([cos(ori_I) sin(ori_I)])*norm(tempV3))));
if isdraw
[~,~,CI(1,:)]=ttest(relativeangle_A(1));
[~,~,CI(2,:)]=ttest(relativeangle_A(3));
[~,~,CI(3,:)]=ttest(relativeangle_A(2));

figure
bar(1:size(relativeangle_A([1 3 2]),2),rad2deg(nanmean(relativeangle_A([1 3 2]),1)),'FaceColor',[.5 .5 .5]);hold on
her = errorbar(1:size(relativeangle_A([1 3 2]),2),rad2deg(nanmean(relativeangle_A([1 3 2]),1)),rad2deg(abs(CI(:,1)'-nanmean(relativeangle_A(:,[1 3 2]),1)))) ;
her.LineStyle = 'none';her.Color = 'k';
xticklabels({'OurModel','CPM','VAM'});gca.FontSize=12;
ylabel('Deviation Angle(бу)','fontsize',12);grid off;box off;
set(gcf,'position',[0,0,700,400]);
MD_CI(rad2deg(relativeangle_A(:,1)));
MD_CI(rad2deg(relativeangle_A(:,3)));
MD_CI(rad2deg(relativeangle_A(:,2)));

end
result = rad2deg(nanmean(relativeangle_A([1 3 2]),1));

end