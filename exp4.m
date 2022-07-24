%% ImportData [Experiment 5]
clear all;clc;
ExpName = 'Exp5a';
tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));
addpath(genpath('./Functions'));

FilesName = dir(['.\Raw Data\' ExpName filesep,'*.txt']);
for subi = 1:size(FilesName,1)
    SubData{subi}= importfile_Exp5a(['.\Raw Data\' ExpName filesep,FilesName(subi).name]);
end

DistanceA = [0.5 1];DistanceB = [0.5 1];ORientaitonA = [20, 40 ,60];ORientaitonB = [20, 40 ,60];
for i = 1:size(SubData,2)
    index_count = 1;
    for DistanceA_i = 1:size(DistanceA,2)
        for DistanceB_i = 1:size(DistanceB,2)
            for ORientaitonA_i = 1:size(ORientaitonA,2)
                for ORientaitonB_i = 1:size(ORientaitonB,2)
                    temp = SubData{i}(find(SubData{i}.z>1),:);
                    index = find(temp.DistanceA == DistanceA(DistanceA_i) &temp.DistanceB == DistanceB(DistanceB_i)...
                        &temp.ORientaitonA == ORientaitonA(ORientaitonA_i)&temp.ORientaitonB == ORientaitonB(ORientaitonB_i));
                    if (index_count == 1)
                        SubMean{i}.x = nanmean(temp.x(index'));
                        SubMean{i}.y = nanmean(temp.y(index'));
                        SubMean{i}.z = nanmean(temp.z(index'));
                        SubMean{i}.Ortation_x = nanmean(temp.Ortation_x(index'));
                        SubMean{i}.Ortation_y = nanmean(temp.Ortation_y(index'));
                        SubMean{i}.Ortation_z = nanmean(temp.Ortation_z(index'));
                        SubMean{i}.DistanceA = DistanceA(DistanceA_i);SubMean{i}.DistanceB = DistanceB(DistanceB_i);
                        SubMean{i}.ORientaitonA = ORientaitonA(ORientaitonA_i);SubMean{i}.ORientaitonB = ORientaitonB(ORientaitonB_i);
                        SubMean{i}.Sub = i;
                    else
                        SubMean{i}.x = [SubMean{i}.x ; nanmean(temp.x(index'))];
                        SubMean{i}.y = [SubMean{i}.y ; nanmean(temp.y(index'))];
                        SubMean{i}.z = [SubMean{i}.z ; nanmean(temp.z(index'))];
                        SubMean{i}.Ortation_x = [SubMean{i}.Ortation_x ; nanmean(temp.Ortation_x(index'))];
                        SubMean{i}.Ortation_y = [SubMean{i}.Ortation_y ; nanmean(temp.Ortation_y(index'))];
                        SubMean{i}.Ortation_z = [SubMean{i}.Ortation_z ; nanmean(temp.Ortation_z(index'))];
                        SubMean{i}.DistanceA = [SubMean{i}.DistanceA;DistanceA(DistanceA_i)];SubMean{i}.DistanceB = [SubMean{i}.DistanceB;DistanceB(DistanceB_i)];
                        SubMean{i}.ORientaitonA = [SubMean{i}.ORientaitonA;ORientaitonA(ORientaitonA_i)];SubMean{i}.ORientaitonB = [SubMean{i}.ORientaitonB;ORientaitonB(ORientaitonB_i)];
                        SubMean{i}.Sub = [SubMean{i}.Sub; i];
                    end
                    index_count = index_count+1;
                end
            end
        end
    end
end
for i = 1:size(SubData,2)
    if (i == 1)
        tempx = SubMean{i}.x;
        tempz = SubMean{i}.z;
        tempr = SubMean{i}.Ortation_y;
        tempDistanceA = SubMean{i}.DistanceA;tempDistanceB = SubMean{i}.DistanceB;
        tempORientaitonA = SubMean{i}.ORientaitonA;tempORientaitonB = SubMean{i}.ORientaitonB;
    else
        tempx = [tempx;SubMean{i}.x];
        tempz = [tempz;SubMean{i}.z];
        tempr = [tempr;SubMean{i}.Ortation_y];
        tempDistanceA = [tempDistanceA SubMean{i}.DistanceA];tempDistanceB = [tempDistanceB SubMean{i}.DistanceB];
        tempORientaitonA = [tempORientaitonA SubMean{i}.ORientaitonA];tempORientaitonB = [tempORientaitonB SubMean{i}.ORientaitonB];
    end
end
for i = 1:size(SubData,2)
    tempx_all(i,:)= SubMean{i}.x;
    tempz_all(i,:) = SubMean{i}.z;
    tempr_all(i,:) = SubMean{i}.Ortation_y;
end
disp('done!')

%% Get Social Potential
close all;clc;counter = 1;
xmax = 2; xmin = -2; ymax = 4;ymin = 0; Dis_step = 0.1;Ori_step = pi/36;

[condi_d1, condi_d2, condi_o1, condi_o2]= ndgrid(DistanceA, -DistanceB, wrapToPi(deg2rad(-ORientaitonA-90)), wrapToPi(deg2rad(ORientaitonB-90)));
condi_d1 = condi_d1(:);
condi_d2 = condi_d2(:);
condi_o1 = condi_o1(:);
condi_o2 = condi_o2(:);

VH_x= [condi_d1,condi_d2];
VH_y= [3.5.*ones(size(condi_d1)), 3.5.*ones(size(condi_d1))];
VH_o = [condi_o1,condi_o2];

for condi = 1:size(VH_x,1)
    [Potential{condi}, max_Ori1{condi}, AR2_record{condi}] = CalculatePotentialField([xmin,xmax;ymin,ymax],Dis_step,Ori_step,[VH_x(condi,:)', VH_y(condi,:)', VH_o(condi,:)'],0,1);
    clc;disp(['----- ' num2str(round(condi./size(VH_x,1)*100)) '% -----']);
end
disp(['done!   ' datestr(now,31)]);
%%
close all
for i = 1:size(SubData,2)
    if (i == 1)
        all_x = SubData{i}.x;
        all_y = SubData{i}.z;
        all_r = SubData{i}.Ortation_y;
        all_SubNum = SubData{i}.Sub_Num;
        all_DistanceA = SubData{i}.DistanceA;
        all_DistanceB = SubData{i}.DistanceB;
        all_OriA = SubData{i}.ORientaitonA;
        all_OriB = SubData{i}.ORientaitonB;
        all_Start = SubData{i}.StartPoint;
        SubNum = SubData{i}.Sub_Num(1);
    else
        all_x = [all_x;SubData{i}.x];
        all_y = [all_y;SubData{i}.z];
        all_r = [all_r;SubData{i}.Ortation_y];
        all_SubNum = [all_SubNum; SubData{i}.Sub_Num];
        all_DistanceA = [all_DistanceA; SubData{i}.DistanceA];
        all_DistanceB = [all_DistanceB; SubData{i}.DistanceB];
        all_OriA = [all_OriA; SubData{i}.ORientaitonA];
        all_OriB = [all_OriB; SubData{i}.ORientaitonB];
        all_Start = [all_Start; SubData{i}.StartPoint];
        SubNum = [SubNum SubData{i}.Sub_Num(1)];
    end
end

counter = 1;
for ORientaitonB_i = 1:size(ORientaitonB,2)
    for ORientaitonA_i = 1:size(ORientaitonA,2)
        for DistanceB_i = 1:size(DistanceB,2)
            for DistanceA_i = 1:size(DistanceA,2)
                index =  find(all_DistanceA == DistanceA(DistanceA_i) &all_DistanceB == DistanceB(DistanceB_i)...
                    &all_OriA == ORientaitonA(ORientaitonA_i)&all_OriB == ORientaitonB(ORientaitonB_i));
                sub_num = unique(all_SubNum);
                for subi = 1:size(SubData,2)
                    index_sub =  find(all_DistanceA == DistanceA(DistanceA_i) &all_DistanceB == DistanceB(DistanceB_i)...
                        &all_OriA == ORientaitonA(ORientaitonA_i)&all_OriB == ORientaitonB(ORientaitonB_i)& all_SubNum == SubNum(subi));
                    all_x_sub(counter,subi) = nanmean(all_x(index_sub));
                    all_y_sub(counter,subi) = nanmean(all_y(index_sub));
                    all_o_v_sub(counter,subi,:) = nanmean([cos(deg2rad((90 - all_r(index_sub)))), sin((deg2rad(90 - all_r(index_sub)) ))]);
                    all_o_sub(counter,subi) = cart2pol( all_o_v_sub(counter,subi,1) , all_o_v_sub(counter,subi,2));
                end
                
                avg_pos(counter,1) = nanmean(all_x_sub(counter,:));
                avg_pos(counter,2) = nanmean(all_y_sub(counter,:));
                mean_r_v = mean([all_o_v_sub(counter,:,1);all_o_v_sub(counter,:,2)],2);
                [avg_o(counter),~]= cart2pol(mean_r_v(1), mean_r_v(2));
                
                avg_pos_x_index = round((avg_pos(counter,1)-xmin)/Dis_step+1);
                avg_pos_y_index = round((avg_pos(counter,2)-ymin)/Dis_step+1);
                [max_x_index,max_y_index] = find(Potential{counter} == max(Potential{counter},[],1:2));
                max_x = (max_x_index-1)*Dis_step+xmin;
                max_y = (max_y_index-1)*Dis_step+ymin;
                [o_dis(counter),max_index] = min(sqrt((max_x - avg_pos(counter,1)).^2 + (max_y - avg_pos(counter,2)).^2));
                max_pos(counter,1) = max_x(max_index);
                max_pos(counter,2) = max_y(max_index);
                
                for subi = 1:size(SubData,2)
                    o_dis_allsub(counter,subi) = sqrt((all_x_sub(counter,subi) - max_pos(counter,1))^2 + (all_y_sub(counter,subi) - max_pos(counter,2))^2);
                    avg_pos_x_index_allsub = round((all_x_sub(counter,subi)-xmin)/Dis_step+1);
                    avg_pos_y_index_allsub = round((all_y_sub(counter,subi)-ymin)/Dis_step+1);
                    avg_raletive_angle_allsub(counter,subi) = wrapToPi(all_o_sub(counter,subi) - max_Ori1{counter}(avg_pos_x_index_allsub,avg_pos_y_index_allsub));
                end
                avg_raletive_angle(counter) = abs(nanmean(avg_raletive_angle_allsub(counter,:)));
                avg_raletive_angle_pn(counter) = (nanmean(avg_raletive_angle_allsub(counter,:)));
                %                 %%%
                %                 IR_VH(counter) = Interaction_Function('IR2',deg2rad(90-ORientaitonB(ORientaitonB_i)),deg2rad(90-ORientaitonA(ORientaitonA_i)),(DistanceA(DistanceA_i)+DistanceB(DistanceB_i)) );
                %                 A(counter) = abs((ORientaitonB(ORientaitonB_i))-(ORientaitonA(ORientaitonA_i)));
                counter = counter +1;
            end
        end
    end
end

%%%%%%%%%%%%% Fig3_c-a  o-distance from maximal   %%%%%%%%%%%%%%%%%
figure
% get shuffled data
for tryi = 1:1000 % turn of shuffles
    randorder = randperm(length(all_x));
    for i=1:size(avg_pos,1)
        avg_pos_rand = [all_x(randorder(i)), all_y(randorder(i))];
        o_dis_rand(i,tryi) = sqrt( (avg_pos_rand(1) - max_pos(i,1)).^2 + (avg_pos_rand(2) - max_pos(i,2)).^2);
    end
end
c = parula(10);
subplot(1,3,1:2)
h(1)=scatter(1:length(o_dis),o_dis,40,c(6,:),'linewidth',2);hold on
he=errorbar(1:length(o_dis),o_dis,nanstd(o_dis_allsub,[],2)./sqrt(size(o_dis_allsub,2)),'k');

h(2)=plot([1,length(o_dis)],[nanmean(o_dis),nanmean(o_dis)],'color',[c(6,:),0.5],'LineWidth',4);
h(3)=plot([1,length(o_dis)],[nanmean(o_dis_rand,1:2),nanmean(o_dis_rand,1:2)],'color',[c(9,:),0.5],'LineWidth',4);
he.LineStyle = 'none';

xlabel({'Conditions'});xlim([0 37]);ylim([-0.08 1.3]);
ylabel('Centrifugal Distance(m)');grid off;box off;
legend(h,{'Experiment','Experiment Average','Scrambled Average'},'Location','best');
subplot(1,3,3)
bar(1:2,[nanmean(o_dis,1:2),nanmean(o_dis_rand,1:2)],0.6,'FaceColor',[.5 .5 .5]);hold on
[~,~,CI1]=ttest(o_dis(:));
[~,~,CI2]=ttest(o_dis_rand(:));
he2 = errorbar(1:2,[nanmean(o_dis),nanmean(o_dis_rand,1:2)],[(CI1(2)-CI1(1))/2 (CI2(2)-CI2(1))/2],'k' );hold on
% he2 = errorbar(1:2,[nanmean(o_dis),nanmean(o_dis_rand,1:2)],[nanstd(o_dis)./sqrt(size(o_dis_allsub,2)),nanstd(mean(o_dis_rand,1))./sqrt(size(o_dis_rand,2))],'k' );hold on
he2.LineStyle = 'none';
xticklabels({'Expeiment','Scrambled'});
ylabel('Centrifugal Distance(m)');grid off;box off;
set(gcf,'position',[0,0,900,400]);

%%%%%%%%%%%%% Fig3_c-b  avg_raletive_angles from the real orientation   %%%%%%%%%%%%%%%%%
% get shuffled data
for tryi = 1:1000 % turn of shuffles
    randorder = randperm(length(all_x));
    for i=1:size(avg_pos,1)
        avg_v_rand = [cosd(90-all_r(randorder(i))), sind(90-all_r(randorder(i)))];
        [rand_o,~]= cart2pol(avg_v_rand(1), avg_v_rand(2));
        rand_pos_x_index = round((all_x(randorder(i))-xmin)/Dis_step+1);
        rand_pos_y_index = round((all_y(randorder(i))-ymin)/Dis_step+1);
        avg_raletive_angle_rand(i,tryi) = abs(wrapToPi(rand_o - max_Ori1{i}(rand_pos_x_index,rand_pos_y_index)));
    end
end

figure
subplot(1,3,1:2)
h2(1)=scatter(1:36,rad2deg(avg_raletive_angle),40,c(6,:),'linewidth',2);hold on
h2(2)=plot([0,37],[rad2deg(nanmean(avg_raletive_angle)),rad2deg(nanmean(avg_raletive_angle))],'color',[c(6,:),0.5],'LineWidth',4);
h2(3)=plot([0,37],[rad2deg(nanmean(avg_raletive_angle_rand,1:2)),rad2deg(nanmean(avg_raletive_angle_rand,1:2))],'color',[c(9,:),0.5],'LineWidth',4);
her = errorbar(1:36,(rad2deg(avg_raletive_angle)),nanstd(rad2deg(abs(avg_raletive_angle_allsub))',[],1)./sqrt(size(o_dis_allsub,2)),'k','CapSize',3) ;
her.LineStyle = 'none';
ylim([-3 40]);
xlabel({'Conditions'});xlim([0 37]);
ylabel('Relative Angle(бу)');grid off;box off;
legend(h2,{'Experiment','Experiment Average','Scrambled Average'},'Location','best');
subplot(1,3,3)
[~,~,CI1]=ttest(rad2deg(avg_raletive_angle(:)));
[~,~,CI2]=ttest(rad2deg(avg_raletive_angle_rand(:)));
h=bar([rad2deg(mean(avg_raletive_angle)) rad2deg(mean(avg_raletive_angle_rand,1:2))],0.6,'FaceColor',[.5 .5 .5]);hold on
her = errorbar(1:2,[mean(rad2deg(avg_raletive_angle)) rad2deg(mean(avg_raletive_angle_rand,1:2))],[(CI1(2)-CI1(1))/2 (CI2(2)-CI2(1))/2],'k') ;
% her = errorbar(1:2,[mean(rad2deg(avg_raletive_angle)) rad2deg(mean(avg_raletive_angle_rand,1:2))],[nanstd(rad2deg(mean((avg_raletive_angle_allsub),1))',[],1)./sqrt(size(o_dis_allsub,2))  nanstd(rad2deg(mean((avg_raletive_angle_rand),1))',[],1)./sqrt(size(avg_raletive_angle_rand,2))],'k') ;
her.LineStyle = 'none';
ylabel('Relative Angle(бу)');grid off;box off;
xticklabels({'Expeiment','Scrambled'});
set(gcf,'position',[0,0,900,400]);

%%%%%%%%%%%%%  statistics  %%%%%%%%%%%%%%%%%
[h,p,ci,stats{1}] = ttest(mean(o_dis_rand,2)',o_dis );
stats{1}.ci = ci;stats{1}.p = p;
stats{1}.d = ES_Cohens_d(o_dis ,mean(o_dis_rand,2)');
disp(stats{1});
[h,p,ci,stats{2}] = ttest(mean(avg_raletive_angle_rand,2)', avg_raletive_angle);
stats{2}.ci = ci;stats{2}.p = p;
stats{2}.d = ES_Cohens_d(avg_raletive_angle, mean(avg_raletive_angle_rand,2)');
disp(stats{2});

%%%%%%%%%%%%%  brief report  %%%%%%%%%%%%%%%%%
MD_CI(o_dis,1);
MD_CI(rad2deg(avg_raletive_angle),1);
MD_CI(o_dis_rand,1);
MD_CI(rad2deg(avg_raletive_angle_rand),1);

% disp(['exp_dis:' num2str(nanmean(o_dis,1:2))  ' SD=' num2str(round(nanstd(o_dis),2))  '       exp_ang:'  num2str(rad2deg(mean(avg_raletive_angle)))  ' SD=' num2str(round(nanstd(rad2deg(avg_raletive_angle)),2))])
% disp(['Shuffle_dis:' num2str(nanmean(o_dis_rand,1:2)) '   Shuffle_ang:' num2str(rad2deg(nanmean(avg_raletive_angle_rand,1:2)))])

%% SI Fig1-a & Fig1-b [all 36 conditions plot]
counter =1;clear VH_x VH_y VH_o
ix_seq = xmin:Dis_step:xmax;
iy_seq = ymin:Dis_step:ymax;
for ORientaitonB_i = 1:size(ORientaitonB,2)
    for ORientaitonA_i = 1:size(ORientaitonA,2)
        for DistanceB_i = 1:size(DistanceB,2)
            for DistanceA_i = 1:size(DistanceA,2)
                index =  find(tempDistanceA == DistanceA(DistanceA_i) &tempDistanceB == DistanceB(DistanceB_i)...
                    &tempORientaitonA == ORientaitonA(ORientaitonA_i)&tempORientaitonB == ORientaitonB(ORientaitonB_i));
                
                x1_pos = DistanceA(DistanceA_i); x2_pos = -DistanceB(DistanceB_i); y1_pos = 3.5; y2_pos = 3.5;
                theta1 = -ORientaitonA(ORientaitonA_i)-90; theta2 = ORientaitonB(ORientaitonB_i)-90;
                VH_x{counter}= [DistanceA(DistanceA_i),-DistanceB(DistanceB_i)];
                VH_y{counter}= [y1_pos,y1_pos];
                VH_o{counter} = [theta1,theta2];
                figure(1)
                subplot(6,6,counter)
                %                 surf(ix_seq,iy_seq,(Potential{counter}./ max(Potential{counter},[],'all') )'-1);hold on
                surf(ix_seq,iy_seq,(Potential{counter})'-1);hold on
                scatter(nanmean(all_x_sub(counter,:)),nanmean(all_y_sub(counter,:)),30,[0 0 0],'filled');hold on
                scatter(VH_x{counter},VH_y{counter},30,'k','filled');hold on
                quiver(VH_x{counter},VH_y{counter},cosd(VH_o{counter}),sind(VH_o{counter}),.5,'k');hold on
                grid on;axis equal;view(0,90);shading interp;
                yticks(0:1:4);xticklabels({'0','1','2','3','4'});xticks(-2:1:2);yticklabels({'0','1','2','3','4'});box off;ylim([0 4]);xlim([-2 2]);
                caxis([0.3  max(Potential{counter},[],'all')]-1);   a=colorbar;a.TickLabels = a.Ticks+1;
                %                 title(['Condition' num2str(counter)]);
                title({['A:' num2str(VH_x{counter}(2)) 'm, ' num2str(VH_o{counter}(2)) 'бу'],[ 'B:' num2str(VH_x{counter}(1)) 'm, ' num2str(VH_o{counter}(1)) 'бу']});
                set(gcf,'position',[0,0,1400,1000]);
                colormap(colorbar_BuOr());
                
                figure(2)
                subplot(6,6,counter)
                temp_Rho = squeeze(AR2_record{counter}(round(((nanmean(all_x_sub(counter,:)))-xmin)/Dis_step+1),round(((nanmean(all_y_sub(counter,:)))-ymin)/Dis_step+1),:));
                for ix = 1:length(ix_seq)
                    for iy = 1:length(iy_seq)
                        [Theta,Rho]= cart2pol((ix_seq(ix)-nanmean(all_x_sub(counter,:))),(iy_seq(iy)-nanmean(all_y_sub(counter,:))));
                        ori_strength(ix,iy,counter)=temp_Rho( round((Theta)/Ori_step + pi/Ori_step + 1) );
                    end
                end
                %                 surf(iy_seq,iy_seq,(ori_strength(:,:,counter)'./max(ori_strength(:,:,counter),[],'all'))-1);hold on;shading interp
                surf(iy_seq,iy_seq,(ori_strength(:,:,counter)')-1);hold on;shading interp
                scatter(nanmean(all_x_sub(counter,:))+2,nanmean(all_y_sub(counter,:)),30,[0 0 0],'filled');hold on
                plot([nanmean(all_x_sub(counter,:))+2*cos(nanmean(wrapToPi((all_o_sub(counter,:)))) ),nanmean(all_x_sub(counter,:))]+2,[(nanmean(all_y_sub(counter,:)))+2*sin(nanmean(wrapToPi((all_o_sub(counter,:))))),(nanmean(all_y_sub(counter,:)))],'--k','linewidth',2);
                grid on;axis equal;view(0,90);
                xticks(0:1:4);xticklabels({'0','1','2','3','4'});yticks(0:1:4);yticklabels({'0','1','2','3','4'});box off;
                ylim([0 4]);xlim([0 4]);caxis([0.2  max(ori_strength(:,:,counter),[],'all')]-1); a=colorbar;a.TickLabels = a.Ticks+1;
                %                 title(['Condition' num2str(counter)]);
                title({['A:' num2str(VH_x{counter}(2)) 'm, ' num2str(VH_o{counter}(2)) 'бу'],[ 'B:' num2str(VH_x{counter}(1)) 'm, ' num2str(VH_o{counter}(1)) 'бу']});
                set(gcf,'position',[0,0,1400,1000]);
                counter = counter+1;
                colormap(colorbar_BuOr());
                
            end
        end
    end
end
%%
close all
% for counter = [1,9,16,33]
for counter = [1]
    figure
    clear Potential_q
    [X,Y]=meshgrid(ix_seq,iy_seq);[Xq,Yq]=meshgrid(xmin:Dis_step/20:xmax,ymin:Dis_step/20:ymax);
    Potential_q = interp2(X,Y,Potential{counter},Xq,Yq,'spline');
    %     surf(Xq,Yq,Potential_q'./ max(Potential_q,[],'all')-1);hold
    surf(Xq,Yq,Potential_q'-1);hold
    scatter(nanmean(all_x_sub(counter,:)),nanmean(all_y_sub(counter,:)),30,[0 0 0],'filled');hold on
    scatter(VH_x{counter},VH_y{counter},30,'k','filled');hold on
    quiver(VH_x{counter},VH_y{counter},cosd(VH_o{counter}),sind(VH_o{counter}),.5,'k');hold on
    
    %     scatter(nanmean(all_x_sub(counter,:)),nanmean(all_y_sub(counter,:)),30,[0 0 0],'filled');hold on
    %     eh1=errorbar(nanmean(all_x_sub(counter,:)),nanmean(all_y_sub(counter,:)),nanstd(all_x_sub(counter,:)),'horizontal');hold on;
    %     eh1.Color = [0 ,0, 0];eh1.LineWidth = 1;
    %     eh2=errorbar(nanmean(all_x_sub(counter,:)),nanmean(all_y_sub(counter,:)),nanstd(all_y_sub(counter,:)));hold on;
    %     eh2.Color = [0 ,0, 0];eh2.LineWidth = 1;
    grid on;axis equal;view(0,90);shading interp;hold off;grid off;
    %     caxis([0.2  .7]-1);
    %     a=colorbar;a.TickLabels = a.Ticks+1; ylim([0 4]);xlim([-2 2]);
    view(0,90);
    %     set(gca,'DataAspectRatio',[1 1 .3]);
    %     lightangle(gca,40,130);material([0.72 0.55 0.3]);
    %     c = colorbar; c.Label.String = 'Normalized Potential ';c.FontSize = 12;
    %     c.TickLabels = c.Ticks+0.8;
    colormap(colorbar_BuOr());
end
%%
close all
Dis_step=0.05
ix_seq = xmin:Dis_step:xmax;
iy_seq = ymin:Dis_step:ymax;
for counter = [1,9,16,33]
    [~, ~, ori_record{counter}] = CalculatePotentialField([xmin,xmax;ymin,ymax],Dis_step,pi/180,[VH_x{counter}', VH_y{counter}', VH_o{counter}'],0,1);
    figure
    temp_Rho = squeeze(ori_record{counter}(round(((nanmean(all_x_sub(counter,:)))-xmin)/Dis_step+1),round(((nanmean(all_y_sub(counter,:)))-ymin)/Dis_step+1),:));
    for ix = 1:length(ix_seq)
        for iy = 1:length(iy_seq)
            [Theta,Rho]= cart2pol((ix_seq(ix)-nanmean(all_x_sub(counter,:))),(iy_seq(iy)-nanmean(all_y_sub(counter,:))));
            ori_strength(ix,iy,counter)=temp_Rho( round((Theta)/(pi/180) + pi/(pi/180) + 1) );
        end
    end
    %     surf(iy_seq,iy_seq,(ori_strength(:,:,counter)'./max(ori_strength(:,:,counter),[],'all'))-1);hold on;shading interp
    %     surf(iy_seq,iy_seq,ori_strength(:,:,counter)');hold on;shading interp
    scatter(nanmean(all_x_sub(counter,:))+2,nanmean(all_y_sub(counter,:)),30,[0 0 0],'filled');hold on
    plot([nanmean(all_x_sub(counter,:))+2*cos(nanmean(wrapToPi((all_o_sub(counter,:)))) ),nanmean(all_x_sub(counter,:))]+2,[(nanmean(all_y_sub(counter,:)))+2*sin(nanmean(wrapToPi((all_o_sub(counter,:))))),(nanmean(all_y_sub(counter,:)))],'--k','linewidth',2);
    grid on;axis equal;view(0,90);
    %     xticks(0:1:4);xticklabels({'0','1','2','3','4(m)'});yticks(0:1:4);yticklabels({'0','1','2','3','4(m)'});box off;
    ylim([0 4]);xlim([0 4]);
    %      caxis([-0.7  0]); a=colorbar;a.TickLabels = a.Ticks+1;
    title(['Condition' num2str(counter)]);
end

%%
figure
bar(1:2,[nanmean(o_dis,1:2),nanmean(o_dis_rand,1:2)],0.6,'FaceColor',[.5 .5 .5]);hold on
h1 = errorbar(1:2,[nanmean(o_dis),nanmean(o_dis_rand,1:2)],[nanstd(o_dis)./sqrt(size(o_dis_allsub,2)),nanstd(mean(o_dis_rand,1))./sqrt(size(o_dis_rand,2))],'k' );hold on
h1.LineStyle = 'none';h1.CapSize = 12;
xticklabels({'Expeiment','Scrambled'});
ylabel('Centrifugal Distance(m)');grid off;box off;ylim([0 0.8]);
set(gcf,'position',[0,0,600,400]);
figure
h=bar([rad2deg(mean(avg_raletive_angle)) rad2deg(mean(avg_raletive_angle_rand,1:2))],0.6,'FaceColor',[.5 .5 .5]);hold on
h2 = errorbar(1:2,[mean(rad2deg(avg_raletive_angle)) rad2deg(mean(avg_raletive_angle_rand,1:2))],[nanstd(rad2deg(mean((avg_raletive_angle_allsub),1))',[],1)./sqrt(size(o_dis_allsub,2))  nanstd(rad2deg(mean((avg_raletive_angle_rand),1))',[],1)./sqrt(size(avg_raletive_angle_rand,2))],'k') ;
h2.LineStyle = 'none';h2.CapSize = 12;
ylabel('Relative Angle(бу)');grid off;box off;ylim([0 20]);
xticklabels({'Expeiment','Scrambled'});
set(gcf,'position',[0,0,600,400]);
%%
h=bar([rad2deg(mean(avg_raletive_angle)) rad2deg(mean(avg_raletive_angle_rand,1:2))],0.6,'FaceColor',[.5 .5 .5]);hold on
her = errorbar(1:2,[mean(rad2deg(avg_raletive_angle)) rad2deg(mean(avg_raletive_angle_rand,1:2))],[nanstd(rad2deg(mean((avg_raletive_angle_allsub),1))',[],1)./sqrt(size(o_dis_allsub,2))  nanstd(rad2deg(mean((avg_raletive_angle_rand),1))',[],1)./sqrt(size(avg_raletive_angle_rand,2))],'k') ;
her.LineStyle = 'none';he.CapSize = 12;
ylabel('Relative Angle(бу)');grid off;box off;
xticklabels({'Expeiment','Scrambled'});
set(gcf,'position',[0,0,600,400]);

bar(1:2,[nanmean(o_dis,1:2),nanmean(o_dis_rand,1:2)],0.6,'FaceColor',[.5 .5 .5]);hold on
he2 = errorbar(1:2,[nanmean(o_dis),nanmean(o_dis_rand,1:2)],[nanstd(o_dis)./sqrt(size(o_dis_allsub,2)),nanstd(mean(o_dis_rand,1))./sqrt(size(o_dis_rand,2))],'k' );hold on
he2.LineStyle = 'none';he2.CapSize = 12;
xticklabels({'Expeiment','Scrambled'});
ylabel('Centrifugal Distance(m)');grid off;box off;
ylim([0 20]);
set(gcf,'position',[0,0,600,400]);