%% ImportData [Experiment 6]
clear all;clc;
ExpName = 'Exp5b';
tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));
addpath(genpath('./Functions'));

%%% Exp Scene Information
Scene_SubData{1}.x = [0.3368168,0.1813235,-0.6292681,-0.60147,1.067015];% Scene_1 5 VHs
Scene_SubData{1}.y = [-0.6848758,0.661832,0.5880659,-0.2623619,0.171985]+3;
Scene_SubData{1}.o = wrapToPi(deg2rad(90-[-15.45,185.55,148.55,72.85001,-77.81001]));

Scene_SubData{2}.x = [-0.3916582,-0.8496783,0.6736049,0.1417348,0.7597367];% Scene_1 5 VHs
Scene_SubData{2}.y = [0.478461,-0.2693477,0.7527134,0.8532834,-0.4642425]+3;
Scene_SubData{2}.o = wrapToPi(deg2rad(90-[128.563,99.563,226.24,-172.437,-82.637]));

Scene_SubData{3}.x = [-0.2410707,0.4437704,-0.5636601,0.6243482];   
Scene_SubData{3}.y = [-0.7449532,-0.4534287,0.2756006,0.1774372]+3;
Scene_SubData{3}.o = wrapToPi(deg2rad(90-[-10.378,-79.91,156.09,-117.718]));

Scene_SubData{4}.x = [-0.81,0.59];
Scene_SubData{4}.y = [0,-0.15]+3;
Scene_SubData{4}.o = wrapToPi(deg2rad(90-[137.495,242.877]));

Scene_SubData{5}.x = [-0.7679286,-0.06869586,0.8654903];
Scene_SubData{5}.y = [-0.09940343,0.6807702,-0.09165483]+3;
Scene_SubData{5}.o = wrapToPi(deg2rad(90-[474.544,172.97,617.113]));

FilesName = dir(['.\Raw Data\' ExpName filesep,'*.txt']);
for subi = 1:size(FilesName,1)
    tempdata = importfile_Exp5b(['.\Raw Data\' ExpName filesep,FilesName(subi).name]);
    oricondition = unique(tempdata.ORientaiton);
    Sub_pos_Data{subi}.ori = wrapToPi(deg2rad(-(tempdata.ORientaiton)));
    Sub_pos_Data{subi}.Scenes = tempdata.Scenes;
    Sub_pos_Data{subi}.x = tempdata.x;
    Sub_pos_Data{subi}.z = tempdata.z;
%     Sub_pos_Data{subi}.heading = (pi/2)-tempdata.Ortation_y;
    Sub_pos_Data{subi}.heading = tempdata.Ortation_y;
end
xmin = -6;xmax = 6;ymin = -3;ymax = 9;Dis_step = 0.1;Ori_step = pi/36;

%%
close all;clc;
walk_dir=[0-pi, pi/2, -pi/2];
for Scenei = 0:4
    VH_pos = [Scene_SubData{Scenei+1}.x', Scene_SubData{Scenei+1}.y', Scene_SubData{Scenei+1}.o'];
    [SaliencyMap(:,:,Scenei+1), max_Ori1{Scenei+1}, AR2_record{Scenei+1}] = CalculatePotentialField([xmin,xmax;ymin,ymax],Dis_step,Ori_step,VH_pos,0,1);
    figure
%     contour(xmin:Dis_step:xmax,ymin:Dis_step:ymax,(SaliencyMap(:,:,Scenei+1)'./max(SaliencyMap(:,:,Scenei+1),[],'all'))-1.1,'k');hold on
%     surf(xmin:Dis_step:xmax,ymin:Dis_step:ymax,(SaliencyMap(:,:,Scenei+1)'./max(SaliencyMap(:,:,Scenei+1),[],'all'))-1.1);hold on;view(0,90);shading interp;
%     contour(xmin:Dis_step:xmax,ymin:Dis_step:ymax,(SaliencyMap(:,:,Scenei+1)'),'k');hold on
    surf(xmin:Dis_step:xmax,ymin:Dis_step:ymax,(SaliencyMap(:,:,Scenei+1)')-1);hold on;view(0,90);shading interp;
    
    VH_pos=[Scene_SubData{Scenei+1}.x;Scene_SubData{Scenei+1}.y;Scene_SubData{Scenei+1}.o];
    scatter(VH_pos(1,:),VH_pos(2,:),400,'.k');hold on;
    quiver(VH_pos(1,:) ,VH_pos(2,:),cos(VH_pos(3,:)),sin(VH_pos(3,:)) ,.3,'.k');hold on;
    for subi = 1:size(Sub_pos_Data,2)
        temp_index = find(Sub_pos_Data{subi}.Scenes == Scenei);
        [~,sortindex]= (sort(Sub_pos_Data{subi}.ori(temp_index)));
        index = temp_index(sortindex);
        index = find(Sub_pos_Data{subi}.Scenes == Scenei);
        AllData{subi,Scenei+1} = [Sub_pos_Data{subi}.x(index).*cos(Sub_pos_Data{subi}.ori(index)) + (Sub_pos_Data{subi}.z(index)-3).*sin(Sub_pos_Data{subi}.ori(index))...
            ,(Sub_pos_Data{subi}.z(index)-3).*cos(Sub_pos_Data{subi}.ori(index))-Sub_pos_Data{subi}.x(index).*sin(Sub_pos_Data{subi}.ori(index))+3];
        AllData_mean(subi,Scenei+1,1:2) = [nanmean(AllData{subi,Scenei+1}(:,1)),nanmean(AllData{subi,Scenei+1}(:,2))];
        AllData_o{subi,Scenei+1} = wrapToPi(pi/2-Sub_pos_Data{subi}.ori(index) - deg2rad(Sub_pos_Data{subi}.heading(index)));
        scatter(AllData{subi,Scenei+1}(:,1),AllData{subi,Scenei+1}(:,2),50,'.k');hold on
%         quiver(AllData{subi,Scenei+1}(:,1),AllData{subi,Scenei+1}(:,2),cos(AllData_o{subi,Scenei+1}),sin(AllData_o{subi,Scenei+1}),.1,'.k');hold on  
        temp = AllData{subi,Scenei+1};
        pos_index = [round((temp(:,1)-xmin)./Dis_step+1)  round((temp(:,2)-ymin)./Dis_step+1)]; 
%         pos_index = [ round((temp(:,2)-ymin)./Dis_step+1) round((temp(:,1)-xmin)./Dis_step+1) ]; 
        for pos_index_i = 1: size(pos_index,1)
            model_o(pos_index_i) = wrapToPi(max_Ori1{Scenei+1}(pos_index(pos_index_i,1),pos_index(pos_index_i,2)));
            RA_o{subi,Scenei+1}(pos_index_i) =wrapToPi(AllData_o{subi,Scenei+1}(pos_index_i))-...
                wrapToPi(max_Ori1{Scenei+1}(pos_index(pos_index_i,1),pos_index(pos_index_i,2)));
        end
%         quiver(AllData{subi,Scenei+1}(:,1),AllData{subi,Scenei+1}(:,2),cos(model_o)',sin(model_o)',.1,'.r');hold on  
        RA_mean_exp5b(subi,Scenei+1) =mean( abs(wrapToPi(RA_o{subi,Scenei+1})));
        
        axis equal;ylim([ymin ymax]);xlim([xmin xmax]);box off;grid off;
%         caxis([-1.1,-0.1]);
    end
    caxis([.2  max(SaliencyMap(:,:,Scenei+1)',[],'all')]-1);
    colormap(colorbar_BuOr());
    xlabel('Distance(m)');ylabel('Distance(m)');
    a=colorbar;a.TickLabels = a.Ticks+1;
    WinSize = get(0);
    set(gcf, 'position', WinSize.ScreenSize.*0.5);
    [peak_x_index,peak_y_index,~] = find(imregionalmax((SaliencyMap(:,:,Scenei+1)./max(SaliencyMap(:,:,Scenei+1),[],'all'))));
    peak_x = (peak_x_index-1).*Dis_step+xmin ;peak_y = (peak_y_index-1).*Dis_step+ymin;
    Data_allsub_x = [];Data_allsub_y=[];
    Localmaximal{Scenei+1} = [peak_x, peak_y];
    for i = 1:length(peak_x)
%         peak_z{Scenei+1}(i) = (SaliencyMap(peak_x_index(i),peak_y_index(i),Scenei+1)./max(SaliencyMap(:,:,Scenei+1),[],'all'));
        peak_z{Scenei+1}(i) = (SaliencyMap(peak_x_index(i),peak_y_index(i),Scenei+1));
        for subi = 1:size(Sub_pos_Data,2)
            Data_allsub_x(subi,:) = AllData{subi,Scenei+1}(:,1);
            Data_allsub_y(subi,:) = AllData{subi,Scenei+1}(:,2);
        end
        peak_odis(i,:) = sqrt( (Data_allsub_x(:)-peak_x(i)).^2 + (Data_allsub_y(:)-peak_y(i)).^2 )';
    end
    
    for i = 1:length(peak_x)
        [~,Group] = min(peak_odis,[],1);
        Groupratio{Scenei+1}(i) = sum(Group==i);
    end
%     Groupratio_original(Scenei+1) = max(Groupratio{Scenei+1})./sum(Groupratio{Scenei+1});
%     Groupratio_original(Scenei+1) = max(Groupratio{Scenei+1})./sum(Groupratio{Scenei+1});
%     Groupratio{Scenei+1} = Groupratio{Scenei+1}./max(Groupratio{Scenei+1});
    clc;disp(['now:' num2str((Scenei)/4*100)  '%']);
end
% rad2deg(mean(RA_mean_exp5b,'all'))

figure
for Scenei = 0:4
    subplot(1,5,Scenei+1)
    bar([peak_z{Scenei+1};Groupratio{Scenei+1}],'FaceColor','flat');hold on
    xticklabels({'Potential','Preference'});ylabel('Normalized Ratio');
    title(['Scene ' num2str(Scenei+1)])
end
set(gcf, 'position', [100 100 1600 300]);

fit_peak = []; fit_ratio = [];
for Scenei = 0:4
    fit_peak = [fit_peak,peak_z{Scenei+1}];
    fit_ratio = [fit_ratio,Groupratio{Scenei+1}];
end

figure
c = parula(10);
f=fit(fit_peak',fit_ratio','poly1');
h=plot(f,fit_peak',fit_ratio','ok');
h(1).LineWidth = 1;
h(2).LineStyle = '--';
h(2).LineWidth = 1;
h(2).Color = c(6,:);
text(0.8,0,['r = ' num2str(corr(fit_peak',fit_ratio'))]);
[r,p] = corr(fit_peak',fit_ratio');
disp(['r = ' num2str(r) '  p = ' num2str(p)]);
legend('off');
% xlim([0,1.02]);
ylim([-0.2,1.02]);
xlabel('Potential');ylabel('Preference');
box('off');
set(gcf,'position',[0,0,600,400]);

%% statistic
clc;
for Scenei = 0:4
    for subi = 1:size(Sub_pos_Data,2)
        temp = [];
        for maxi = 1:size(Localmaximal{Scenei+1},1)
            temp   = [temp, sqrt((AllData{subi,Scenei+1}(:,1) - Localmaximal{Scenei+1}(maxi,1)).^2 + (AllData{subi,Scenei+1}(:,2) - Localmaximal{Scenei+1}(maxi,2)).^2)];
        end
        deviation_dis(subi,Scenei+1) = mean(min(temp,[],2));
        deviation_dis_all(subi,:,Scenei+1) = min(temp,[],2);
        deviation_ang_all(subi,:,Scenei+1) = RA_o{subi,Scenei+1};
    end
end

% disp(['M:' num2str(round(mean(deviation_dis_all,1:3),3)) '   SD:' num2str(round(std(mean(deviation_dis_all,1:2)),3))]);
[tbl,~]=simple_mixed_anova(deviation_dis_all,[],{'Ori','Scenes'});
partial_Eta_squared1 = ES_Eta_squared(tbl,2);
[tbl,~]=simple_mixed_anova(rad2deg(deviation_ang_all));
partial_Eta_squared1 = ES_Eta_squared(tbl,2);

deviation_dis_exp5b = deviation_dis;
figure
h=bar([mean(deviation_dis_exp5b,'all')],0.6,'FaceColor',[.5 .5 .5]);hold on
her = errorbar(1,[mean(deviation_dis_exp5b,'all')],[nanstd(mean(deviation_dis_exp5b,1))./sqrt(size(mean(deviation_dis_exp5b,1),2))],'k') ;
her.LineStyle = 'none';
ylim([0 0.8]);
ylabel('Deviation Distance (m)');grid off;box off;
xticklabels({'VR','Real'});
set(gcf,'position',[0,0,600,400]);

figure
h=bar(1,[rad2deg(mean(RA_mean_exp5b,'all'))],0.6,'FaceColor',[.5 .5 .5]);hold on
her = errorbar(1,[rad2deg(mean(RA_mean_exp5b,'all'))],[nanstd(mean(RA_mean_exp5b,1))./sqrt(size(mean(RA_mean_exp5b,1),2))],'k') ;
her.LineStyle = 'none';
ylim([0 23]);
ylabel('Relative Angle(бу)');grid off;box off;
xticklabels({'VR','Real'});
set(gcf,'position',[0,0,600,400]);

%%%%%%%%%%%%%  brief report  %%%%%%%%%%%%%%%%%
MD_CI(deviation_dis_exp5b,1);
MD_CI(rad2deg(RA_mean_exp5b),1);
% disp(['exp_dis:' num2str(nanmean(deviation_dis_exp5b,1:2))  ' SD=' num2str(round(nanstd(deviation_dis_exp5b,[],'all'),2))...
%     '       exp_ang:'  num2str(rad2deg(mean(RA_mean_exp5b,'all')))  ' SD=' num2str(round(nanstd(rad2deg(RA_mean_exp5b),[],'all'),2))])

%%
all_x = [];
for Scenei = 0:4
    for subi = 1:size(Sub_pos_Data,2)
        if isempty(all_x)
            all_x = AllData{subi,Scenei+1}(:,1);
            all_y = AllData{subi,Scenei+1}(:,2);
            all_o = wrapToPi(AllData_o{subi,Scenei+1})
        else
            all_x = [all_x; AllData{subi,Scenei+1}(:,1)];
            all_y = [all_y; AllData{subi,Scenei+1}(:,2)];
            all_o = [all_o; wrapToPi(AllData_o{subi,Scenei+1})];
        end
    end
end

for tryi = 1:1000 % turn of shuffles
    for Scenei = 0:4
        [peak_x_index,peak_y_index,~] = find(imregionalmax((SaliencyMap(:,:,Scenei+1)./max(SaliencyMap(:,:,Scenei+1),[],'all'))));
        peak_x = (peak_x_index-1).*Dis_step+xmin ;peak_y = (peak_y_index-1).*Dis_step+ymin;
        Data_allsub_x = [];Data_allsub_y=[];
        Localmaximal{Scenei+1} = [peak_x, peak_y];
        randindex = randi([1 length(all_x)]);
        for i = 1:length(peak_x)
            peak_z{Scenei+1}(i) = (SaliencyMap(peak_x_index(i),peak_y_index(i),Scenei+1));
            avg_pos_rand = [all_x(randindex), all_y(randindex)];

            o_dis_rand(Scenei+1,tryi) = sqrt( (avg_pos_rand(1)-peak_x(i)).^2 + (avg_pos_rand(2)-peak_y(i)).^2 )';
        end
        avg_v_rand = [cosd(90-all_o(randindex)), sind(90-all_o(randindex))];
        [rand_o,~]= cart2pol(avg_v_rand(1), avg_v_rand(2));
        rand_pos_x_index = round((all_x(randindex)-xmin)/Dis_step+1);
        rand_pos_y_index = round((all_y(randindex)-ymin)/Dis_step+1);
        avg_raletive_angle_rand(Scenei+1,tryi) = abs(wrapToPi(rand_o - max_Ori1{Scenei+1}(rand_pos_x_index,rand_pos_y_index)));
    end
end
figure
h=bar(1:2,[mean(deviation_dis_exp5b,'all') mean(o_dis_rand,'all')],0.6,'FaceColor',[.5 .5 .5]);hold on
[~,~,CI1]=ttest(deviation_dis_exp5b(:));
[~,~,CI2]=ttest(o_dis_rand(:));
her = errorbar(1:2,[mean(deviation_dis_exp5b,'all') mean(o_dis_rand,'all')],[(CI1(2)-CI1(1))/2 (CI2(2)-CI2(1))/2],'k') ;
her.LineStyle = 'none';
% ylim([0 0.8]);
ylabel('Deviation Distance (m)');grid off;box off;
xticklabels({'Expeiment','Scrambled'});
set(gcf,'position',[0,0,600,400]);

figure
h=bar(1:2,[rad2deg(mean(RA_mean_exp5b,'all')) rad2deg(mean(avg_raletive_angle_rand,'all'))],0.6,'FaceColor',[.5 .5 .5]);hold on
[~,~,CI1]=ttest(rad2deg(RA_mean_exp5b(:)));
[~,~,CI2]=ttest(rad2deg(avg_raletive_angle_rand(:)));
her = errorbar(1:2,[rad2deg(mean(RA_mean_exp5b,'all')) rad2deg(mean(avg_raletive_angle_rand,'all'))],[(CI1(2)-CI1(1))/2 (CI2(2)-CI2(1))/2],'k') ;
her.LineStyle = 'none';
% ylim([0 23]);
ylabel('Relative Angle(бу)');grid off;box off;
xticklabels({'Expeiment','Scrambled'});
set(gcf,'position',[0,0,600,400]);

%%%%%%%%%%%%%  brief report  %%%%%%%%%%%%%%%%%
MD_CI(deviation_dis_exp5b,1);
MD_CI(rad2deg(RA_mean_exp5b),1);

%%
figure
% subplot(1,3,1:2)
h2_(1)=scatter(1:5,squeeze(mean(deviation_dis_all(:,1,:),1)),40,c(6,:),'linewidth',2);hold on
h2_(1)=scatter(1:5,squeeze(mean(deviation_dis_all(:,2,:),1)),40,c(6,:),'linewidth',2);hold on
h2_(1)=scatter(1:5,squeeze(mean(deviation_dis_all(:,3,:),1)),40,c(6,:),'linewidth',2);hold on
h2_(1)=scatter(1:5,squeeze(mean(deviation_dis_all(:,4,:),1)),40,c(6,:),'linewidth',2);hold on

h2_(2)=plot([0.5,5.5],[mean(deviation_dis_all,'all'),mean(deviation_dis_all,'all')],'color',[c(6,:),0.5],'LineWidth',4);

her1 = errorbar(1:5,squeeze(mean(deviation_dis_all(:,1,:),1)),squeeze(nanstd(deviation_dis_all(:,1,:)))./sqrt(size(deviation_dis_all,1)),'k','CapSize',3) ;hold on
her1.LineStyle = 'none';
her2 = errorbar(1:5,squeeze(mean(deviation_dis_all(:,2,:),1)),squeeze(nanstd(deviation_dis_all(:,1,:)))./sqrt(size(deviation_dis_all,1)),'k','CapSize',3) ;hold on
her2.LineStyle = 'none';
her3 = errorbar(1:5,squeeze(mean(deviation_dis_all(:,3,:),1)),squeeze(nanstd(deviation_dis_all(:,1,:)))./sqrt(size(deviation_dis_all,1)),'k','CapSize',3) ;hold on
her3.LineStyle = 'none';
her4 = errorbar(1:5,squeeze(mean(deviation_dis_all(:,4,:),1)),squeeze(nanstd(deviation_dis_all(:,1,:)))./sqrt(size(deviation_dis_all,1)),'k','CapSize',3) ;hold on
her4.LineStyle = 'none';
xticks(1:5);xlabel("Scene");ylabel("deviation distance");

% ylim([-25 25]);
% xlabel({'Conditions'});xlim([0 37]);
% ylabel('Relative Angle(бу)');grid off;box off;
% legend(h2_,{'Experiment','Experiment Average'},'Location','best');
% subplot(1,3,3)
% h=bar([rad2deg(mean(avg_raletive_angle)) rad2deg(mean(avg_raletive_angle_rand_abs,1:2))],0.6,'FaceColor',[.5 .5 .5]);hold on
% her = errorbar(1:2,[mean(rad2deg(avg_raletive_angle)) rad2deg(mean(avg_raletive_angle_rand_abs,1:2))],[nanstd(rad2deg(mean((avg_raletive_angle_allsub),1))',[],1)./sqrt(size(o_dis_allsub,2))  nanstd(rad2deg(mean((avg_raletive_angle_rand),1))',[],1)./sqrt(size(o_dis_rand,2))],'k') ;
% her.LineStyle = 'none';
% [h,p,ci] = ttest([avg_raletive_angle mean(avg_raletive_angle_rand_abs,2)']);
% ylabel('Relative Angle(бу)');grid off;box off;
% xticklabels({'Expeiment','Scrambled'});
% set(gcf,'position',[0,0,900,400]);
%%
xmin = -6;xmax = 6;ymin = -3;ymax = 9;Dis_step = 0.1;Ori_step = pi/36;

for Scenei = 0:4
    
    VH_pos = [Scene_SubData{Scenei+1}.x', Scene_SubData{Scenei+1}.y', Scene_SubData{Scenei+1}.o'];
    [SaliencyMap(:,:,Scenei+1), max_Ori1{Scenei+1}, AR2_record{Scenei+1}] = CalculatePotentialField([xmin,xmax;ymin,ymax],Dis_step,Ori_step,VH_pos,0,1);
    
    figure
    contourf(xmin:Dis_step:xmax,ymin:Dis_step:ymax,(SaliencyMap(:,:,Scenei+1)'./max(SaliencyMap(:,:,Scenei+1),[],'all'))-1.1,'k');hold on
    VH_pos=[Scene_SubData{Scenei+1}.x;Scene_SubData{Scenei+1}.y;Scene_SubData{Scenei+1}.o];
    for subi = 1:size(Sub_pos_Data,2)
        index = find(Sub_pos_Data{subi}.Scenes == Scenei);
        AllData{subi,Scenei+1} = [Sub_pos_Data{subi}.x(index).*cos(Sub_pos_Data{subi}.ori(index)) + (Sub_pos_Data{subi}.z(index)-3).*sin(Sub_pos_Data{subi}.ori(index))...
            ,(Sub_pos_Data{subi}.z(index)-3).*cos(Sub_pos_Data{subi}.ori(index))-Sub_pos_Data{subi}.x(index).*sin(Sub_pos_Data{subi}.ori(index))+3];
        scatter(AllData{subi,Scenei+1}(:,1),AllData{subi,Scenei+1}(:,2),100,'.k');hold on
    end
    scatter(VH_pos(1,:),VH_pos(2,:),400,'.','MarkerEdgeColor', [.75 0 0],'MarkerFaceColor', [.75 0 0]);hold on;
    quiver(VH_pos(1,:) ,VH_pos(2,:),cos(VH_pos(3,:)),sin(VH_pos(3,:)) ,.3,'.','color', [.75 0 0]);hold on;
    axis equal;box off;grid off;caxis([-1.1,-0.1]);
    xlabel('Distance(m)');ylabel('Distance(m)');
    a=colorbar;a.TickLabels = a.Ticks+1.1;
    WinSize = get(0);
    set(gcf, 'position', WinSize.ScreenSize.*0.8);ylim([ymin+2 ymax-2]);xlim([xmin+2 xmax-2]);
end

%%
% for Scenei = 0:4
xmin = -4;xmax = 4;ymin = -1;ymax = 7;Dis_step = 0.1;Ori_step = pi/180;

Scenei = 2;
VH_pos = [Scene_SubData{Scenei+1}.x', Scene_SubData{Scenei+1}.y', Scene_SubData{Scenei+1}.o'];
figure
[U, Ori, Ori_record] = CalculatePotentialField([xmin,xmax;ymin,ymax],Dis_step,Ori_step,VH_pos,0,1);
surf(xmin:Dis_step:xmax,ymin:Dis_step:ymax,U');hold on;view(0,90);shading interp;
axis equal;ylim([ymin ymax]);xlim([xmin xmax]);box off;grid off;
colormap jet;


AllData{subi,Scenei+1}(1,1), AllData{subi,Scenei+1}(1,2)
pos_index = [round((temp(:,1)-xmin)./Dis_step+1)  round((temp(:,2)-ymin)./Dis_step+1)];
ix_seq = xmin:Dis_step:xmax;
iy_seq = ymin:Dis_step:ymax;
temp_Rho = squeeze(Ori_record(round((AllData{subi,Scenei+1}(1,1)-xmin)/Dis_step+1),round(((AllData{subi,Scenei+1}(1,2))-ymin)/Dis_step+1),:));

for ix = 1:length(ix_seq)
    for iy = 1:length(iy_seq)
        [Theta,Rho]= cart2pol((ix_seq(ix)-AllData{subi,Scenei+1}(1,1)),(iy_seq(iy)-AllData{subi,Scenei+1}(1,2)));
        ori_strength(ix,iy)=temp_Rho( round((Theta)/(Ori_step) + pi/Ori_step + 1) );
    end
end

surf(ix_seq,iy_seq,ori_strength(:,:)');hold on;shading interp
% scatter(AllData{subi,Scenei+1}(1,1),AllData{subi,Scenei+1}(1,2),50,'.k');hold on
% plot([AllData{subi,Scenei+1}(1,1)+5*cos(nanmean(wrapToPi(AllData_o{subi,Scenei+1}(1)))), AllData{subi,Scenei+1}(1,1)],[AllData{subi,Scenei+1}(1,2)+5*sin(nanmean(wrapToPi(AllData_o{subi,Scenei+1}(1)))),AllData{subi,Scenei+1}(1,2)],'--k','linewidth',2);
grid off;axis equal;view(0,90);ylim([ymin ymax]);xlim([xmin xmax]);
caxis([0 max(ori_strength,[],'all')]);
colormap jet;