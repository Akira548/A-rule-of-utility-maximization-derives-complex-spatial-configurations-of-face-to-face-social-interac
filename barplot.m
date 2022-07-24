tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));
addpath(genpath('./Functions'));addpath(genpath('./Data'));
load('.\Data\exp3a.mat');
clear odis orid exp3a Utility;
xmax = 6; xmin = -6; ymax = 6;ymin = -6; Dis_step = 0.1;Ori_step = pi/36;

Contrast = CalculatePotentialContrast([xmin,xmax;ymin,ymax],Dis_step,Ori_step,[0, 0, 0]);
for subi = 1:size(SubData.Distance,1)
    pos_I_all(subi,1,:) = [mean(SubData.z(subi,(conditon_d(subi,:) == 3))),mean(SubData.x(subi,(conditon_d(subi,:) == 3)))];
    pos_I_all(subi,2,:) = [mean(SubData.z(subi,(conditon_d(subi,:) == 5))),mean(SubData.x(subi,(conditon_d(subi,:) == 5)))];
    pos_I_all(subi,3,:) = [mean(SubData.z(subi,(conditon_d(subi,:) == 7))),mean(SubData.x(subi,(conditon_d(subi,:) == 7)))];
    ori_I_all(subi,1,:) = mean([cos(pi-(((wrapToPi(deg2rad(SubData.o(subi,(conditon_d(subi,:) == 3)))) ))))', sin(pi-(((wrapToPi(deg2rad(SubData.o(subi,(conditon_d(subi,:) == 3)))) ))))']);
    ori_I_all(subi,2,:) = mean([cos(pi-(((wrapToPi(deg2rad(SubData.o(subi,(conditon_d(subi,:) == 5)))) ))))', sin(pi-(((wrapToPi(deg2rad(SubData.o(subi,(conditon_d(subi,:) == 5)))) ))))']);
    ori_I_all(subi,3,:) = mean([cos(pi-(((wrapToPi(deg2rad(SubData.o(subi,(conditon_d(subi,:) == 7)))) ))))', sin(pi-(((wrapToPi(deg2rad(SubData.o(subi,(conditon_d(subi,:) == 7)))) ))))']);
end
pos_I_x= pos_I_all(:,:,1);
pos_I_y= pos_I_all(:,:,2);
ori_I_x= ori_I_all(:,:,1);
ori_I_y= ori_I_all(:,:,2);

for i = 1:length(pos_I_x(:))
    exp3a.pos(i,:) = [pos_I_x(i),pos_I_y(i)];
    v = [ori_I_x(i) ori_I_y(i)];
    [exp3a.ori(i),~] = cart2pol(v(1),v(2));
    exp3a.odis(i,:) = ExpContrast_o_Distance(Contrast,exp3a.pos(i,:),[xmin,xmax;ymin,ymax],Dis_step,[0, 0, 0]);
    exp3a.orid(i,:) = ExpContrast_Angle(Contrast,exp3a.pos(i,:),exp3a.ori(i),[xmin,xmax;ymin,ymax],Dis_step,[0, 0, 0]);
end

disp([mean(exp3a.odis,1); mean(exp3a.orid,1)]);
exp3a.Utility{1} = Contrast;

%%
tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));
addpath(genpath('./Functions'));addpath(genpath('./Data'));
load('.\Data\exp3b.mat')
xmax = 6; xmin = -6; ymax = 6;ymin = -6; Dis_step = 0.1;Ori_step = pi/36;
VH_pos = [0,3];
clear odis orid trials exp3b Utility;
trials = 1;
for orii = 1:length(oricondition)
    clear tempI tempv pos_I_x pos_I_y ori_I_x ori_I_y ;
    Contrast = CalculatePotentialContrast([xmin,xmax;ymin,ymax],Dis_step,Ori_step,[0, 3,wrapToPi(-pi/2-deg2rad(oricondition(orii)))]);
    exp3b.Utility{orii} = Contrast;
    for subi = 1:size(FilesName,1)
        index = find((Sub_pos_Data{subi}.ori) == oricondition(orii));
        tempv(subi,:,:) = mean([cos(Sub_pos_Data{subi}.heading(index)) sin(Sub_pos_Data{subi}.heading(index))]);
        tempI(subi,:,:) = mean([Sub_pos_Data{subi}.x(index) Sub_pos_Data{subi}.z(index)]);
    end
    pos_I_x = tempI(:,:,1);
    pos_I_y = tempI(:,:,2);
    ori_I_x = tempv(:,:,1);
    ori_I_y = tempv(:,:,2);
    
    for i = 1:length(pos_I_x(:))
        exp3b.pos(trials,:) = [pos_I_x(i) pos_I_y(i)];
        [exp3b.ori(trials),~]=  cart2pol(ori_I_x(i),ori_I_y(i));
        exp3b.odis(trials,:) = ExpContrast_o_Distance(Contrast,exp3b.pos(trials,:),[xmin,xmax;ymin,ymax],Dis_step,[0, 3,wrapToPi(-pi/2-deg2rad(oricondition(orii)))]);
        exp3b.orid(trials,:) = ExpContrast_Angle(Contrast,exp3b.pos(trials,:),exp3b.ori(trials),[xmin,xmax;ymin,ymax],Dis_step,[0, 3,wrapToPi(-pi/2-deg2rad(oricondition(orii)))]);
        trials = trials + 1;
    end
end
disp([nanmean(exp3b.odis,1); nanmean(exp3b.orid,1)]);
% %%
% [exp3ab.Shuffled_dis,exp3ab.Shuffled_Angle,exp3ab.dis_CI,exp3ab.Angle_CI] = Shuffled_level(exp3b.Utility,[[exp3a.pos(:,1), exp3a.pos(:,2)+3];exp3b.pos],[exp3a.ori, exp3b.ori],[-3,3;-2,6],0.1);
% disp([mean(exp3ab.Shuffled_dis,2); mean(exp3ab.Shuffled_Angle,2)]);
% 
% figure
% bar_odis(1,1)=mean(exp3a.odis(:,1),'all');
% bar_odis(2,1)=nanmean(exp3b.odis(:,1),'all');
% bar_odis(3,1)=nanmean(exp3ab.Shuffled_dis,'all');
% bar_orid(1,1)=mean(exp3a.orid(:,1),'all');
% bar_orid(2,1)=nanmean(exp3b.orid(:,1),'all');
% bar_orid(3,1)=nanmean(exp3ab.Shuffled_Angle,'all');
% 
% [~, exp3a.odis_CI]=MD_CI(exp3a.odis(:,1) ,0);
% [~, exp3a.Angle_CI]=MD_CI(exp3b.orid(:,1) ,0);
% [~, exp3b.odis_CI]=MD_CI(exp3b.odis(:,1) ,0);
% [~, exp3b.Angle_CI]=MD_CI(exp3b.orid(:,1) ,0);
% 
% dCI(1) = abs(exp3a.odis_CI(1)-exp3a.odis_CI(2))/2;
% dCI(2) = abs(exp3b.odis_CI(1)-exp3b.odis_CI(2))/2;
% dCI(3) = abs(exp3ab.dis_CI(1)-exp3ab.dis_CI(2))/2;
% oCI(1) = abs(exp3a.Angle_CI(1)-exp3a.Angle_CI(2))/2;
% oCI(2) = abs(exp3b.Angle_CI(1)-exp3b.Angle_CI(2))/2;
% oCI(3) = abs(exp3ab.Angle_CI(1)-exp3ab.Angle_CI(2))/2;
% 
% figure(1)
% h=bar(1:size(bar_odis,1),bar_odis);hold on;
% errorbar((1:size(bar_odis,1)),bar_odis(:,1)',dCI,'k','LineStyle','none','Marker','none');hold on;
% box off;set(gcf,'Position',[0,0,800,400]);
% 
% figure(2)
% h=bar(1:size(bar_orid,1),bar_orid);hold on;
% errorbar((1:size(bar_orid,1)),bar_orid(:,1)',oCI,'k','LineStyle','none','Marker','none');hold on;
% box off;set(gcf,'Position',[0,0,800,400]);

%%
tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));
addpath(genpath('./Functions'));addpath(genpath('./Data'));
load('.\Data\exp5a.mat')
xmax = 6; xmin = -6; ymax = 6;ymin = -6; Dis_step = 0.1;Ori_step = pi/36;
clear odis orid trials Utility exp5a dCI oCI;
close all;clc;counter = 1;

[condi_d1, condi_d2, condi_o1, condi_o2]= ndgrid(DistanceA, -DistanceB, wrapToPi(deg2rad(-ORientaitonA-90)), wrapToPi(deg2rad(ORientaitonB-90)));
condi_d1 = condi_d1(:);
condi_d2 = condi_d2(:);
condi_o1 = condi_o1(:);
condi_o2 = condi_o2(:);

VH_x= [condi_d1,condi_d2];
VH_y= [3.5.*ones(size(condi_d1)), 3.5.*ones(size(condi_d1))];
VH_o = [condi_o1,condi_o2];

for Obi = 1:size(ORientaitonB,2)
    for Oai = 1:size(ORientaitonA,2)
        for Dbi = 1:size(DistanceB,2)
            for Dai = 1:size(DistanceA,2)
                index =  find(all_DistanceA == DistanceA(Dai) &all_DistanceB == DistanceB(Dbi)...
                    &all_OriA == ORientaitonA(Oai)&all_OriB == ORientaitonB(Obi));
                sub_num = unique(all_SubNum);
                for subi = 1:size(SubData,2)
                    index_sub =  find(all_DistanceA == DistanceA(Dai) &all_DistanceB == DistanceB(Dbi)...
                        &all_OriA == ORientaitonA(Oai)&all_OriB == ORientaitonB(Obi)& all_SubNum == SubNum(subi));
                    all_x_sub(counter,subi) = nanmean(all_x(index_sub));
                    all_y_sub(counter,subi) = nanmean(all_y(index_sub));
                    exp5a.pos(counter,subi,:)=[nanmean(all_x(index_sub)),nanmean(all_y(index_sub))];
                    all_o_v_sub(counter,subi,:) = nanmean([cos(deg2rad((90 - all_r(index_sub)))), sin((deg2rad(90 - all_r(index_sub)) ))]);
                    exp5a.ori(counter,subi) = cart2pol( all_o_v_sub(counter,subi,1) , all_o_v_sub(counter,subi,2));
                end
                counter = counter+1;
            end
        end
    end
end
for condi = 1:size(VH_x,1)
    Contrast = CalculatePotentialContrast([xmin,xmax;ymin,ymax],Dis_step,Ori_step,[VH_x(condi,:)', VH_y(condi,:)', VH_o(condi,:)']);
    exp5a.Utility{condi} = Contrast;
    for subi = 1:size(SubData,2)
        exp5a.odis(subi,condi,:) = ExpContrast_o_Distance(Contrast,[all_x_sub(condi,subi)', all_y_sub(condi,subi)'],[xmin,xmax;ymin,ymax],Dis_step,[VH_x(condi,:)', VH_y(condi,:)', VH_o(condi,:)']);
        exp5a.orid(subi,condi,:) = ExpContrast_Angle(Contrast,[all_x_sub(condi,subi)', all_y_sub(condi,subi)'],exp5a.ori(condi,subi)',[xmin,xmax;ymin,ymax],Dis_step,[VH_x(condi,:)', VH_y(condi,:)', VH_o(condi,:)']);
    end
    clc;disp(['----- ' num2str(round(condi./size(VH_x,1)*100)) '% -----']);
end
disp([nanmean(exp5a.odis,'all'); nanmean(exp5a.orid,'all')]);
disp(['done!   ' datestr(now,31)]);

MD_CI(exp5a.odis(:,:,1));
[tbl2,stats2]=simple_mixed_anova(exp5a.odis(:,:,1));
partial_Eta_squared = ES_Eta_squared(tbl2,2);

% deviation distance analysis
MD_CI(exp5a.orid(:,:,1));
[tbl3,stats3]=simple_mixed_anova(exp5a.orid(:,:,1));
partial_Eta_squared = ES_Eta_squared(tbl3,2);

% %%
% tempx = exp5a.pos(:,:,1);tempy = exp5a.pos(:,:,2);
% [exp5a.Shuffled_dis,exp5a.Shuffled_Angle,exp5a.Shuffled_dis_CI,exp5a.Shuffled_Angle_CI] = Shuffled_level(exp5a.Utility,[tempx(:), tempy(:)],exp5a.ori(:),[xmin,xmax;ymin,ymax],Dis_step);
% disp([mean(exp5a.Shuffled_dis,2); mean(exp5a.Shuffled_Angle,2)]);
% 
% figure
% exp5a.bar_odis(1,1)=nanmean(exp5a.odis(:,:,1),'all');
% exp5a.bar_odis(2,1)=nanmean(exp5a.Shuffled_dis,'all');
% exp5a.bar_orid(1,1)=mean(exp5a.orid(:,:,1),'all');
% exp5a.bar_orid(2,1)=nanmean(exp5a.Shuffled_Angle,'all');
% 
% [~, exp5a.odis_CI]=MD_CI(exp5a.odis(:,:,1) ,0);
% [~, exp5a.Angle_CI]=MD_CI(exp5a.orid(:,:,1) ,0);
% 
% dCI(1) = abs(exp5a.odis_CI(1)-exp5a.odis_CI(2))/2;
% oCI(1) = abs(exp5a.Angle_CI(1)-exp5a.Angle_CI(2))/2;
% dCI(2) = abs(exp5a.Shuffled_dis_CI(1)-exp5a.Shuffled_dis_CI(2))/2;
% oCI(2) = abs(exp5a.Shuffled_Angle_CI(1)-exp5a.Shuffled_Angle_CI(2))/2;
% 
% figure(1)
% h=bar(1:size(exp5a.bar_odis,1),exp5a.bar_odis);hold on;
% errorbar((1:size(exp5a.bar_odis,1)),exp5a.bar_odis,dCI,'k','LineStyle','none','Marker','none');hold on;
% box off;set(gcf,'Position',[0,0,800,400]);
% 
% figure(2)
% h=bar(1:size(exp5a.bar_orid,1),exp5a.bar_orid');hold on;
% errorbar((1:size(exp5a.bar_orid,1)),exp5a.bar_orid',oCI,'k','LineStyle','none','Marker','none');hold on;
% box off;set(gcf,'Position',[0,0,800,400]);
%%
% clear all;clc;
clear odis orid AllData AllData_o  Utility exp5b dCI oCI;
ExpName = 'Exp5b';
tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));
addpath(genpath('./Functions'));
load('.\Data\exp5b.mat');
xmax = 6; xmin = -6; ymax = 6;ymin = -6; Dis_step = 0.1;Ori_step = pi/36;
for Scenei = 0:4
    VH_pos = [Scene_SubData{Scenei+1}.x', Scene_SubData{Scenei+1}.y', Scene_SubData{Scenei+1}.o'];
    Contrast = CalculatePotentialContrast([xmin,xmax;ymin,ymax],Dis_step,Ori_step,VH_pos);
    exp5b.Utility{Scenei+1} = Contrast;
    for subi = 1:size(Sub_pos_Data,2)
        temp_index = find(Sub_pos_Data{subi}.Scenes == Scenei);
        [~,sortindex]= (sort(Sub_pos_Data{subi}.ori(temp_index)));
        index = temp_index(sortindex);
        exp5b.pos(subi,Scenei+1,:,:) = [Sub_pos_Data{subi}.x(index).*cos(Sub_pos_Data{subi}.ori(index)) + (Sub_pos_Data{subi}.z(index)-3).*sin(Sub_pos_Data{subi}.ori(index))...
            ,(Sub_pos_Data{subi}.z(index)-3).*cos(Sub_pos_Data{subi}.ori(index))-Sub_pos_Data{subi}.x(index).*sin(Sub_pos_Data{subi}.ori(index))+3];
        exp5b.ori(subi,Scenei+1,:) = wrapToPi(pi/2-Sub_pos_Data{subi}.ori(index) - deg2rad(Sub_pos_Data{subi}.heading(index)));
        for triali = 1:size(exp5b.pos,3)
            if exp5b.pos(subi,Scenei+1,triali,2)<=6
                exp5b.odis(subi,Scenei+1,triali,:) = ExpContrast_o_Distance(Contrast,[exp5b.pos(subi,Scenei+1,triali,1) exp5b.pos(subi,Scenei+1,triali,2)],[xmin,xmax;ymin,ymax],Dis_step,VH_pos);
                exp5b.orid(subi,Scenei+1,triali,:) = ExpContrast_Angle(Contrast,[exp5b.pos(subi,Scenei+1,triali,1) exp5b.pos(subi,Scenei+1,triali,2)],exp5b.ori(subi,Scenei+1,triali),[xmin,xmax;ymin,ymax],Dis_step,VH_pos);
            else
                exp5b.odis(subi,Scenei+1,triali,:) = nan;
                exp5b.orid(subi,Scenei+1,triali,:) = nan;
            end
        end
    end
end
disp([squeeze(nanmean(exp5b.odis,1:3))'; squeeze(nanmean(exp5b.orid,1:3))']);
disp(['done!   ' datestr(now,31)]);

MD_CI(exp5b.odis(:,:,:,1));
[tbl2,stats2]=simple_mixed_anova(exp5b.odis(:,:,:,1));
partial_Eta_squared = ES_Eta_squared(tbl2,2);

% deviation distance analysis
MD_CI(exp5b.orid(:,:,:,1));
[tbl3,stats3]=simple_mixed_anova(exp5b.orid(:,:,:,1));
partial_Eta_squared = ES_Eta_squared(tbl3,2);
% %%
% tempx = exp5b.pos(:,:,:,1);tempy = exp5b.pos(:,:,:,2);tempo = exp5b.ori(:,:,:);
% [exp5b.Shuffled_dis,exp5b.Shuffled_Angle,exp5b.Shuffled_dis_CI,exp5b.Shuffled_Angle_CI] = Shuffled_level(exp5b.Utility,[tempx(:), tempy(:)],tempo(:),[xmin,xmax;ymin,ymax],Dis_step);
% disp([mean(exp5b.Shuffled_dis,2); mean(exp5b.Shuffled_Angle,2)]);
% 
% figure
% exp5b.bar_odis(1,1)=nanmean(exp5b.odis(:,:,:,1),'all');
% exp5b.bar_odis(2,1)=nanmean(exp5b.Shuffled_dis,'all');
% exp5b.bar_orid(1,1)=mean(exp5b.orid(:,:,:,1),'all');
% exp5b.bar_orid(2,1)=nanmean(exp5b.Shuffled_Angle,'all');
% 
% [~, exp5b.odis_CI]=MD_CI(exp5b.odis(:,:,:,1) ,0);
% [~, exp5b.Angle_CI]=MD_CI(exp5b.orid(:,:,:,1) ,0);
% 
% dCI(1) = abs(exp5b.odis_CI(1)-exp5b.odis_CI(2))/2;
% oCI(1) = abs(exp5b.Angle_CI(1)-exp5b.Angle_CI(2))/2;
% dCI(2) = abs(exp5b.Shuffled_dis_CI(1)-exp5b.Shuffled_dis_CI(2))/2;
% oCI(2) = abs(exp5b.Shuffled_Angle_CI(1)-exp5b.Shuffled_Angle_CI(2))/2;
% 
% figure(1)
% h=bar(1:size(exp5b.bar_odis,1),exp5b.bar_odis);hold on;
% errorbar((1:size(exp5b.bar_odis,1)),exp5b.bar_odis,dCI,'k','LineStyle','none','Marker','none');hold on;
% box off;set(gcf,'Position',[0,0,800,400]);
% 
% figure(2)
% h=bar(1:size(exp5b.bar_orid,1),exp5b.bar_orid');hold on;
% errorbar((1:size(exp5b.bar_orid,1)),exp5b.bar_orid',oCI,'k','LineStyle','none','Marker','none');hold on;
% box off;set(gcf,'Position',[0,0,800,400]);
%%
clear odis orid AllData AllData_o;
ExpName = 'Exp6';
tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));
addpath(genpath('./Functions'));
load('.\Data\exp6.mat')
xmax = 6; xmin = -6; ymax = 6;ymin = -6; Dis_step = 0.1;Ori_step = pi/36;
JoinSeqenceName = [2:4,6:8,10:12,14:16];
FormSeqenceName = [23:24,2,5,9,13,21,22,25,26];
JoinTimeLabel = [6,10,8, 6,7,6, 7,8,5, 5,8,7];
JoinWalkerLabel = [4,5,6, 3,2,6, 4,2,1, 6,2,1];

GS=[4 6 4 5 6 4 5 6 4 6];
GroupMember = 1:6;
for Seqence = 1:length(JoinSeqenceName)
    seq = (0:300:JoinTimeLabel(Seqence)*300)+1;
    q = [];x = [];y = [];o = [];
    for member = 1:6
        if ~isempty(SeqenceData{JoinSeqenceName(Seqence),member}) && member ~= JoinWalkerLabel(Seqence)
            x = [x,nanmean(SeqenceData{JoinSeqenceName(Seqence),member}.X(JoinTimeLabel(Seqence)*300+1))];
            y = [y,nanmean(SeqenceData{JoinSeqenceName(Seqence),member}.Z(JoinTimeLabel(Seqence)*300+1))];
            v_of_o = mean([cos(SeqenceData{JoinSeqenceName(Seqence),member}.Rotation_Y(JoinTimeLabel(Seqence)*300+1)),sin(SeqenceData{JoinSeqenceName(Seqence),member}.Rotation_Y(JoinTimeLabel(Seqence)*300+1))],1);
            [mean_o,~]=cart2pol(v_of_o(1),v_of_o(2));
            o = [o,mean_o];
        end
    end
    for member = 1:6
        if member == JoinWalkerLabel(Seqence)
            Testx = SeqenceData{JoinSeqenceName(Seqence),member}.X(JoinTimeLabel(Seqence)*300+1);
            Testy = SeqenceData{JoinSeqenceName(Seqence),member}.Z(JoinTimeLabel(Seqence)*300+1);
            Testo = SeqenceData{JoinSeqenceName(Seqence),member}.Rotation_Y(JoinTimeLabel(Seqence)*300+1);
            Testv = [cos(SeqenceData{JoinSeqenceName(Seqence),member}.Rotation_Y(JoinTimeLabel(Seqence)*300+1)),sin(SeqenceData{JoinSeqenceName(Seqence),member}.Rotation_Y(JoinTimeLabel(Seqence)*300+1))];
            Contrast = CalculatePotentialContrast([xmin,xmax;ymin,ymax],.1,pi/36,[x',y',o']);
            exp6.Utility{Seqence} = Contrast;
            exp6.pos(Seqence,:)=[Testx,Testy];
            exp6.ori(Seqence) = Testo;
            exp6.odis(Seqence,:) = ExpContrast_o_Distance(Contrast,[Testx Testy],[xmin,xmax;ymin,ymax],Dis_step,[x',y',o']);
            exp6.orid(Seqence,:) = ExpContrast_Angle(Contrast,[Testx Testy],Testo,[xmin,xmax;ymin,ymax],Dis_step,[x',y',o']);
        end
    end
end
disp([squeeze(nanmean(exp6.odis,1)); squeeze(nanmean(exp6.orid,1))]);
MD_CI(exp6.odis(:,1));
MD_CI(exp6.orid(:,1));


%%

tempx_exp3a = exp3a.pos(:,1);
tempy_exp3a = exp3a.pos(:,2);
tempo_exp3a = exp3a.ori;

tempx_exp3b = exp3b.pos(:,1);
tempy_exp3b = exp3b.pos(:,2);
tempo_exp3b = exp3b.ori;

tempx_exp5a = exp5a.pos(:,:,1);
tempy_exp5a = exp5a.pos(:,:,2);
tempo_exp5a = exp5a.ori(:,:);

tempx_exp5b = exp5b.pos(:,:,:,1);
tempy_exp5b = exp5b.pos(:,:,:,2);
tempo_exp5b = exp5b.ori(:,:,:); 
tempx_exp5b2 = tempy_exp5b(:);% rule out one bad trial
tempy_exp5b2 = tempy_exp5b(:);
tempo_exp5b2 = tempo_exp5b(:);
tempx_exp5b = tempx_exp5b2(tempy_exp5b<=6);
tempy_exp5b = tempy_exp5b2(tempy_exp5b<=6);
tempo_exp5b = tempo_exp5b2(tempy_exp5b<=6);

tempx_exp6 = exp6.pos(:,1);
tempy_exp6 = exp6.pos(:,2);
tempo_exp6 = exp6.ori;

trials_x = [tempx_exp3a(:); tempx_exp3b(:); tempx_exp5a(:); tempx_exp5b(:); tempx_exp6(:)];
trials_y = [tempy_exp3a(:); tempy_exp3b(:); tempy_exp5a(:); tempy_exp5b(:); tempy_exp6(:)];
trials_o = [tempo_exp3a(:); tempo_exp3b(:); tempo_exp5a(:); tempo_exp5b(:); tempo_exp6(:)];

Utility = [exp3a.Utility,exp3b.Utility,exp5a.Utility,exp5b.Utility,exp6.Utility];

[Shuffled_dis,Shuffled_Angle,Shuffled_dis_CI,Shuffled_Angle_CI] = Shuffled_level(Utility,[trials_x(:), trials_y],trials_o,[xmin,xmax;ymin,ymax],Dis_step);
disp([mean(Shuffled_dis,2); mean(Shuffled_Angle,2)]);

%%
bar_odis(1,:)=mean(exp3a.odis,1);
bar_odis(2,:)=nanmean(exp3b.odis,1);
bar_odis(3,:)=nanmean(exp5a.odis,1:2);
bar_odis(4,:)=squeeze(nanmean(exp5b.odis,1:3));
bar_odis(5,:)=nanmean(exp6.odis,1);

[~, bar_odis_CI(1,1,:)]=MD_CI(exp3a.odis(:,1) ,0);
[~, bar_odis_CI(2,1,:)]=MD_CI(exp3b.odis(:,1) ,0);
[~, bar_odis_CI(3,1,:)]=MD_CI(exp5a.odis(:,:,1) ,0);
[~, bar_odis_CI(4,1,:)]=MD_CI(exp5b.odis(:,:,:,1) ,0);
[~, bar_odis_CI(5,1,:)]=MD_CI(exp6.odis(:,1) ,0);

bar_orid(1,:)=mean(exp3a.orid,1);
bar_orid(2,:)=nanmean(exp3b.orid,1);
bar_orid(3,:)=nanmean(exp5a.orid,1:2);
bar_orid(4,:)=squeeze(nanmean(exp5b.orid,1:3));
bar_orid(5,:)=nanmean(exp6.orid,1);

[~, bar_orid_CI(1,:)]=MD_CI(exp3a.orid(:,1) ,0);
[~, bar_orid_CI(2,:)]=MD_CI(exp3b.orid(:,1) ,0);
[~, bar_orid_CI(3,:)]=MD_CI(exp5a.orid(:,:,1) ,0);
[~, bar_orid_CI(4,:)]=MD_CI(exp5b.orid(:,:,:,1) ,0);
[~, bar_orid_CI(5,:)]=MD_CI(exp6.orid(:,1) ,0);

figure(1)
h=bar(1:size(bar_odis,1),bar_odis(:,1));hold on;
errorbar((1:size(bar_odis,1)),bar_odis(:,1),abs(bar_odis_CI(:,1)'-bar_odis_CI(:,2)')/2,'k','LineStyle','none','Marker','none');hold on;
plot([0 size(bar_odis,1)+1],[mean(Shuffled_dis,2) mean(Shuffled_dis,2)],'--k');
box off;set(gcf,'Position',[0,0,800,400]);

figure(2)
h=bar(1:size(bar_orid,1),bar_orid(:,1));hold on;
errorbar((1:size(bar_orid,1)),bar_orid(:,1),abs(bar_orid_CI(:,1)'-bar_orid_CI(:,2)')/2,'k','LineStyle','none','Marker','none');hold on;
plot([0 size(bar_orid,1)+1],[mean(Shuffled_Angle,2) mean(Shuffled_Angle,2)],'--k');
box off;set(gcf,'Position',[0,0,800,400]);