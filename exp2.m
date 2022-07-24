%% Exp2a
clc;close all;clear all;
ExpName = 'Exp2';
tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));
addpath(genpath('./Functions'));

FilesName = dir(['.\Raw Data\' ExpName filesep '*.txt']);
for i = 1:size(FilesName,1)
    tempFilename{i} = FilesName(i).name(1:end-6);
end
newFilename=unique(tempFilename);
BlockFilesName = dir(['.\Raw Data\' ExpName filesep,char(newFilename(1)),'_*.txt']);
temp = importfile_Exp1(['.\Raw Data\' ExpName filesep,BlockFilesName(1).name]);
temp.VisualAngle = temp.VisualAngle.*(-temp.VisualAngleDirection+0.5)*2;
Condition_Distance = unique(temp.Distance);
Condition_Theta1 = unique(temp.VisualAngle);
Condition_Theta2 = unique(temp.A_Angle);
clear temp
for subi = 1:size(newFilename,2)
    BlockFilesName = dir(['.\Raw Data\' ExpName filesep,char(newFilename(subi)),'_*.txt']);
    disp(['processing: sub' num2str(subi)])
    index_count = 1;
    for blocki = 1:size(BlockFilesName,1)
        if blocki == 1
            temp = importfile_Exp1(['.\Raw Data\' ExpName filesep,BlockFilesName(blocki).name]);
            temp.VisualAngle = temp.VisualAngle.*(-temp.VisualAngleDirection+0.5)*2;
        else
            temp2 = importfile_Exp1(['.\Raw Data\' ExpName filesep,BlockFilesName(blocki).name]);
            temp2.VisualAngle = temp2.VisualAngle.*(-temp2.VisualAngleDirection+0.5)*2;
            temp = [temp;temp2];
        end
    end
    temp_all{subi} = temp;
    for  Distance_i = 1:length(Condition_Distance)
        for  Theta1_i = 1:length(Condition_Theta1)
            index = find(temp.Distance == Condition_Distance(Distance_i) & temp.VisualAngle == Condition_Theta1(Theta1_i));
            if isempty(index)
                SubAll_a{subi}.Distance = [SubAll_a{subi}.Distance Condition_Distance(Distance_i)];
                SubAll_a{subi}.VisualAngle = [SubAll_a{subi}.VisualAngle Condition_Theta1(Theta1_i)];
                SubAll_a{subi}.Answer = [SubAll_a{subi}.Answer nan];
                SubAll_a{subi}.Sub = [SubAll_a{subi}.Sub; subi];
            else
                if (index_count == 1)
                    SubAll_a{subi}.Distance = Condition_Distance(Distance_i);
                    SubAll_a{subi}.VisualAngle = Condition_Theta1(Theta1_i);
                    SubAll_a{subi}.Answer = nanmean(temp.Answer(index));
                    SubAll_a{subi}.Sub = subi;
                else
                    SubAll_a{subi}.Distance = [SubAll_a{subi}.Distance Condition_Distance(Distance_i)];
                    SubAll_a{subi}.VisualAngle = [SubAll_a{subi}.VisualAngle Condition_Theta1(Theta1_i)];
                    SubAll_a{subi}.Answer = [SubAll_a{subi}.Answer nanmean(temp.Answer(index))];
                    SubAll_a{subi}.Sub = [SubAll_a{subi}.Sub; subi];
                end
            end
            index_count = index_count+1;
        end
    end
end
for subi = 1:size(newFilename,2)
    exp2a.Distance(:,subi) = SubAll_a{subi}.Distance';
    exp2a.VisualAngle(:,subi) = SubAll_a{subi}.VisualAngle';
    exp2a.Answer(:,subi) = SubAll_a{subi}.Answer';
    exp2a.Sub(:,subi) = SubAll_a{subi}.Sub';
end
%%
close all
dis_condition = round(unique(exp2a.Distance),1);
c = gray(length(dis_condition)+2);
temp_va=mean(exp2a.VisualAngle,2);
temp_dis=round(mean(exp2a.Distance,2),1);
temp_p=mean(exp2a.Answer,2);
temp_std=nanstd(exp2a.Answer,[],2);
for disi = 1:length(dis_condition)
    b(disi) = plot(temp_va(temp_dis == dis_condition(disi)),temp_p(temp_dis == dis_condition(disi)),'linewidth',1,'Color',c(disi,:));hold on
    errorbar(temp_va(temp_dis == dis_condition(disi)),...
         temp_p(temp_dis == dis_condition(disi)),temp_std(temp_dis == dis_condition(disi))./sqrt(size(exp2a.Answer,2)),...
        '.','MarkerSize',1,'CapSize',10,'linewidth',1,'Color',c(disi,:));hold on;
    ylabel('Probability');xlabel('ж╚(бу)');
    set(gca,'XTick',-60:15:60);xlim([-63 63]);ylim([0 1]);box off;
end
set(gcf,'position',[200,200,600,550]);

figure
[Theta,Rho,p] = griddata(deg2rad(temp_va),temp_dis,temp_p,linspace(-pi/3,pi/3)',linspace(0,max(dis_condition)),'linear');
[x,y,z] = pol2cart(Theta,Rho,p); surf(x,y,z);hold on
[Theta2,Rho2,p2] = griddata(deg2rad(temp_va),temp_dis,temp_p,linspace(-pi/3,pi/3)',linspace(0,min(dis_condition)+0.1),'V4');
[x2,y2,z2] = pol2cart(Theta2,Rho2,p2);surf(x2,y2,z2);hold on;
axis equal;view(90,-90);shading interp;caxis([0 1]);
for rhoi = -pi/6:deg2rad(30):pi/3
    [texty textx] = pol2cart(rhoi,10+.05);
    text(texty+0.05,textx-0.04,[num2str(abs(round(rad2deg(rhoi),0))) 'бу'],'FontSize',10);
end
box off;grid off;axis off;
set(gca,'YTick',-1.2:0.4:1.2);set(gca,'YTicklabel',{'-1.2','-0.8','-0.4','0','0.4','0.8','1.2(m)'})
set(gca,'XTick',0:0.4:1.2);set(gca,'XTicklabel',{'0','0.4','0.8','1.2(m)'})

figure
[Theta,Rho,p] = griddata(deg2rad(temp_va),temp_dis,Interaction_Function(deg2rad(temp_va),[],temp_dis,'Prep'),linspace(-pi/3,pi/3)',linspace(0,max(dis_condition)),'linear');
[x,y,z] = pol2cart(Theta,Rho,p); surf(x,y,z);hold on
[Theta2,Rho2,p2] = griddata(deg2rad(temp_va),temp_dis,Interaction_Function(deg2rad(temp_va),[],temp_dis,'Prep'),linspace(-pi/3,pi/3)',linspace(0,min(dis_condition)+0.1),'V4');
[x2,y2,z2] = pol2cart(Theta2,Rho2,p2);surf(x2,y2,z2);hold on;
axis equal;view(90,-90);shading interp;caxis([0 1]);
box off;grid off;
axis off;

%% Exp2b

FilesName = dir(['.\Raw Data\' ExpName filesep '*.txt']);
for i = 1:size(FilesName,1)
    tempFilename{i} = FilesName(i).name(1:end-6);
end
newFilename=unique(tempFilename);
BlockFilesName = dir(['.\Raw Data\' ExpName filesep,char(newFilename(1)),'_*.txt']);
temp = importfile_Exp1(['.\Raw Data\' ExpName filesep,BlockFilesName(1).name]);
temp.A_Angle = temp.A_Angle.*(-temp.VisualAngleDirection+0.5)*2;
Condition_Distance = unique(temp.Distance);
Condition_Theta1 = unique(temp.VisualAngle);
Condition_Theta2 = unique(temp.A_Angle);

for subi = 1:size(newFilename,2)
    BlockFilesName = dir(['.\Raw Data\' ExpName filesep,char(newFilename(subi)),'_*.txt']);
    disp(['processing: sub' num2str(subi)])
    index_count = 1;
    for blocki = 1:size(BlockFilesName,1)
        if blocki == 1
            temp = importfile_Exp1(['.\Raw Data\' ExpName filesep,BlockFilesName(blocki).name]);
            temp.A_Angle = temp.A_Angle.*(-temp.VisualAngleDirection+0.5)*2;
        else
            temp2 = importfile_Exp1(['.\Raw Data\' ExpName filesep,BlockFilesName(blocki).name]);
            temp2.A_Angle = temp2.A_Angle.*(-temp2.VisualAngleDirection+0.5)*2;
            temp = [temp;temp2];
        end
    end
    for  Distance_i = 1:length(Condition_Distance)
        for  Theta2_i = 1:length(Condition_Theta2)
            index = find(temp.Distance == Condition_Distance(Distance_i) & temp.A_Angle == Condition_Theta2(Theta2_i));
            
            if isempty(index)
                SubAll_b{subi}.Distance = [SubAll_b{subi}.Distance Condition_Distance(Distance_i)];
                SubAll_b{subi}.A_Angle = [SubAll_b{subi}.A_Angle Condition_Theta2(Theta2_i)];
                SubAll_b{subi}.Answer = [SubAll_b{subi}.Answer nan];
                SubAll_b{subi}.Sub = [SubAll_b{subi}.Sub; subi];
            else
                if (index_count == 1)
                    SubAll_b{subi}.Distance = Condition_Distance(Distance_i);
                    SubAll_b{subi}.A_Angle = Condition_Theta2(Theta2_i);
                    SubAll_b{subi}.Answer = nanmean(temp.Answer(index));
                    SubAll_b{subi}.Sub = subi;
                else
                    SubAll_b{subi}.Distance = [SubAll_b{subi}.Distance Condition_Distance(Distance_i)];
                    SubAll_b{subi}.A_Angle = [SubAll_b{subi}.A_Angle Condition_Theta2(Theta2_i)];
                    SubAll_b{subi}.Answer = [SubAll_b{subi}.Answer nanmean(temp.Answer(index))];
                    SubAll_b{subi}.Sub = [SubAll_b{subi}.Sub; subi];
                end
            end
            index_count = index_count+1;
        end
    end
end
for subi = 1:size(newFilename,2)
    exp2b.Distance(:,subi) = SubAll_b{subi}.Distance';
    exp2b.A_Angle(:,subi) = SubAll_b{subi}.A_Angle';
    exp2b.Answer(:,subi) = SubAll_b{subi}.Answer';
    exp2b.Sub(:,subi) = SubAll_b{subi}.Sub';
end

%%

for subi = 1:size(newFilename,2)
    Sub_Anova.Distance(:,subi) = SubAll_b{subi}.Distance';
    Sub_Anova.A_Angle(:,subi) = SubAll_b{subi}.A_Angle';
    Sub_Anova.Answer(:,subi) = SubAll_b{subi}.Answer';
    Sub_Anova.Sub(:,subi) = SubAll_b{subi}.Sub';
end
disp('done!')

%% statistic
clc
o1_cond = unique((Sub_Anova.VisualAngle(:,1)));
o2_cond = unique((Sub_Anova.A_Angle(:,1)));
d_cond = unique((Sub_Anova.Distance(:,1)));
sub_cond = unique((Sub_Anova.Sub));

for subi = 1:size(sub_cond,1)
    for o1i = 1:size(o1_cond,1)
        for o2i = 1:size(o2_cond,1)
            for di = 1:size(d_cond,1)
                exp2(subi,o1i,o2i,di) = Sub_Anova.Answer(Sub_Anova.VisualAngle(:,subi)==o1_cond(o1i) &...
                    Sub_Anova.A_Angle(:,subi)==o2_cond(o2i) & Sub_Anova.Distance(:,subi)==d_cond(di),subi);
            end
        end
    end
end

[tbl,~]=simple_mixed_anova(exp2,[],{'Theta1','Theta2','Distance'});
partial_Eta_squared1 = ES_Eta_squared(tbl,2);

%%
exp2 = table();
exp2.Distance = nanmean(Sub_allmean.Distance,2);
exp2.VisualAngle = (nanmean(Sub_allmean.VisualAngle,2));
exp2.A_Angle = (nanmean(Sub_allmean.A_Angle,2));
exp2.Probability = nanmean(Sub_allmean.Answer,2);

% theta2_condition = [0,15,30,45,60,-15,-30,-45,-60];
theta2_condition = [0,15,-15,30,-30,45,-45,60,-60];

dis_condition = unique(exp2.Distance);
c = gray(10);
for ang_i = 1:9
    subplot(3,3,ang_i)
    temp_va=deg2rad((exp2.A_Angle(find(abs(exp2.VisualAngle == theta2_condition(ang_i))) )));
    temp_dis=exp2.Distance(find(abs(exp2.VisualAngle == theta2_condition(ang_i))));
    temp_p=exp2.Probability(find(abs(exp2.VisualAngle == theta2_condition(ang_i)) ));
    for disi = 1:size(dis_condition,1)
        b(disi) = plot(rad2deg(temp_va(temp_dis == dis_condition(disi))),temp_p(temp_dis == dis_condition(disi)),'linewidth',1,'Color',c(disi,:));hold on
        plot(rad2deg(temp_va(temp_dis == dis_condition(disi))),temp_p(temp_dis == dis_condition(disi)),'o','MarkerSize',5,'linewidth',1,'Color',c(disi,:));hold on
        ylabel('Probability');xlabel('ж╚_2(бу)');
        title(['ж╚_1 = ', num2str(theta2_condition(ang_i)) 'бу']);
        set(gca,'XTick',-150:60:180);xlim([-153 183]);ylim([0 1]);box off;
    end
end
legend(b,{'0.2m','0.3m','0.4m','0.5m','0.6m','0.7m','0.9m','1.2m'},'location','best');
set(gcf,'position',[0,0,1000,800]);

%% figure 1-b
% theta2_condition = [0,15,30,45,60,-15,-30,-45,-60];
theta2_condition = [0,15,-15,30,-30,45,-45,60,-60];

dis_cond = 1.2;warning off

for ang_i = 1:9
    figure(1)
    set(gcf,'position',[0,0,1250,1125]);
    subplot(3,3,ang_i)
    temp_va=deg2rad((exp2.VisualAngle(find(abs(exp2.A_Angle == theta2_condition(ang_i)) & exp2.Distance <= dis_cond) )));
    temp_dis=exp2.Distance(find(abs(exp2.A_Angle == theta2_condition(ang_i)) & exp2.Distance <= dis_cond));
    temp_p=exp2.Probability(find(abs(exp2.A_Angle == theta2_condition(ang_i)) & exp2.Distance <= dis_cond));
    
    [Theta,Rho,p] = griddata(temp_va,temp_dis,temp_p,linspace(-pi/3,pi/3)',linspace(0,dis_cond),'linear');
    [x,y,z] = pol2cart(Theta,Rho,p); surf(x,y,z);hold on
    [Theta2,Rho2,p2] = griddata(temp_va,temp_dis,temp_p,linspace(-pi/3,pi/3)',linspace(0,0.5),'V4');[x2,y2,z2] = pol2cart(Theta2,Rho2,p2);surf(x2,y2,z2);hold on;
    
    % axis
    [Pol_polaraxis, Rho_polaraxis]=meshgrid(-pi/3:deg2rad(.1):pi/3,0:.2:dis_cond);
    [x_polaraxis, y_polaraxis] = pol2cart(Pol_polaraxis,Rho_polaraxis);
    h = scatter3(x_polaraxis(:),y_polaraxis(:),repmat(-5,size(x_polaraxis(:))),.25,'filled','MarkerEdgeColor','none','MarkerFaceColor',[.25 .25 .25]);hold on
    [Pol_polaraxi_line, Rho_polaraxis_line] = meshgrid(-pi/3:deg2rad(30):pi/3,0:0.01:dis_cond);
    [x_polaraxis_line, y_polaraxis_line] = pol2cart(Pol_polaraxi_line,Rho_polaraxis_line);
    h2 = scatter3(x_polaraxis_line(:),y_polaraxis_line(:),repmat(-5,size(x_polaraxis_line(:))),.25,'filled','MarkerEdgeColor','none','MarkerFaceColor',[.25 .25 .25]);hold on
    axis equal;view(90,-90);shading interp;caxis([0 1]);
    %     title({'Angle of View',[]},'FontSize',12);
    box off;grid off;
    set(gca,'YTick',-1.2:0.4:1.2);set(gca,'YTicklabel',{'-1.2','-0.8','-0.4','0','0.4','0.8','1.2(m)'})
    set(gca,'XTick',0:0.4:1.2);set(gca,'XTicklabel',{'0','0.4','0.8','1.2(m)'})
    for rhoi = -pi/6:deg2rad(30):pi/3
        [texty textx] = pol2cart(rhoi,dis_cond+.05);
        text(texty+0.05,textx-0.04,[num2str(abs(round(rad2deg(rhoi),0))) 'бу'],'FontSize',10);
    end
    
    %%%
    figure(2)
    set(gcf,'position',[0,0,1250,1125]);
    subplot(3,3,ang_i)
    [surf_Theta,surf_Rho] = meshgrid(-pi/3:deg2rad(1):pi/3,0:.1:dis_cond);
    [Theta,Rho,p] = griddata(surf_Theta(:),surf_Rho(:),Interaction_Function('R',surf_Theta(:)', repmat(deg2rad(theta2_condition(ang_i)),size(surf_Theta(:)))' ,surf_Rho(:)'),linspace(-pi/3,pi/3)',linspace(0,dis_cond),'linear');
    [x,y,z] = pol2cart(Theta,Rho,p);surf(x,y,z);hold on;
    % axis
    [Pol_polaraxis, Rho_polaraxis]=meshgrid(-pi/3:deg2rad(.1):pi/3,0:.2:dis_cond);
    [x_polaraxis, y_polaraxis] = pol2cart(Pol_polaraxis,Rho_polaraxis);
    h = scatter3(x_polaraxis(:),y_polaraxis(:),repmat(-5,size(x_polaraxis(:))),.25,'filled','MarkerEdgeColor','none','MarkerFaceColor',[.25 .25 .25]);hold on
    [Pol_polaraxi_line, Rho_polaraxis_line] = meshgrid(-pi/3:deg2rad(30):pi/3,0:0.1:dis_cond);
    [x_polaraxis_line, y_polaraxis_line] = pol2cart(Pol_polaraxi_line,Rho_polaraxis_line);
    h2 = scatter3(x_polaraxis_line(:),y_polaraxis_line(:),repmat(-5,size(x_polaraxis_line(:))),.25,'filled','MarkerEdgeColor','none','MarkerFaceColor',[.25 .25 .25]);hold on
    axis equal;view(90,-90);shading interp;caxis([0 1]);grid off
    box off
    set(gca,'YTick',-1.2:0.4:1.2);set(gca,'YTicklabel',{'-1.2','-0.8','-0.4','0','0.4','0.8','1.2(m)'})
    set(gca,'XTick',0:0.4:1.2);set(gca,'XTicklabel',{'0','0.4','0.8','1.2(m)'})
    for rhoi = -pi/6:deg2rad(30):pi/3
        [texty textx] = pol2cart(rhoi,dis_cond+.05);
        text(texty+0.05,textx-0.04,[num2str(abs(round(rad2deg(rhoi),0))) 'бу'],'FontSize',10);
    end
end
%%
warning off
x(:,1) = wrapToPi(deg2rad(Sub_allmean.VisualAngle(:,1)));
x(:,2) =  wrapToPi(deg2rad(Sub_allmean.A_Angle(:,1)));
x(:,3) = (Sub_allmean.Distance(:,1));
y = nanmean(Sub_allmean.Answer,2);

fun_F = @(beta,x) (1-exp(-( beta(4)./  ((beta(1) + abs(x(:,2)).^beta(2) ).*x(:,3) )).^abs(beta(5)) ));

para1 = nlinfit(x,y,fun_F,[randi([0,1]),randi([0,1]),randi([0,1]),randi([0,1]), randi([0,1]),randi([0,1])]);

ft=fitlm(y,fun_F(para1,x));
ft.Rsquared
format short g
disp(para1)
%%
exp2 = table();
exp2.Distance = nanmean(Sub_allmean.Distance,2);
exp2.VisualAngle = (nanmean(Sub_allmean.VisualAngle,2));
exp2.A_Angle = (nanmean(Sub_allmean.A_Angle,2));
exp2.Probability = nanmean(Sub_allmean.Answer,2);

theta2_condition = [0,15,-15,30,-30,45,-45,60,-60];

dis_condition = unique(exp2.Distance);
c = gray(11);

temp_va=deg2rad((exp2.VisualAngle(find(abs(exp2.A_Angle == theta2_condition(1))) )));
temp_dis=exp2.Distance(find(abs(exp2.A_Angle == theta2_condition(1))));
temp_p=exp2.Probability(find(abs(exp2.A_Angle == theta2_condition(1)) ));
for disi = 1:size(dis_condition,1)
    b(disi) = plot(rad2deg(temp_va(temp_dis == dis_condition(disi))),temp_p(temp_dis == dis_condition(disi)),'linewidth',1,'Color',c(disi,:));hold on
    plot(rad2deg(temp_va(temp_dis == dis_condition(disi))),temp_p(temp_dis == dis_condition(disi)),'.','MarkerSize',1,'linewidth',1,'Color',c(disi,:));hold on
    %Errorbar
    temp_indexmat = (round(Sub_allmean.Distance,2) == round(dis_condition(disi),2)&Sub_allmean.A_Angle == theta2_condition(1));
    temp_indexmat2=temp_indexmat;
    temp_indexmat2(all(temp_indexmat2==0,2),:)=[];
    temp_mat = reshape(Sub_allmean.Answer(temp_indexmat),size(temp_indexmat2,1),size(temp_indexmat2,2));
    errorbar(rad2deg(temp_va(temp_dis == dis_condition(disi))),...
        nanmean(temp_mat,2),nanstd(temp_mat,[],2)./sqrt(size(temp_mat,2)),...
        '.','MarkerSize',1,'linewidth',1,'Color',c(disi,:));
    ylabel('Probability');xlabel('ж╚_2(бу)');
    title(['ж╚_1 = ', num2str(theta2_condition(1)) 'бу']);
    set(gca,'XTick',-60:15:60);xlim([-63 63]);ylim([0 1]);box off;
end

%legend(b,{'0.2m','0.3m','0.4m','0.5m','0.6m','0.7m','0.9m','1.2m'},'location','best');
set(gcf,'position',[0,0,123,107]);
%%
clc;close all;clear all;
ExpName = 'Exp2';
tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));
addpath(genpath('./Functions'));

FilesName = dir(['.\Raw Data\' ExpName filesep '*.txt']);
for i = 1:size(FilesName,1)
    tempFilename{i} = FilesName(i).name(1:end-6);
end
newFilename=unique(tempFilename);
BlockFilesName = dir(['.\Raw Data\' ExpName filesep,char(newFilename(1)),'_*.txt']);
temp = importfile_Exp1(['.\Raw Data\' ExpName filesep,BlockFilesName(1).name]);
temp.VisualAngle = temp.VisualAngle.*(-temp.VisualAngleDirection+0.5)*2;
Condition_Distance = unique(temp.Distance);
Condition_Theta = unique(abs(temp.VisualAngle));
clear temp
for subi = 1:size(newFilename,2)
    BlockFilesName = dir(['.\Raw Data\' ExpName filesep,char(newFilename(subi)),'_*.txt']);
    disp(['processing: sub' num2str(subi)])
    index_count = 1;
    for blocki = 1:size(BlockFilesName,1)
        if blocki == 1
            temp = importfile_Exp1(['.\Raw Data\' ExpName filesep,BlockFilesName(blocki).name]);
            temp.VisualAngle = temp.VisualAngle.*(-temp.VisualAngleDirection+0.5)*2;
        else
            temp2 = importfile_Exp1(['.\Raw Data\' ExpName filesep,BlockFilesName(blocki).name]);
            temp2.VisualAngle = temp2.VisualAngle.*(-temp2.VisualAngleDirection+0.5)*2;
            temp = [temp;temp2];
        end
    end
    temp_all{subi} = temp;
    for  Distance_i = 1:length(Condition_Distance)
        for  Theta1_i = 1:length(Condition_Theta)
            index = find(temp.Distance == Condition_Distance(Distance_i) & abs(temp.VisualAngle) == Condition_Theta(Theta1_i));
            if isempty(index)
                SubAll_a{subi}.Distance = [SubAll_a{subi}.Distance Condition_Distance(Distance_i)];
                SubAll_a{subi}.VisualAngle = [SubAll_a{subi}.VisualAngle Condition_Theta(Theta1_i)];
                SubAll_a{subi}.Answer = [SubAll_a{subi}.Answer nan];
                SubAll_a{subi}.Sub = [SubAll_a{subi}.Sub; subi];
            else
                if (index_count == 1)
                    SubAll_a{subi}.Distance = Condition_Distance(Distance_i);
                    SubAll_a{subi}.VisualAngle = Condition_Theta(Theta1_i);
                    SubAll_a{subi}.Answer = nanmean(temp.Answer(index));
                    SubAll_a{subi}.Sub = subi;
                else
                    SubAll_a{subi}.Distance = [SubAll_a{subi}.Distance Condition_Distance(Distance_i)];
                    SubAll_a{subi}.VisualAngle = [SubAll_a{subi}.VisualAngle Condition_Theta(Theta1_i)];
                    SubAll_a{subi}.Answer = [SubAll_a{subi}.Answer nanmean(temp.Answer(index))];
                    SubAll_a{subi}.Sub = [SubAll_a{subi}.Sub; subi];
                end
            end
            index_count = index_count+1;
        end
    end
end
for subi = 1:size(newFilename,2)
    exp2_Anova.Distance(:,subi) = SubAll_a{subi}.Distance';
    exp2_Anova.VisualAngle(:,subi) = SubAll_a{subi}.VisualAngle';
    exp2_Anova.Answer(:,subi) = SubAll_a{subi}.Answer';
    exp2_Anova.Sub(:,subi) = SubAll_a{subi}.Sub';
end
%%% statistic
clc
o1_cond = abs(unique((exp2_Anova.VisualAngle(:,1))));
d_cond = unique((exp2_Anova.Distance(:,1)));
sub_cond = unique((exp2_Anova.Sub));

for subi = 1:size(sub_cond,1)
    for o1i = 1:size(o1_cond,1)
            for di = 1:size(d_cond,1)
                exp2(subi,o1i,di) = mean(exp2_Anova.Answer(exp2_Anova.VisualAngle(:,subi)==o1_cond(o1i) &...
                    exp2_Anova.Distance(:,subi)==d_cond(di),subi));
            end
    end
end

[tbl,~]=simple_mixed_anova(exp2,[],{'Theta1','Distance'});
partial_Eta_squared1 = ES_Eta_squared(tbl,2);
