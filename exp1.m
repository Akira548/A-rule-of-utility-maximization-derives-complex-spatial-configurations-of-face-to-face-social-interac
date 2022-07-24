%% ImportData a
clc;close all;clear all;
ExpName = 'Exp1';
tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));
addpath(genpath('./Functions'));

%%% theta1 potential fitting
FilesName = dir(['.\Raw Data\' ExpName filesep '*.txt']);
for i = 1:size(FilesName,1)
    tempFilename{i} = FilesName(i).name(1:end-6);
end
newFilename=unique(tempFilename);
BlockFilesName = dir(['.\Raw Data\' ExpName filesep,char(newFilename(1)),'_*.txt']);
temp = importfile_Exp1(['.\Raw Data\' ExpName filesep,BlockFilesName(1).name]);
Condition_Distance = unique(temp.Distance);
Condition_Theta1 = unique(temp.VisualAngle);
Condition_Theta2 = unique(temp.A_Angle);

for subi = 1:size(newFilename,2)
    BlockFilesName = dir(['.\Raw Data\' ExpName filesep,char(newFilename(subi)),'_*.txt']);
    index_count = 1;
    for blocki = 1:size(BlockFilesName,1)
        if blocki == 1
            temp = importfile_Exp1(['.\Raw Data\' ExpName filesep,BlockFilesName(blocki).name]);
        else
            temp = [temp;importfile_Exp1(['.\Raw Data\' ExpName filesep,BlockFilesName(blocki).name])];
        end
    end
    for  Distance_i = 1:length(Condition_Distance)
        for  Theta1_i = 1:length(Condition_Theta1)
            index = find(temp.Distance == Condition_Distance(Distance_i) & temp.VisualAngle == Condition_Theta1(Theta1_i) );
            if (index_count == 1)
                SubAll_a{subi}.Distance = unique(temp.Distance(index));
                SubAll_a{subi}.VisualAngle = unique(temp.VisualAngle(index));
                SubAll_a{subi}.Answer = nanmean(temp.Answer(index));
                SubAll_a{subi}.Sub = subi;
            else
                SubAll_a{subi}.Distance = [SubAll_a{subi}.Distance unique(temp.Distance(index))];
                SubAll_a{subi}.VisualAngle = [SubAll_a{subi}.VisualAngle unique(temp.VisualAngle(index))];
                SubAll_a{subi}.Answer = [SubAll_a{subi}.Answer nanmean(temp.Answer(index))];
                SubAll_a{subi}.Sub = [SubAll_a{subi}.Sub; subi];
            end
            index_count = index_count+1;
        end
    end
end
for subi = 1:size(newFilename,2)
    exp1a.Distance(:,subi) = SubAll_a{subi}.Distance';
    exp1a.VisualAngle(:,subi) = SubAll_a{subi}.VisualAngle';
    exp1a.Answer(:,subi) = SubAll_a{subi}.Answer';
    exp1a.Sub(:,subi) = SubAll_a{subi}.Sub';
end
disp('done!')

%%
max_dis=10;
dis_condition = round(unique(exp1a.Distance),1);
c = gray(length(dis_condition)+2);
temp_va=mean(exp1a.VisualAngle,2);
temp_dis=mean(exp1a.Distance,2);
temp_p=mean(exp1a.Answer,2);
temp_std=nanstd(exp1a.Answer,[],2);
for disi = 1:length(dis_condition)
    b(disi) = plot(temp_va(temp_dis == dis_condition(disi)),temp_p(temp_dis == dis_condition(disi)),'linewidth',1,'Color',c(disi,:));hold on
    errorbar(temp_va(temp_dis == dis_condition(disi)),temp_p(temp_dis == dis_condition(disi)),...
        temp_std(temp_dis == dis_condition(disi))./sqrt(size(exp1a.Answer,2)),...
        '.','MarkerSize',1,'CapSize',10,'linewidth',1,'Color',c(disi,:));hold on;
    ylabel('Probability');xlabel('ж╚(бу)');
    set(gca,'XTick',-60:15:60);xlim([-63 63]);ylim([0 1]);box off;
end
set(gcf,'position',[200,200,600,550]);

figure
[Theta,Rho,p] = griddata(deg2rad(temp_va),temp_dis,temp_p,linspace(-pi/3,pi/3)',linspace(0,max_dis),'linear');
[x,y,z] = pol2cart(Theta,Rho,p); surf(x,y,z);hold on
[Theta2,Rho2,p2] = griddata(deg2rad(temp_va),temp_dis,temp_p,linspace(-pi/3,pi/3)',linspace(0,min(dis_condition)+0.1),'V4');
[x2,y2,z2] = pol2cart(Theta2,Rho2,p2);surf(x2,y2,z2);hold on;
axis equal;view(90,-90);shading interp;caxis([0 1]);
box off;grid off;
axis off;

set(gca,'YTick',-10:2:10);set(gca,'YTicklabel',{'-10','-8','-6','-4','-2','0','2','4','6','8','10(m)'})
set(gca,'XTick',0:2:10);set(gca,'XTicklabel',{'0','2','4','6','8','10(m)'})
xlim([0 10]);

figure
[Theta,Rho,p] = griddata(deg2rad(temp_va),temp_dis,Interaction_Function(deg2rad(temp_va),[],temp_dis,'Patt'),linspace(-pi/3,pi/3)',linspace(0,max_dis),'linear');
[x,y,z] = pol2cart(Theta,Rho,p); surf(x,y,z);hold on
[Theta2,Rho2,p2] = griddata(deg2rad(temp_va),temp_dis,Interaction_Function(deg2rad(temp_va),[],temp_dis,'Patt'),linspace(-pi/3,pi/3)',linspace(0,min(dis_condition)+0.1),'V4');
[x2,y2,z2] = pol2cart(Theta2,Rho2,p2);surf(x2,y2,z2);hold on;
axis equal;view(90,-90);shading interp;caxis([0 1]);
box off;grid off;
axis off;

set(gca,'YTick',-10:2:10);set(gca,'YTicklabel',{'-10','-8','-6','-4','-2','0','2','4','6','8','10(m)'})
set(gca,'XTick',0:2:10);set(gca,'XTicklabel',{'0','2','4','6','8','10(m)'})
xlim([0 10]);
%%
step = 0.01; 
Potentialmaptest = PotentialFieldPlot([-5,5;-5,5],step,pi/90,[0, -1.5, pi/2],'Psoc').*PotentialFieldPlot([-5,5;-5,5],step,pi/90,[0, 1.5,  -pi/3],'Psoc');
% Potentialmaptest = PotentialFieldPlot([-5,5;-5,5],step,pi/90,[0, 1.5,  -pi/6],'Psoc');
% Potentialmaptest = PotentialFieldPlot([-5,5;-5,5],step,pi/90,[0, 0, 0],'Psoc');
% Potentialmaptest = PotentialFieldPlot([-5,5;-5,5],step,pi/90,[0, 0, 0],'Patt');
% Potentialmaptest = PotentialFieldPlot([-5,5;-5,5],step,pi/90,[0, 0, 0],'Prep');
box off;grid off
[s_x,s_y] = meshgrid(-5:step:5,-5:step:5);
surf(s_x,s_y,(Potentialmaptest)');hold on;shading interp;
axis equal
xlim([-5,5]);ylim([-5,5]);

view(0,90);
set(gca,'DataAspectRatio',[1 1 1]);
% lightangle(gca,0,10);material([0.8 1 0.6]);
c = colorbar;
c.Label.String = 'Social Potential ';box off;axis off
%%
step = 0.1; 
Potentialmaptest = PotentialFieldPlot([-5,5;-5,5],step,pi/90,[0, 0, 0],'Psoc');
% Potentialmaptest = CalculatePotentialField([-5,5;-5,5],step,pi/90,[0, 0, 0],0);
% Potentialmaptest = PotentialFieldPlot([-5,5;-5,5],step,pi/90,[0, 0, 0],'Patt');
% Potentialmaptest = 1-PotentialFieldPlot([-5,5;-5,5],step,pi/90,[0, 0, 0],'Prep');
box off;grid off
[s_x,s_y] = meshgrid(-5:step:5,-5:step:5);
% mesh(s_x,s_y,(Potentialmaptest)','edgecolor',[0 0 0],'facecolor','none');hold on;
surf(s_x,s_y,(Potentialmaptest)');hold on;
shading interp;


axis equal
xlim([-5,5]);ylim([-5,5]);

% view(0,90);
view(25,75);
set(gca,'DataAspectRatio',[1 1 .1]);
lightangle(gca,0,00);material([0.8 1 0.1]);
c = colorbar;
c.Label.String = 'Social Potential ';box off;axis off;alpha(.8)
%%
step = 0.1; 
% Potentialmaptest = PotentialFieldPlot([-5,5;-5,5],step,pi/90,[0, 0, 0],'Psoc');
% Potentialmaptest = CalculatePotentialField([-5,5;-5,5],step,pi/90,[0, 0, 0],0);
% Potentialmaptest = PotentialFieldPlot([-5,5;-5,5],step,pi/90,[0, 0, 0],'Patt');
Potentialmaptest = PotentialFieldPlot([-5,5;-5,5],step,pi/90,[0, 0, 0],'Prep');
box off;grid off
[s_x,s_y] = meshgrid(-5:step:5,-5:step:5);
% mesh(s_x,s_y,(Potentialmaptest)','edgecolor',[0 0 0],'facecolor','none');hold on;
surf(s_x,s_y,(Potentialmaptest)');hold on;
shading interp;

axis equal
xlim([-5,5]);ylim([-5,5]);

view(0,90);
% view(25,75);
% set(gca,'DataAspectRatio',[1 1 .1]);
% lightangle(gca,0,00);material([0.8 1 0.1]);
c = colorbar;
c.Label.String = 'Social Potential ';box off;axis off;
%%
ExpName = 'Exp1';
tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));
addpath(genpath('./Functions'));

%%% theta1 potential fitting
FilesName = dir(['.\Raw Data\' ExpName filesep '*.txt']);
for i = 1:size(FilesName,1)
    tempFilename{i} = FilesName(i).name(1:end-6);
end
newFilename=unique(tempFilename);
BlockFilesName = dir(['.\Raw Data\' ExpName filesep,char(newFilename(1)),'_*.txt']);
temp = importfile_Exp1(['.\Raw Data\' ExpName filesep,BlockFilesName(1).name]);
Condition_Distance = unique(temp.Distance);
Condition_Theta1 = unique(abs(temp.VisualAngle));

for subi = 1:size(newFilename,2)
    BlockFilesName = dir(['.\Raw Data\' ExpName filesep,char(newFilename(subi)),'_*.txt']);
    index_count = 1;
    for blocki = 1:size(BlockFilesName,1)
        if blocki == 1
            temp = importfile_Exp1(['.\Raw Data\' ExpName filesep,BlockFilesName(blocki).name]);
        else
            temp = [temp;importfile_Exp1(['.\Raw Data\' ExpName filesep,BlockFilesName(blocki).name])];
        end
    end
    for  Distance_i = 1:length(Condition_Distance)
        for  Theta1_i = 1:length(Condition_Theta1)
            index = find(temp.Distance == Condition_Distance(Distance_i) & abs(temp.VisualAngle) == Condition_Theta1(Theta1_i) );
            if (index_count == 1)
                SubAll_a{subi}.Distance = unique(temp.Distance(index));
                SubAll_a{subi}.VisualAngle = unique(abs(temp.VisualAngle(index)));
                SubAll_a{subi}.Answer = nanmean(temp.Answer(index));
                SubAll_a{subi}.Sub = subi;
            else
                SubAll_a{subi}.Distance = [SubAll_a{subi}.Distance unique(temp.Distance(index))];
                SubAll_a{subi}.VisualAngle = [SubAll_a{subi}.VisualAngle unique(abs(temp.VisualAngle(index)))];
                SubAll_a{subi}.Answer = [SubAll_a{subi}.Answer nanmean(temp.Answer(index))];
                SubAll_a{subi}.Sub = [SubAll_a{subi}.Sub; subi];
            end
            index_count = index_count+1;
        end
    end
end

for subi = 1:size(newFilename,2)
    Sub_Anova.Distance(:,subi) = SubAll_a{subi}.Distance';
    Sub_Anova.VisualAngle(:,subi) = SubAll_a{subi}.VisualAngle';
    Sub_Anova.Answer(:,subi) = SubAll_a{subi}.Answer';
    Sub_Anova.Sub(:,subi) = SubAll_a{subi}.Sub';
end

%%% statistic
clc
o1_cond = abs(unique((Sub_Anova.VisualAngle(:,1))));
d_cond = unique((Sub_Anova.Distance(:,1)));
sub_cond = unique((Sub_Anova.Sub));

for subi = 1:size(sub_cond,1)
    for o1i = 1:size(o1_cond,1)
            for di = 1:size(d_cond,1)
                exp1(subi,o1i,di) = mean(Sub_Anova.Answer(Sub_Anova.VisualAngle(:,subi)==o1_cond(o1i) &...
                    Sub_Anova.Distance(:,subi)==d_cond(di),subi));
            end
    end
end

[tbl,~]=simple_mixed_anova(exp1,[],{'Theta1','Distance'});
partial_Eta_squared1 = ES_Eta_squared(tbl,2);