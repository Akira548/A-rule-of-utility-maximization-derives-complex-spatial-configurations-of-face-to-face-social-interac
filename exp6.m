%% ImportData
clc;close all;addpath(genpath('./'));clear all;
ForderName = [1:18,20:26];
tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));
addpath(genpath('./Function'));

for Seqence = 1:size(ForderName,2)
    File = dir(['.\Raw Data\Exp6\trajectory\' num2str(ForderName(Seqence)) '\*.txt']);
    clear  FileName;
    for File_i = 1:length(File)
        FileName{File_i} = File(File_i).name;
    end
    for subi = 1:6
        if ismember( {[num2str(subi) '.txt']},FileName)
            SeqenceData{Seqence,subi} = importfile_trace(['.\Raw Data\Exp6\trajectory\' num2str(ForderName(Seqence)) filesep num2str(subi) '.txt']);
            SeqenceData{Seqence,subi}.Rotation_Y = wrapToPi(deg2rad(90-SeqenceData{Seqence,subi}.Rotation_Y));
        else
            SeqenceData{Seqence,subi} = [];
        end
    end
end



%% plot

JoinSeqenceName = [2:4,6:8,10:12,14:16];
FormSeqenceName = [23:24,2,5,9,13,21,22,25,26];
JoinTimeLabel = [6,10,8, 6,7,6, 7,8,5, 5,8,7];
JoinWalkerLabel = [4,5,6, 3,2,6, 4,2,1, 6,2,1];

GS=[4 6 4 5 6 4 5 6 4 6];

GroupMember = 1:6;

c = prism(6);
for Seqence = 1:length(JoinSeqenceName)
    seq = (0:300:JoinTimeLabel(Seqence)*300)+1;
    seq_full = (0:1:JoinTimeLabel(Seqence)*300)+1;
    q = [];x = [];y = [];o = [];
    for member = 1:6
        if ~isempty(SeqenceData{JoinSeqenceName(Seqence),member}) && member ~= JoinWalkerLabel(Seqence)
            x = [x,nanmean(SeqenceData{JoinSeqenceName(Seqence),member}.X(JoinTimeLabel(Seqence)*300+1))];
            y = [y,nanmean(SeqenceData{JoinSeqenceName(Seqence),member}.Z(JoinTimeLabel(Seqence)*300+1))];
            v_of_o = mean([cos(SeqenceData{JoinSeqenceName(Seqence),member}.Rotation_Y(JoinTimeLabel(Seqence)*300+1)),sin(SeqenceData{JoinSeqenceName(Seqence),member}.Rotation_Y(JoinTimeLabel(Seqence)*300+1))],1);
            [mean_o,~]=cart2pol(v_of_o(1),v_of_o(2));
            o = [o,mean_o];
            Seqence_Path{Seqence}.x{member} = SeqenceData{JoinSeqenceName(Seqence),member}.X(JoinTimeLabel(Seqence)*300+1);
            Seqence_Path{Seqence}.y{member} = SeqenceData{JoinSeqenceName(Seqence),member}.Z(JoinTimeLabel(Seqence)*300+1);
            Seqence_Path{Seqence}.o{member} = SeqenceData{JoinSeqenceName(Seqence),member}.Rotation_Y(JoinTimeLabel(Seqence)*300+1);
        else
            Seqence_Path{Seqence}.x{member} = [];
            Seqence_Path{Seqence}.y{member} = [];
            Seqence_Path{Seqence}.o{member} = [];
        end
    end
    [Potential{Seqence}, max_Ori1{Seqence}, AR2_record{Seqence}] = CalculatePotentialField([-4,4;-3,3],.1,pi/36,[x',y',o'],0,1);
    Seqence_mean_Path{Seqence} = [x',y',o'];
end
%%%%

for Seqence = 1:length(JoinSeqenceName)
    seq = (0:300:JoinTimeLabel(Seqence)*300)+1;
    q = [];x = [];y = [];o = [];
    for member = 1:6
        if member == JoinWalkerLabel(Seqence)
            Test{Seqence}.x = SeqenceData{JoinSeqenceName(Seqence),member}.X(JoinTimeLabel(Seqence)*300+1);
            Test{Seqence}.y = SeqenceData{JoinSeqenceName(Seqence),member}.Z(JoinTimeLabel(Seqence)*300+1);
            Test{Seqence}.o = SeqenceData{JoinSeqenceName(Seqence),member}.Rotation_Y(JoinTimeLabel(Seqence)*300+1);
            Test{Seqence}.v = [cos(SeqenceData{JoinSeqenceName(Seqence),member}.Rotation_Y(JoinTimeLabel(Seqence)*300+1)),sin(SeqenceData{JoinSeqenceName(Seqence),member}.Rotation_Y(JoinTimeLabel(Seqence)*300+1))];
            [max_x_index,max_y_index] = find(imregionalmax( Potential{Seqence} ));
            max_x = (max_x_index-1)/10-4;
            max_y = (max_y_index-1)/10-3;
            [o_dis(Seqence),max_index] = min(sqrt((max_x - Test{Seqence}.x).^2 + (max_y - Test{Seqence}.y).^2));
            ro(Seqence) = abs(wrapToPi(Test{Seqence}.o - max_Ori1{Seqence}( round((Test{Seqence}.x+4)*10+1),round((Test{Seqence}.y+3)*10+1)) ));
        end
    end
end

figure
c = parula(10);
h(1)=scatter(1:length(o_dis),o_dis,40,c(6,:),'linewidth',2);hold on
h(2)=plot([0,length(o_dis)],[nanmean(o_dis),nanmean(o_dis)],'color',[c(6,:),0.5],'LineWidth',4);
xlabel({'Seqence'});xlim([0 length(o_dis)]);ylim([-0.08 .5]);
ylabel('Centrifugal Distance(m)');grid off;box off;
legend(h,{'Experiment','Experiment Average'},'Location','best');

figure
h2(1)=scatter(1:length(o_dis),rad2deg(ro),40,c(6,:),'linewidth',2);hold on
h2(2)=plot([0,length(o_dis)],[rad2deg(nanmean(ro)),rad2deg(nanmean(ro))],'color',[c(6,:),0.5],'LineWidth',4);
ylim([-3 40]);
xlabel({'Seqence'});xlim([0 length(o_dis)]);
ylabel('Relative Angle(бу)');grid off;box off;
legend(h2,{'Experiment','Experiment Average'},'Location','best');

%%%%%%%%%%%%%  brief report  %%%%%%%%%%%%%%%%%
% disp(['exp_dis:' num2str(nanmean(o_dis,1:2))  ' SD=' num2str(round(nanstd(o_dis),2))  '       exp_ang:'  num2str(rad2deg(mean(ro)))  ' SD=' num2str(round(nanstd(rad2deg(ro)),2))])
MD_CI(o_dis,1);
MD_CI(rad2deg(ro),1);

%%
all_x = [];
xmin=-4;xmax=4;ymin=-3;ymax=3,Dis_step=.1,Ori_step=pi/36;
for Seqence = 1:length(JoinSeqenceName)
    if isempty(all_x)
        all_x = Test{Seqence}.x;
        all_y = Test{Seqence}.y;
        all_o = wrapToPi(Test{Seqence}.o);
    else
        all_x = [all_x; Test{Seqence}.x];
        all_y = [all_y; Test{Seqence}.y];
        all_o = [all_o; wrapToPi(Test{Seqence}.o)];
    end
end

for tryi = 1:1000 % turn of shuffles
    for Seqence = 1:length(JoinSeqenceName)
        [peak_x_index,peak_y_index,~] = find(imregionalmax( Potential{Seqence} ));
        peak_x = (peak_x_index-1).*Dis_step+xmin ;peak_y = (peak_y_index-1).*Dis_step+ymin;
        Data_allsub_x = [];Data_allsub_y=[];
        Localmaximal{Seqence} = [peak_x, peak_y];
        randindex = randi([1 length(all_x)]);
        for i = 1:length(peak_x)
            peak_z{Seqence}(i) = Potential{Seqence}(peak_x_index(i),peak_y_index(i));
            avg_pos_rand = [all_x(randindex), all_y(randindex)];
            o_dis_rand(Seqence,tryi) = sqrt( (avg_pos_rand(1)-peak_x(i)).^2 + (avg_pos_rand(2)-peak_y(i)).^2 )';
        end
        avg_v_rand = [cosd(90-all_o(randindex)), sind(90-all_o(randindex))];
        [rand_o,~]= cart2pol(avg_v_rand(1), avg_v_rand(2));
        rand_pos_x_index = round((all_x(randindex)-xmin)/Dis_step+1);
        rand_pos_y_index = round((all_y(randindex)-ymin)/Dis_step+1);
        avg_raletive_angle_rand(Seqence,tryi) = abs(wrapToPi(rand_o - max_Ori1{Seqence}(rand_pos_x_index,rand_pos_y_index)));
    end
end
figure
h=bar(1:2,[mean(o_dis,'all') mean(o_dis_rand,'all')],0.6,'FaceColor',[.5 .5 .5]);hold on
[~,~,CI1]=ttest(o_dis(:));
[~,~,CI2]=ttest(o_dis_rand(:));
her = errorbar(1:2,[mean(o_dis,'all') mean(o_dis_rand,'all')],[(CI1(2)-CI1(1))/2 (CI2(2)-CI2(1))/2],'k') ;
her.LineStyle = 'none';
% ylim([0 0.8]);
ylabel('Deviation Distance (m)');grid off;box off;
xticklabels({'Expeiment','Scrambled'});
set(gcf,'position',[0,0,600,400]);

figure
h=bar(1:2,[rad2deg(mean(ro,'all')) rad2deg(mean(avg_raletive_angle_rand,'all'))],0.6,'FaceColor',[.5 .5 .5]);hold on
[~,~,CI1]=ttest(rad2deg(ro(:)));
[~,~,CI2]=ttest(rad2deg(avg_raletive_angle_rand(:)));
her = errorbar(1:2,[rad2deg(mean(ro,'all')) rad2deg(mean(avg_raletive_angle_rand,'all'))],[(CI1(2)-CI1(1))/2 (CI2(2)-CI2(1))/2],'k') ;
her.LineStyle = 'none';
% ylim([0 23]);
ylabel('Relative Angle(бу)');grid off;box off;
xticklabels({'Expeiment','Scrambled'});
set(gcf,'position',[0,0,600,400]);

%%%%%%%%%%%%%  brief report  %%%%%%%%%%%%%%%%%
MD_CI(o_dis_rand,1);
MD_CI(rad2deg(avg_raletive_angle_rand),1);


%%
load('.\Data\exp5b.mat');
figure
h=bar([mean(deviation_dis_exp5b,'all') mean(o_dis)],0.6,'FaceColor',[.5 .5 .5]);hold on
her = errorbar(1:2,[mean(deviation_dis_exp5b,'all') mean(o_dis)],[nanstd(mean(deviation_dis_exp5b,1))./sqrt(size(mean(deviation_dis_exp5b,1),2)) nanstd(o_dis)./sqrt(size(o_dis,2))],'k') ;
her.LineStyle = 'none';
ylim([0 0.8]);
ylabel('Deviation Distance (m)');grid off;box off;
xticklabels({'VR','Real'});
set(gcf,'position',[0,0,600,400]);

figure
h=bar([rad2deg(mean(RA_mean_exp5b,'all')) rad2deg(mean(ro))],0.6,'FaceColor',[.5 .5 .5]);hold on
her = errorbar(1:2,[rad2deg(mean(RA_mean_exp5b,'all')) rad2deg(mean(ro))],[nanstd(mean(RA_mean_exp5b,1))./sqrt(size(mean(RA_mean_exp5b,1),2)) nanstd(rad2deg(ro))./sqrt(size(ro,2))],'k') ;
her.LineStyle = 'none';
ylim([0 23]);
ylabel('Relative Angle(бу)');grid off;box off;
xticklabels({'VR','Real'});
set(gcf,'position',[0,0,600,400]);

%% Generate Path
close all

JoinSeqenceName = [2:4,6:8,10:12,14:16];
FormSeqenceName = [23:24,2,5,9,13,21,22,25,26];
JoinTimeLabel = [6,10,8, 6,7,6, 7,8,5, 5,8,7];
JoinWalkerLabel = [4,5,6, 3,2,6, 4,2,1, 6,2,1];

GroupMember = 1:6;
for Seqence = 1:length(JoinSeqenceName)
% for Seqence = 6:6
    seq = (0:100:JoinTimeLabel(Seqence)*300)+1;
    q = [];x = [];y = [];o = [];
    for member = 1:6
        if member == JoinWalkerLabel(Seqence)
            Seqence_Path{Seqence}.x{member} = SeqenceData{JoinSeqenceName(Seqence),member}.X(seq);
            Seqence_Path{Seqence}.y{member} = SeqenceData{JoinSeqenceName(Seqence),member}.Z(seq);
            Seqence_Path{Seqence}.o{member} = SeqenceData{JoinSeqenceName(Seqence),member}.Rotation_Y(seq);
            StartPoint = [SeqenceData{JoinSeqenceName(Seqence),member}.X(1),SeqenceData{JoinSeqenceName(Seqence),member}.Z(1)];
            ModelPath{Seqence}=Generate_Path(Potential{Seqence},[-4,4;-3,3],.1,StartPoint,.4);
            f = fit(Seqence_Path{Seqence}.y{member},Seqence_Path{Seqence}.x{member},'linearinterp');
            ft=fitlm(ModelPath{Seqence}(:,1),f(ModelPath{Seqence}(:,2)));
%             Rsquared_of_Path(Seqence) = double(ft.Rsquared.Ordinary);
           [Rsquared_of_Path(Seqence),p_of_Path(Seqence),F_of_Path(Seqence),frechet_dis_of_Path(Seqence) ]= Path_compare([Seqence_Path{Seqence}.x{member},Seqence_Path{Seqence}.y{member}],ModelPath{Seqence},StartPoint);
            %%% path plot
            subplot(2,6,Seqence)
            surf(-4:.1:4,-3:.1:3,Potential{Seqence}'-1);hold on;shading interp
            plot(Seqence_Path{Seqence}.x{member},Seqence_Path{Seqence}.y{member},'.r');hold on
            plot(ModelPath{Seqence}(:,1),(ModelPath{Seqence}(:,2)),'.c');
%             plot(f(Seqence_Path{Seqence}.y{member}),Seqence_Path{Seqence}.y{member},'-','color',[1,0,0]);hold on
            scatter(Seqence_mean_Path{Seqence}(:,1),Seqence_mean_Path{Seqence}(:,2),50,'.','MarkerEdgeColor','k');hold on
            temp_q=quiver(Seqence_mean_Path{Seqence}(:,1),Seqence_mean_Path{Seqence}(:,2),cos(Seqence_mean_Path{Seqence}(:,3)),sin(Seqence_mean_Path{Seqence}(:,3)),0.5,'Color',[.5,.5,.5],'linewidth',.5);hold on
%             q = [q,temp_q];
%             scale_quivers(q,.02);
            
            title(['Seqence: ' num2str(Seqence)]);
%             xlabel(['R^2 = ' num2str(round(double(ft.Rsquared.Ordinary),3))]);
            xlim([-4 4]);ylim([-3 3]);axis equal;view(0,90)
            set(gcf,'position',[0,200,1400,400]);
            legend off
            
        end
    end
end

%%%%%%%%%%%%%  brief report  %%%%%%%%%%%%%%%%%
disp(['Rsquared:' num2str(round(mean(Rsquared_of_Path),2) )  ' SD=' num2str(round(std(Rsquared_of_Path),2))  '       frechet_dis:'  num2str(round(mean(frechet_dis_of_Path),2))  ' SD=' num2str(round(std(frechet_dis_of_Path),2))])
MD_CI(Rsquared_of_Path,1);
MD_CI(frechet_dis_of_Path,1);
%%
xmin = -5;xmax = 3;ymin = -3;ymax = 5;Dis_step = 0.1;Ori_step = pi/36;
Seqence = 6
seq = (0:300:JoinTimeLabel(Seqence)*300)+1;
seq_full = (0:1:JoinTimeLabel(Seqence)*300)+1;
q = [];x = [];y = [];o = [];
for member = 1:6
    if ~isempty(SeqenceData{JoinSeqenceName(Seqence),member}) && member ~= JoinWalkerLabel(Seqence)
        x = [x,nanmean(SeqenceData{JoinSeqenceName(Seqence),member}.X(JoinTimeLabel(Seqence)*300+1))];
        y = [y,nanmean(SeqenceData{JoinSeqenceName(Seqence),member}.Z(JoinTimeLabel(Seqence)*300+1))];
        v_of_o = mean([cos(SeqenceData{JoinSeqenceName(Seqence),member}.Rotation_Y(JoinTimeLabel(Seqence)*300+1)),sin(SeqenceData{JoinSeqenceName(Seqence),member}.Rotation_Y(JoinTimeLabel(Seqence)*300+1))],1);
        [mean_o,~]=cart2pol(v_of_o(1),v_of_o(2));
        o = [o,mean_o];
        Seqence_Path{Seqence}.x{member} = SeqenceData{JoinSeqenceName(Seqence),member}.X(JoinTimeLabel(Seqence)*300+1);
        Seqence_Path{Seqence}.y{member} = SeqenceData{JoinSeqenceName(Seqence),member}.Z(JoinTimeLabel(Seqence)*300+1);
        Seqence_Path{Seqence}.o{member} = SeqenceData{JoinSeqenceName(Seqence),member}.Rotation_Y(JoinTimeLabel(Seqence)*300+1);
    else
        Seqence_Path{Seqence}.x{member} = [];
        Seqence_Path{Seqence}.y{member} = [];
        Seqence_Path{Seqence}.o{member} = [];
    end
    if member == JoinWalkerLabel(Seqence)
        Test{Seqence}.x = SeqenceData{JoinSeqenceName(Seqence),member}.X(JoinTimeLabel(Seqence)*300+1);
        Test{Seqence}.y = SeqenceData{JoinSeqenceName(Seqence),member}.Z(JoinTimeLabel(Seqence)*300+1);
        Test{Seqence}.o = SeqenceData{JoinSeqenceName(Seqence),member}.Rotation_Y(JoinTimeLabel(Seqence)*300+1);
        Test{Seqence}.v = [cos(SeqenceData{JoinSeqenceName(Seqence),member}.Rotation_Y(JoinTimeLabel(Seqence)*300+1)),sin(SeqenceData{JoinSeqenceName(Seqence),member}.Rotation_Y(JoinTimeLabel(Seqence)*300+1))];
        [max_x_index,max_y_index] = find(imregionalmax( Potential{Seqence} ));
        max_x = (max_x_index-1)/10-4;
        max_y = (max_y_index-1)/10-3;
        [o_dis(Seqence),max_index] = min(sqrt((max_x - Test{Seqence}.x).^2 + (max_y - Test{Seqence}.y).^2));
        ro(Seqence) = abs(wrapToPi(Test{Seqence}.o - max_Ori1{Seqence}( round((Test{Seqence}.x+4)*10+1),round((Test{Seqence}.y+3)*10+1)) ));
    end

end

figure
[U, Ori, Ori_record] = CalculatePotentialField([xmin,xmax;ymin,ymax],Dis_step,Ori_step,[x',y',o'],0,1);
surf(xmin:Dis_step:xmax,ymin:Dis_step:ymax,U');hold on;view(0,90);shading interp;
axis equal;ylim([ymin ymax]);xlim([xmin xmax]);box off;grid off;
colormap(colorbar_BuOr());

scatter(Test{Seqence}.x,Test{Seqence}.y,50,'xk');hold on
axis  equal;ylim([ymin ymax]);xlim([xmin xmax]);

Test{Seqence}.x = SeqenceData{JoinSeqenceName(Seqence),member}.X(JoinTimeLabel(Seqence)*300+1);
Test{Seqence}.y = SeqenceData{JoinSeqenceName(Seqence),member}.Z(JoinTimeLabel(Seqence)*300+1);
%%
for tryi = 1:1000 % turn of shuffles
    randorder = randperm(length(all_x));
    for i=1:size(avg_pos,1)
        avg_pos_rand = [all_x(randorder(i)), all_y(randorder(i))];
        o_dis_rand(i,tryi) = sqrt( (avg_pos_rand(1) - max_pos(i,1)).^2 + (avg_pos_rand(2) - max_pos(i,2)).^2);
    end
end

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


%%

pos_index = [round((Test{Seqence}.x-xmin)./Dis_step+1)  round((Test{Seqence}.y-ymin)./Dis_step+1)];
ix_seq = xmin:Dis_step:xmax;
iy_seq = ymin:Dis_step:ymax;
temp_Rho = squeeze(Ori_record(round((Test{Seqence}.x-xmin)/Dis_step+1),round(((Test{Seqence}.y)-ymin)/Dis_step+1),:));

for ix = 1:length(ix_seq)
    for iy = 1:length(iy_seq)
        [Theta,Rho]= cart2pol((ix_seq(ix)-Test{Seqence}.x),(iy_seq(iy)-Test{Seqence}.y));
        ori_strength(ix,iy)=temp_Rho( round((Theta)/(Ori_step) + pi/Ori_step + 1) );
    end
end

scatter(Test{Seqence}.x,Test{Seqence}.y,50,'.k');hold on
plot([Test{Seqence}.x+5*cos(nanmean(wrapToPi(Test{Seqence}.o))), Test{Seqence}.x],[Test{Seqence}.y+5*sin(nanmean(wrapToPi(Test{Seqence}.o))),Test{Seqence}.y],'--k','linewidth',2);

% surf(ix_seq,iy_seq,ori_strength(:,:)');hold on;shading interp
grid off;axis equal;view(0,90);ylim([ymin ymax]);xlim([xmin xmax]);
caxis([0 max(ori_strength,[],'all')]);
colormap jet;