%% ImportData [Experiment 4]
clc;close all;clear all;
ExpName = 'Exp3b';
tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));
addpath(genpath('./Functions'));
FilesName = dir(['.\Raw Data\' ExpName filesep '*.txt']);

for subi = 1:size(FilesName,1)
    tempdata = importfile_Exp3b(['.\Raw Data\' ExpName filesep FilesName(subi).name]);
    oricondition = unique(tempdata.ORientaiton);
    abs_oricondition = unique(abs(tempdata.ORientaiton));
    Sub_pos_Data{subi}.ori = tempdata.ORientaiton;
    Sub_pos_Data{subi}.x = tempdata.x;
    Sub_pos_Data{subi}.z = tempdata.z;
    Sub_pos_Data{subi}.heading = deg2rad(90-tempdata.Ortation_y);
    Sub_pos_Data{subi}.odis = sqrt((tempdata.x).^2+(tempdata.z - 3).^2);
    for orii = 1:length(oricondition)
        index = find(tempdata.ORientaiton == oricondition(orii))';
        for index_i = 1:length(index)
            SubData{subi, index_i, orii} = importfile_path(['.\Raw Data\' ExpName '\Log\' FilesName(subi).name(isstrprop(FilesName(subi).name, 'digit')) filesep num2str(index(index_i)) '.txt']);
        end
    end
end

%%% Interaction Position Analysis [Experiment 4]
for orii = 1:length(oricondition)
    for subi = 1:size(FilesName,1)
        index = find(Sub_pos_Data{subi}.ori == oricondition(orii));
        mean_x(subi,orii) = mean(Sub_pos_Data{subi}.x(index));
        mean_y(subi,orii) = mean(Sub_pos_Data{subi}.z(index));
        AllData(subi,orii,:) =mean( [Sub_pos_Data{subi}.x(index).*cosd(-Sub_pos_Data{subi}.ori(index)) + (Sub_pos_Data{subi}.z(index)-3).*sind(-Sub_pos_Data{subi}.ori(index))...
            ,(Sub_pos_Data{subi}.z(index)-3).*cosd(-Sub_pos_Data{subi}.ori(index))-Sub_pos_Data{subi}.x(index).*sind(-Sub_pos_Data{subi}.ori(index))],1);
        AllData2(subi,orii,:,:) =( [Sub_pos_Data{subi}.x(index).*cosd(-Sub_pos_Data{subi}.ori(index)) + (Sub_pos_Data{subi}.z(index)-3).*sind(-Sub_pos_Data{subi}.ori(index))...
            ,(Sub_pos_Data{subi}.z(index)-3).*cosd(-Sub_pos_Data{subi}.ori(index))-Sub_pos_Data{subi}.x(index).*sind(-Sub_pos_Data{subi}.ori(index))]);
        if (subi == 1)
            all_respond{orii} = Sub_pos_Data{subi}.x(index)';
        else
            all_respond{orii} = [all_respond{orii} ,Sub_pos_Data{subi}.x(index)'];
        end
    end
end
 temp_x=[AllData2(:,:,1,1),AllData2(:,:,2,1)];
 temp_y=[AllData2(:,:,1,2),AllData2(:,:,2,2)];
%%
VH_pos = [0,3];
for orii = 1:length(oricondition)
    [Potentialmap{orii}, max_Ori1{orii}, AR2_record{orii}] = CalculatePotentialField([-3,3;0,6],0.1,pi/36,[0, 3,wrapToPi(-pi/2-deg2rad(oricondition(orii)))],0,1);
      for subi = 1:size(FilesName,1)
        index = find((Sub_pos_Data{subi}.ori) == oricondition(orii));
        pos(subi,orii,:) = [nanmean((Sub_pos_Data{subi}.x(index))) nanmean(Sub_pos_Data{subi}.z(index))];
        v(subi,orii,:) = mean([cos(Sub_pos_Data{subi}.heading(index)) sin(Sub_pos_Data{subi}.heading(index))]);
      end
    exp3b.d(:,orii) = deviation_of_position(pos(:,orii,1),pos(:,orii,2),Potentialmap{orii},[-3,3;0,6],0.1);
    exp3b.o(:,orii) = rad2deg(deviation_of_orientation(pos(:,orii,1),pos(:,orii,2),v(:,orii,1),v(:,orii,2),max_Ori1{orii},[-3,3;0,6],0.1));
end
%%
absoricondition = unique(abs(oricondition));
for i = 1:length(absoricondition)
   exp3b_anaova.d(:,i) = mean(exp3b.d(:,find(abs(oricondition) == absoricondition(i))),2);
   exp3b_anaova.o(:,i) = mean(exp3b.o(:,find(abs(oricondition) == absoricondition(i))),2);
end

MD_CI(exp3b.d);
[tbl2,stats2]=simple_mixed_anova(exp3b_anaova.d);
partial_Eta_squared = ES_Eta_squared(tbl2,2);

% deviation distance analysis
MD_CI(exp3b.o);
[tbl3,stats3]=simple_mixed_anova(exp3b_anaova.o);
partial_Eta_squared = ES_Eta_squared(tbl3,2);

%%
VH_pos = [0,3];
for orii = 1:length(oricondition)
    [Potentialmap{orii}, max_Ori1{orii}, AR2_record{orii}] = CalculatePotentialField([-3,3;0,6],0.1,pi/36,[0, 3,wrapToPi(-pi/2-deg2rad(oricondition(orii)))],0,1);
    [max_x,max_y] = find(Potentialmap{orii} == max(Potentialmap{orii},[],'all'));
    for subi = 1:size(FilesName,1)
        index = find((Sub_pos_Data{subi}.ori) == oricondition(orii));
        if isempty(index)
            deviation_dis(subi,orii) =  nan;
            deviation_angle(subi,orii) =  nan;
        else
        deviation_dis(subi,orii) = min(abs(sqrt(( nanmean((Sub_pos_Data{subi}.x(index))) - ((max_x - 1)/10-3)).^2+(nanmean(Sub_pos_Data{subi}.z(index)) - (max_y - 1)/10).^2)));
        v = mean([cos(Sub_pos_Data{subi}.heading(index)) sin(Sub_pos_Data{subi}.heading(index))]);
        MOri=max_Ori1{orii}( round(nanmean((Sub_pos_Data{subi}.x(index)))*10+31) ,round(nanmean(Sub_pos_Data{subi}.z(index))*10+1));
        deviation_angle(subi,orii) = acos(dot(v,[cos(MOri),sin(MOri)]));
        end
    end    
end

for orii2 = 1:length(abs_oricondition)
   index2 =  find(abs(oricondition) == abs_oricondition(orii2));
   deviation_dis_ano(:,orii2) = nanmean(deviation_dis(:,index2),2);
   deviation_angle_ano(:,orii2) = nanmean(deviation_angle(:,index2),2);
end
exp3b.d = (deviation_dis_ano);
exp3b.o = rad2deg(deviation_angle_ano);
% deviation distance analysis
MD_CI(exp3b.d,1);
[tbl2,stats2]=simple_mixed_anova(exp3b.d);
partial_Eta_squared = ES_Eta_squared(tbl2,2);

% deviation angle analysis
MD_CI((exp3b.o(:)),1);
[tbl3,stats3]=simple_mixed_anova(exp3b.o);
partial_Eta_squared = ES_Eta_squared(tbl3,2);

%% 
Potentialmaptest = CalculatePotentialField([-1,7;-4,4],0.05,pi/90,[0, 0, 0],0);
temp_x = mean(AllData2(:,:,:,1),[2 3]);
temp_y = mean(AllData2(:,:,:,2),[2 3]);

scatter3(-temp_y(:),temp_x(:),ones(size(temp_y(:)))*1.1,50,'xk');hold on

box off;grid off
[s_x,s_y] = meshgrid(-1:0.05:7,-4:0.05:4);
% surf(s_x,s_y,(Potentialmaptest)');hold on;shading interp;
axis equal
xlim([-1,7]);ylim([-4,4]);

view(0,90);
set(gca,'DataAspectRatio',[1 1 .7]);
lightangle(gca,0,20);material([0.8 1 0.6]);
c = colorbar;
c.Label.String = 'Social Potential ';box off;axis on

%%
Potentialmaptest = CalculatePotentialField([-6,6;-6,6],0.1,pi/36,[0, 0, 0],0);
[max_x,max_y] = find(Potentialmaptest' == max(Potentialmaptest',[],'all'));
c2 = parula(10);
bar(1:length(abs_oricondition),nanmean(deviation_dis_ano,1),0.8,'FaceColor',[.5 .5 .5],'linewidth',1);hold on
her = errorbar(1:length(abs_oricondition),nanmean(deviation_dis_ano,1),nanstd(deviation_dis_ano,1)'/sqrt(size(deviation_dis_ano,1)),'linewidth',1) ;
her.LineStyle = 'none';her.Color = 'k';
xticklabels({ '0' ,   '30',    '60',    '90' ,  '120' ,  '150' ,  '180(бу)'});
xlabel('Start Direction');
ylim([0,1]);xlim([0,length(abs_oricondition)+1]);box off
ylabel('Deviation Distance');
dis_Range = 0:0.01:3;
plot([0 length(abs_oricondition)+1],[mean(deviation_dis_ano,1:2),mean(deviation_dis_ano,1:2)],'--r');
set(gcf,'position',[0,200,750,400]);

%%  
[Potentialmaptest, max_Ori1, ori_distribution] = CalculatePotentialField([-1,7;-4,4],0.05,pi/180,[0, 0, 0],0,1);
temp_x = mean(AllData2(:,:,:,1),[1 2 3]);
temp_y = mean(AllData2(:,:,:,2),[1 2 3]);

% scatter3(-temp_y,temp_x,ones(size(temp_y))*1.1,50,'xk');hold on
ix_seq = -1:0.05:7;
iy_seq = -4:0.05:4;
Ori_step=pi/180;

for subi = 1:length(Sub_pos_Data)
    if subi == 1
        average_v = mean([cos(Sub_pos_Data{subi}.heading) sin(Sub_pos_Data{subi}.heading)]);
    else
        average_v =[average_v; mean([cos(Sub_pos_Data{subi}.heading) sin(Sub_pos_Data{subi}.heading)])];
    end
end
average_v_all = mean(average_v);
[average_ori,~]=  cart2pol(average_v_all(1),average_v_all(2));

o=[round(-temp_y(:)*20+21);round(temp_x(:)*20+81)];
temp_Rho=squeeze(ori_distribution(o(1),o(2),:));
for ix = 1:length(ix_seq)
    for iy = 1:length(iy_seq)
        [Theta,Rho]=  cart2pol((ix_seq(ix)--temp_y),(iy_seq(iy)-temp_x));
        ori_strength(ix,iy)=temp_Rho( round((Theta)/Ori_step + pi/Ori_step + 1) );
    end
end

surf(ix_seq,iy_seq,ori_strength(:,:)');hold on;shading interp;
% plot([-temp_y+2*cos(pi/2+average_ori),-temp_y],[temp_x+2*sin(pi/2+average_ori),temp_x],'--k','linewidth',2);

view(90,90);
axis equal
xlim([-1,7]);ylim([-4,4]);caxis([0 max(Potentialmaptest,[],'all')]);
c = colorbar;
box off;grid off
% axis off;
colormap(colorbar_BuOr());

%%
close all
Potential_p1 = CalculatePotentialField([-4,8;-2.5,2.5],0.05,pi/180,[0, 0, 0],0);

[max_index_y,max_index_x] =find(Potential_p1== max(Potential_p1,[],'all'));
[dx,dy]=gradient(Potential_p1');
[s_x,s_y] = meshgrid(-4:0.05:8,-2.5:0.05:2.5);
contour3(s_x,s_y,real(Potential_p1)',[0,0.1,0.2,0.3:0.05:0.6 0.625:0.025:0.8],'-k','LineWidth',1);hold on
surf(s_x,s_y,(Potential_p1)');hold on;shading interp;hold on
view(50,45);
axis equal;
xlim([-4,8]);ylim([-2.5,2.5]);
set(gca,'ZDir','reverse','DataAspectRatio',[1 1 .3]);
material([0.65 1 0.6]);
axis off
box off;grid off
lightangle(gca,47,-60)
lighting gouraud
set(gcf,'position',[0,0,1000,800]);

%% Walking Path Analysis: all path & average path [Experiment 5]
Ori_condition = unique(Sub_pos_Data{1}.ori);
[~,order_index] = sort(abs(Ori_condition));
Ori_condition_sort = Ori_condition(order_index);
step = 0.4
close all;clear trace temp_x range subject_mean_x subject_mean_y trace
end_point = [mean(mean_x,1)' ,mean(mean_y,1)' ];
for i_ori = 1:length(oricondition)
    clear range; fitting_trace_x{i_ori} = [];fitting_trace_y{i_ori} = [];
    range = 0:step:round(end_point(i_ori,2),1);
    for subi = 1:size(FilesName,1)
        clear mean_trace_Y mean_trace_X fitting_trace;
        mean_trace_Y = [SubData{subi,1,i_ori}.Z; SubData{subi,2,i_ori}.Z];
        mean_trace_X = [SubData{subi,1,i_ori}.X; SubData{subi,2,i_ori}.X];
        fitting_trace_x{i_ori} = [fitting_trace_x{i_ori};mean_trace_X];
        fitting_trace_y{i_ori} = [fitting_trace_y{i_ori};mean_trace_Y];
        for rangei = 1:length(range)-1
            temp_x{i_ori,rangei} = [];
            subject_mean_x{i_ori}(subi,rangei) = nanmean(mean_trace_X(find((mean_trace_Y >= range(rangei))  & (mean_trace_Y < range(rangei+1)))));
            subject_mean_y{i_ori}(subi,rangei) = (range(rangei)+0.1);
            temp_x{i_ori,rangei} =[temp_x{i_ori,rangei} ;mean_trace_X(find((mean_trace_Y >= range(rangei))  & (mean_trace_Y < range(rangei+1))) )];
            trace{i_ori,rangei} = [nanmean(temp_x{i_ori,rangei}), range(rangei)+0.1];
        end
    end
end

[~,order_index] = sort(abs(oricondition));
for i_ori = 1:length(oricondition)
    %     figure
    subplot(3,4,i_ori);
    scatter(fitting_trace_x{order_index(i_ori)}, fitting_trace_y{order_index(i_ori)},0.1,'.k');hold on;axis equal;title(oricondition(order_index(i_ori)));
    scatter(0,3,50,'.r');hold on;
    plot([0,cosd(-Ori_condition(order_index(i_ori))-90)],[3,3+sind(-Ori_condition(order_index(i_ori))-90)],'-r');hold on
    title(['original path: ' num2str(Ori_condition(order_index(i_ori))) 'бу']);
    axis equal;xlim([-3,3]);ylim([-0.5,5]);
end

%%% average path
figure
clear trace
for i_ori = 1:length(oricondition)
    clear range; fitting_trace_x{i_ori} = [];fitting_trace_y{i_ori} = [];
    range = 0:step:round(end_point(i_ori,2),1);
    for subi = 1:size(FilesName,1)% reshape 2 trials in one dimension
        clear mean_trace_Y mean_trace_X fitting_trace;
        mean_trace_Y1 = SubData{subi,1,i_ori}.Z;
        mean_trace_X1 = SubData{subi,1,i_ori}.X;
        mean_trace_Y2 = SubData{subi,2,i_ori}.Z;
        mean_trace_X2 = SubData{subi,2,i_ori}.X;
        for rangei = 1:length(range)-1
            temp_x{i_ori,rangei} = [];
            subject_mean_x{i_ori}(2*subi-1,rangei) = nanmean(mean_trace_X1(find((mean_trace_Y1 >= range(rangei))  & (mean_trace_Y1 < range(rangei+1)))));
            subject_mean_x{i_ori}(2*subi,rangei) = nanmean(mean_trace_X2(find((mean_trace_Y2 >= range(rangei))  & (mean_trace_Y2 < range(rangei+1)))));
            subject_mean_y{i_ori}(2*subi-1,rangei) = (range(rangei)+0.1);
            subject_mean_y{i_ori}(2*subi,rangei) = (range(rangei)+0.1);
        end
    end
end

for i_ori = 1:length(Ori_condition_sort)
    clear temp_trace
    temp_trace_x  = subject_mean_x{order_index(i_ori)};
    if i_ori == 12
        [~,index] = max(abs(temp_trace_x(:,:)),[],2);
        for subi = 1:size(FilesName,1)*2
            max_value(subi) = temp_trace_x(subi,index(subi));
        end
        index_P = find(max_value>=0);    index_N = find(max_value<0);
    end
    
    for lengthi = 1:size(temp_trace_x,2)  %%% exclude abnormal path
        if i_ori == 12
            temp_mean_left =  nanmean(temp_trace_x(index_N,lengthi));
            temp_std_left = nanstd(temp_trace_x(index_N,lengthi));
            index = find(temp_trace_x(index_N,lengthi)>=(temp_mean_left + 3*temp_std_left) | temp_trace_x(index_N,lengthi) <= (temp_mean_left - 3*temp_std_left));
            temp_trace_x(index_N(index),lengthi) = nan;
            
            temp_mean_right =  nanmean(temp_trace_x(index_P,lengthi));
            temp_std_right = nanstd(temp_trace_x(index_P,lengthi));
            index = find(temp_trace_x(index_P,lengthi)>=(temp_mean_right + 3*temp_std_right) | temp_trace_x(index_P,lengthi) <= (temp_mean_right - 3*temp_std_right));
            temp_trace_x(index_P(index),lengthi) = nan;
        else
            temp_mean =  nanmean(temp_trace_x(:,lengthi));
            temp_std = nanstd(temp_trace_x(:,lengthi));
            index = find(temp_trace_x(:,lengthi)>=(temp_mean + 3*temp_std) | temp_trace_x(:,lengthi) <= (temp_mean - 3*temp_std));
            temp_trace_x(index,lengthi) = nan;
        end
    end
    
    temp_trace_y  = subject_mean_y{order_index(i_ori)};
    subplot(3,4,i_ori);
    scatter(0,3,50,'.r');hold on
    plot([0,cosd(-Ori_condition_sort(i_ori)-90)],[3,3+sind(-Ori_condition_sort(i_ori)-90)],'-r');hold on
    if i_ori == 12 % all subjects mean
        scatter(nanmean(temp_trace_x(index_N,:),1),nanmean(temp_trace_y(index_N,:),1),1,'.k');hold on
        scatter(nanmean(temp_trace_x(index_P,:),1),nanmean(temp_trace_y(index_P,:),1),1,'.k');hold on
        mean_Trace{12} = [nanmean(temp_trace_x(index_N,:),1)',nanmean(temp_trace_y(index_N,:),1)'];
        mean_Trace{13} = [nanmean(temp_trace_x(index_P,:),1)',nanmean(temp_trace_y(index_P,:),1)'];
    else
        scatter(nanmean(temp_trace_x,1),nanmean(temp_trace_y,1),1,'.k');hold on
        mean_Trace{i_ori} = [nanmean(temp_trace_x,1)',nanmean(temp_trace_y,1)'];
    end
    title(['average path: ' num2str(Ori_condition_sort(i_ori)) 'бу']);
    axis equal;xlim([-3,3]);ylim([-0.5,5]);
end

%% path
close all;
Ori_condition = unique(Sub_pos_Data{1}.ori);
[~,order_index] = sort(abs(Ori_condition));
Ori_condition_sort = Ori_condition(order_index);
Dis_step = 0.01;

for i_ori = 1:length(Ori_condition_sort)+1
    %%% 
    switch i_ori
        case 12
            startpoint = [-0.1,0];
            VH_condi = [0, 3, wrapToPi(deg2rad(-Ori_condition_sort(i_ori)-90))];
        case 13
            startpoint = [0.1,0];
            VH_condi = [0, 3, wrapToPi(deg2rad(-Ori_condition_sort(12)-90))];
        otherwise
            startpoint = [0,0];
            VH_condi = [0, 3, wrapToPi(deg2rad(-Ori_condition_sort(i_ori)-90))];
    end
    Potential{i_ori} = CalculatePotentialField([-3,3;-0.5,5],Dis_step,pi/90,[VH_condi(1),VH_condi(2),VH_condi(3)],0);
%     ModelPath{i_ori} = Generate_Path(Potential{i_ori},[-3,3;-0.5,5],.1,startpoint,.4,'Combination');
    ModelPath{i_ori} = Generate_Path(Potential{i_ori},[-3,3;-0.5,5],Dis_step,startpoint,step);
    [Rsquared_of_Path(i_ori),p_of_Path(i_ori),F_of_Path(i_ori),frechet_dis(i_ori)] = Path_compare(mean_Trace{i_ori},ModelPath{i_ori},startpoint);
    
    %%% path plot
    if i_ori ~= 13
        subplot(2,6,i_ori)
        contour(-3:Dis_step:3,-0.5:Dis_step:5,Potential{i_ori}'-1,'k');hold on;
        surf(-3:Dis_step:3,-0.5:Dis_step:5,Potential{i_ori}'-1);hold on;shading interp
        plot(mean_Trace{i_ori}(:,1),mean_Trace{i_ori}(:,2),'.r','Markersize',10);hold on
        plot(ModelPath{i_ori}(:,1),ModelPath{i_ori}(:,2),'.c','Markersize',10);hold on
        title(['Orientation: ' num2str(Ori_condition_sort(i_ori))  'бу'] );
    else
        subplot(2,6,12)
        plot(mean_Trace{i_ori}(:,1),mean_Trace{i_ori}(:,2),'.r','Markersize',10);hold on
        plot(ModelPath{i_ori}(:,1),ModelPath{i_ori}(:,2),'.c','Markersize',10);
    end
    view(0,90);
    legend off;box off;axis equal;xlim([-3,3]);ylim([-0.5,5]);
    set(gcf,'position',[0,200,1200,400]);
end

figure
b=bar(1:length(Ori_condition_sort)+1,Rsquared_of_Path,.75,'FaceColor',[.5,.5,.5]);
xlim([0 length(Ori_condition_sort)+2]);ylim([0 1.09]);box off;
xlabel('Condition');
ylabel('R^2 between Gradient and Human Path');
set(gcf,'position',[0,200,800,300]);
xtips = b.XData;
ytips = b.YData;
labels = string(round(b.YData,3));
text(xtips,ytips,labels,'HorizontalAlignment','center','VerticalAlignment','bottom');

figure
b=bar(1:length(Ori_condition_sort)+1,frechet_dis,.75,'FaceColor',[.5,.5,.5]);hold on
plot([0 length(Ori_condition_sort)+2],[mean(frechet_dis) mean(frechet_dis)],'color',[.5 .5 .5],'linewidth',2);hold on
xlim([0 length(Ori_condition_sort)+2]);box off;
xlabel('Condition');
ylabel('Frechet Distance');
set(gcf,'position',[0,200,800,300]);
xtips = b.XData;
ytips = b.YData;
labels = string(round(b.YData,3));
text(xtips,ytips,labels,'HorizontalAlignment','center','VerticalAlignment','bottom');

MD_CI(Rsquared_of_Path,1);
MD_CI(frechet_dis,1);
%% exp3a vs exp3b
figure
bar(1:2,[nanmean(exp3a.d,1:2),nanmean(exp3b.d,1:2)],0.6,'FaceColor',[.5 .5 .5]);hold on
he = errorbar(1:2,[nanmean(exp3a.d,1:2),nanmean(exp3b.d,1:2)],[nanstd(nanmean(exp3a.d,2))./sqrt(size(exp3a.d,1)),nanstd(nanmean(exp3b.d,2))./sqrt(size(exp3b.d,1))],'k' );hold on
he.LineStyle = 'none';he.CapSize = 12;
xticklabels({'Exp3a','Exp3b'});
ylabel('Distance(m)');grid off;box off;ylim([0 0.8]);
set(gcf,'position',[0,0,600,400]);
figure
bar(1:2,[nanmean(exp3a.o,1:2),nanmean(exp3b.o,1:2)],0.6,'FaceColor',[.5 .5 .5]);hold on
he = errorbar(1:2,[nanmean(exp3a.o,1:2),nanmean(exp3b.o,1:2)],[nanstd(nanmean(exp3a.o,2))./sqrt(size(exp3a.o,1)),nanstd(nanmean(exp3b.o,2))./sqrt(size(exp3b.o,1))],'k' );hold on
he.LineStyle = 'none';he.CapSize = 12;
ylabel('Relative Angle(бу)');grid off;box off;
xticklabels({'Exp3a','Exp3b'});ylim([0 20]);
set(gcf,'position',[0,0,600,400]);
