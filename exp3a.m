%% ImportData
clc;close all;clear all;
ExpName = 'Exp3a';
tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));
addpath(genpath('./Functions'));

FilesName = dir(['.\Raw Data\' ExpName filesep '*.txt']);
for subi = 1:size(FilesName,1)
    tempdata = importfile_Exp3a(['.\Raw Data\' ExpName filesep FilesName(subi).name]);
    SubData.Distance(subi,:) = tempdata.Distance;
    SubData.z(subi,:) = tempdata.Distance - tempdata.z;
    SubData.x(subi,:) = tempdata.x;
    SubData.o(subi,:) = tempdata.Ortation_y;
end

%% bar 
close all
conditon_d = SubData.Distance;
[Potentialmap, max_Ori1, AR2_record] = CalculatePotentialField([-6,6;-6,6],0.1,pi/36,[0, 0, 0],0,1);
[max_x,max_y] = find(Potentialmap' == max(Potentialmap',[],'all'));

for subi = 1:size(SubData.Distance,1)
    pos_all(subi,1,:) = [mean(SubData.z(subi,(conditon_d(subi,:) == 3))),mean(SubData.x(subi,(conditon_d(subi,:) == 3)))];
    pos_all(subi,2,:) = [mean(SubData.z(subi,(conditon_d(subi,:) == 5))),mean(SubData.x(subi,(conditon_d(subi,:) == 5)))];
    pos_all(subi,3,:) = [mean(SubData.z(subi,(conditon_d(subi,:) == 7))),mean(SubData.x(subi,(conditon_d(subi,:) == 7)))];
    ori_all(subi,1,:) = mean([cos(pi-(((wrapToPi(deg2rad(SubData.o(subi,(conditon_d(subi,:) == 3)))) ))))', sin(pi-(((wrapToPi(deg2rad(SubData.o(subi,(conditon_d(subi,:) == 3)))) ))))']);
    ori_all(subi,2,:) = mean([cos(pi-(((wrapToPi(deg2rad(SubData.o(subi,(conditon_d(subi,:) == 5)))) ))))', sin(pi-(((wrapToPi(deg2rad(SubData.o(subi,(conditon_d(subi,:) == 5)))) ))))']);
    ori_all(subi,3,:) = mean([cos(pi-(((wrapToPi(deg2rad(SubData.o(subi,(conditon_d(subi,:) == 7)))) ))))', sin(pi-(((wrapToPi(deg2rad(SubData.o(subi,(conditon_d(subi,:) == 7)))) ))))']);
end
exp3a.d = deviation_of_position(pos_all(:,:,1),pos_all(:,:,2),Potentialmap,[-6,6;-6,6],0.1);
exp3a.o = rad2deg(deviation_of_orientation(pos_all(:,:,1),pos_all(:,:,2),ori_all(:,:,1),ori_all(:,:,2),max_Ori1,[-6,6;-6,6],0.1));

% for subi = 1:size(SubData.Distance,1)
%     deviation_dis3(subi) = (sqrt((-mean(SubData.x(subi,(conditon_d(subi,:) == 3)))+((max_x - 1)/10 - 6)).^2 + (-mean(SubData.z(subi,(conditon_d(subi,:) == 3)))+((max_y - 1)/10 - 6)).^2));
%     deviation_dis5(subi) = (sqrt((-mean(SubData.x(subi,(conditon_d(subi,:) == 5)))+((max_x - 1)/10 - 6)).^2 + (-mean(SubData.z(subi,(conditon_d(subi,:) == 5)))+((max_y - 1)/10 - 6)).^2));
%     deviation_dis7(subi) = (sqrt((-mean(SubData.x(subi,(conditon_d(subi,:) == 7)))+((max_x - 1)/10 - 6)).^2 + (-mean(SubData.z(subi,(conditon_d(subi,:) == 7)))+((max_y - 1)/10 - 6)).^2));
% 
%  for j = 1:5
%     o3=-pi/2+(((wrapToPi(deg2rad(SubData.o(subi,(conditon_d(subi,:) == 3)))) )));
%     o5=-pi/2+(((wrapToPi(deg2rad(SubData.o(subi,(conditon_d(subi,:) == 5)))) )));
%     o7=-pi/2+(((wrapToPi(deg2rad(SubData.o(subi,(conditon_d(subi,:) == 7)))) )));
%     ot3=[round(SubData.x(subi,(conditon_d(subi,:) == 3))*10+61);round(SubData.z(subi,(conditon_d(subi,:) == 3))*10+61)];
%     ot5=[round(SubData.x(subi,(conditon_d(subi,:) == 5))*10+61);round(SubData.z(subi,(conditon_d(subi,:) == 5))*10+61)];
%     ot7=[round(SubData.x(subi,(conditon_d(subi,:) == 7))*10+61);round(SubData.z(subi,(conditon_d(subi,:) == 7))*10+61)];
%     temp_o_dis3(subi,j) = abs(o3(j)-max_Ori1(ot3(1,j),ot3(2,j)));
%     temp_o_dis5(subi,j) = abs(o5(j)-max_Ori1(ot5(1,j),ot5(2,j)));
%     temp_o_dis7(subi,j) = abs(o7(j)-max_Ori1(ot7(1,j),ot7(2,j)));
% end
%     o_dis3(subi) = mean(temp_o_dis3(subi,:),2);
%     o_dis5(subi) = mean(temp_o_dis5(subi,:),2);
%     o_dis7(subi) = mean(temp_o_dis7(subi,:),2);
%     
%     x_mean(subi,:) = [mean(SubData.x(subi,(conditon_d(subi,:) == 3))),mean(SubData.x(subi,(conditon_d(subi,:) == 5))),mean(SubData.x(subi,(conditon_d(subi,:) == 7)))];
%     y_mean(subi,:) = [mean(SubData.z(subi,(conditon_d(subi,:) == 3))),mean(SubData.z(subi,(conditon_d(subi,:) == 5))),mean(SubData.z(subi,(conditon_d(subi,:) == 7)))];
% end
% 
% exp3a.d = [deviation_dis3',deviation_dis5',deviation_dis7'];
% exp3a.o = rad2deg([o_dis3',o_dis5',o_dis7']);

MD_CI(exp3a.d,1);
[tbl2,stats2]=simple_mixed_anova(exp3a.d);
partial_Eta_squared = ES_Eta_squared(tbl2,2);

% deviation distance analysis
MD_CI(exp3a.o,1);
[tbl3,stats3]=simple_mixed_anova(exp3a.o);
partial_Eta_squared = ES_Eta_squared(tbl3,2);

%%
Potentialmaptest = CalculatePotentialField([-1,7;-4,4],0.05,pi/180,[0, 0, 0],0);
scatter3(mean(SubData.z,2),mean(SubData.x,2),ones(size(mean(SubData.x,2)))*1.1,200,'.k');hold on
scatter3(mean(SubData.z,2),mean(SubData.x,2),ones(size(mean(SubData.x,2)))*1.1,200,'.k');hold on
colorbar
% box off;
[s_x,s_y] = meshgrid(-1:0.05:7,-4:0.05:4);
surf(s_x,s_y,(Potentialmaptest)');hold on;shading interp;
axis equal
xlim([-1,7]);ylim([-4,4]);

view(90,90);
% set(gca,'DataAspectRatio',[1 1 .7]);
% lightangle(gca,0,20);material([0.8 1 0.6]);
c = colorbar;
% c.Label.String = 'Social Potential ';
box off;axis off;grid off
colormap(colorbar_BuOr());
%%
[Potentialmaptest, max_Ori1, ori_distribution] = CalculatePotentialField([-1,7;-4,4],0.05,pi/180,[0, 0, 0],0,1);
scatter(mean(SubData.z,1:2),mean(SubData.x,1:2),ones(size(mean(SubData.x,1:2)))*1.1,50,'xk');hold on
ix_seq = -1:0.05:7;
iy_seq = -4:0.05:4;
Ori_step=pi/180;

o=[round(mean(SubData.z,1:2)*20+21);round(mean(SubData.x,1:2)*20+81)];
temp_Rho=squeeze(ori_distribution(o(1),o(2),:));
for ix = 1:length(ix_seq)
    for iy = 1:length(iy_seq)
        [Theta,Rho]= cart2pol((ix_seq(ix)-mean(SubData.z,1:2)),(iy_seq(iy)-mean(SubData.x,1:2)));
        average_v= mean([cosd(SubData.o(:)),sind(SubData.o(:))],1);
        [average_ori,~]=  cart2pol(average_v(1),average_v(2));
        [Theta,Rho]=  cart2pol((ix_seq(ix)-mean(SubData.z,1:2)),(iy_seq(iy)-mean(SubData.x,1:2)));
        ori_strength(ix,iy)=temp_Rho( round((Theta)/Ori_step + pi/Ori_step + 1) );
    end
end
surf(ix_seq,iy_seq,ori_strength(:,:)');hold on;shading interp;
view(90,90);
axis equal
xlim([-1,7]);ylim([-4,4]);caxis([0 max(Potentialmaptest,[],'all')]);
c = colorbar;
axis off;box off;grid off

plot([mean(SubData.z,1:2)+2*cos(pi+average_ori),mean(SubData.z,1:2)],[mean(SubData.x,1:2)+2*sin(pi+average_ori),mean(SubData.x,1:2)],'--k','linewidth',2);

%%
close all
Potential_p1 = CalculatePotentialField([-4,8;-2.5,2.5],0.05,pi/90,[0, 0, 0],0);

[max_index_y,max_index_x] =find(Potential_p1== max(Potential_p1,[],'all'));
% disp((max_index_y-1)*20 -2.5  - VA_pos(2));
[dx,dy]=gradient(Potential_p1');
[s_x,s_y] = meshgrid(-4:0.05:8,-2.5:0.05:2.5);
contour3(s_x,s_y,real(Potential_p1)',[0,0.1,0.2,0.3:0.05:0.6 0.625:0.025:0.8],'-k','LineWidth',.5);hold on
surf(s_x,s_y,(Potential_p1)');hold on;shading interp;hold on
view(50,45);
axis equal;
xlim([-4,8]);ylim([-2.5,2.5]);
set(gca,'ZDir','reverse','DataAspectRatio',[1 1 .3]);
material([0.65 1 0.6]);
axis off
box off;grid off
lightangle(gca,40,-50)
lighting gouraud
