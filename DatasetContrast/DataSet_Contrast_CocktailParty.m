%% Import Data [DataSet:'CocktailParty']
clc;clear all;warning off;
DataSet = 'CocktailParty';
tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));
addpath(genpath('../Functions'));

if (~exist(['./ImportedData' filesep DataSet filesep 'SubData.mat']))
    if ~exist(['./ImportedData/' DataSet])
        mkdir(['./ImportedData/' DataSet]);
    end
    if ~exist(['./Log/' DataSet])
        mkdir(['./Log/' DataSet]);
    end
    SubData = import_SubData_CocktailParty(['./Dataset_rawdata' filesep DataSet filesep 'tracking_hp.log']);
    Import_group = import_group(['./Dataset_rawdata' filesep DataSet filesep 'grouping.ref']);
    for scenei = 1:size(Import_group,1)
        index = find(SubData{1}.Time == Import_group.Time(scenei));
        for Gi = 1:6
            eval(sprintf('tempGname = Import_group.G%d(%d);',Gi,scenei));
            tempGnindex = isstrprop(tempGname{1},'digit');
            Group = tempGname{1}(find(double(tempGnindex & tempGname{1}~='0')));
            for gsizei = 1:length(Group)
                SubData{str2num(Group(gsizei))}.GT(index) = Gi;
            end
        end
    end
    save(['./ImportedData' filesep DataSet filesep 'SubData.mat'],'SubData');
    save(['./ImportedData' filesep DataSet filesep 'Import_group.mat'],'Import_group');
else
    load(['./ImportedData' filesep DataSet filesep  'SubData.mat']);
    load(['./ImportedData' filesep DataSet filesep  'Import_group.mat']);
end

%% DataSet Analysis [DataSet:'CocktailParty']
clc;Start_time = datestr(now,31);
counter = 1;DataSet = 'CocktailParty';
xmax = 8; xmin = 0; ymax = 8;ymin = 0; Dis_step = 0.2;Ori_step = pi/36;

% Testing from every scene and every person
for i = 1:size(Import_group,1)-1 % extract  samples every 6s (every 90 frames)
% for i = 1:10
    clear Subinform
    for j = 1:length(SubData)
        index = find(SubData{j}.Time == Import_group.Time(i));
        Subinform.X(i,j) = nanmean(SubData{j}.x(index-45:index+45));% mean X. timewindow = 6s
        Subinform.Y(i,j) = nanmean(SubData{j}.y(index-45:index+45));% mean Y. timewindow = 6s
        v = nanmean([cos(SubData{j}.o(index-45:index+45)),sin(SubData{j}.o(index-45:index+45))],1);
        Subinform.O(i,j) = cart2pol(v(:,1),v(:,2)); % mean orientation. timewindow = 6s
        Subinform.GT(i,j) = SubData{j}.GT(index);
    end
    for k = 1:max(Subinform.GT(i,:))
        Groupindex = find(Subinform.GT(i,:)==k);
        if (size(Groupindex,2)>2 )
            for iGroup = 1:size(Groupindex,2)
                I = Groupindex(iGroup);
                Members = Groupindex((Groupindex ~= I));
                currentI = (I);
                clear VH_x VH_y VH_o ;
                VH_x = Subinform.X(i,Members);
                VH_y = Subinform.Y(i,Members);
                VH_o = Subinform.O(i,Members);
                VH_pos{counter} = [VH_x',VH_y',VH_o'];
                SaliencyMap{:,counter} = CalculatePotentialFieldContrast([xmin,xmax;ymin,ymax],Dis_step,Ori_step,[VH_x',VH_y',VH_o'],0,1);
                pos_I(:,:,counter) = [Subinform.X(i,I),Subinform.Y(i,I)];
                ori_I(counter) = Subinform.O(i,I);
                
                counter = counter+1;
            end
        end
    end
    clc;disp([ 'now:   ', num2str(i) '/319' ]);
end

%%%
disp('Analyzing...');
close all

save(['./DataSetAnaLysis/' DataSet filesep DataSet]);
DoP = Contrast_o_Distance(SaliencyMap,pos_I,[xmin,xmax;ymin,ymax],Dis_step,VH_pos);title(DataSet);
saveas(gca,['./DataSetAnaLysis/' DataSet filesep  'DeviationDis.eps']);ylim([0 1.2]);
saveas(gca,['./DataSetAnaLysis/' DataSet filesep  'DeviationDis.png']);close all
DoO = Contrast_Angle(SaliencyMap,pos_I,ori_I,[xmin,xmax;ymin,ymax],Dis_step,VH_pos);title(DataSet);
saveas(gca,['./DataSetAnaLysis/' DataSet filesep  'DeviationAngle.png']);
saveas(gca,['./DataSetAnaLysis/' DataSet filesep  'DeviationAngle.eps']);close all

csvwrite('./DataSetAnaLysis/C_DoP.csv',[[1,2,3];DoP]);
csvwrite('./DataSetAnaLysis/C_DoO.csv',[[1,2,3];DoO]);

% DataAnalysis_odis_scatter(SaliencyMap,pos_I,[xmin,xmax;ymin,ymax],Dis_step,[0,3]);title(DataSet);
% DataAnalysis_Angle_scatter(SaliencyMap,pos_I,ori_I,[xmin,xmax;ymin,ymax],Dis_step);

disp(['Time: ' Start_time '   to   ' datestr(now,31)]);
disp('done!');