%% Import Data [DataSet:'SLASA']
clc;clear all;warning off;
DataSet = 'SLASA';
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
    path = 'D:\work\dataset\SALAS\SALSA_Annotation_cpp\Annotation\salsa_cpp\geometryGT\';
    FilesName = dir([path,'*.csv']);
    for i = 1:length(FilesName)
        Annotation{i} = SALSA_ImportAnnotation([path FilesName(i).name]);
    end
    path = 'D:\work\dataset\SALAS\SALSA_Annotation_cpp\Annotation\salsa_cpp\';
    FilesName2 = dir([path,'*.csv']);
    GT_Annotation = SALSA_ImportGT([path FilesName2.name]);

    for i =1:size(FilesName,1)
        temp_SubData{i}.x = Annotation{i}.VarName2;
        temp_SubData{i}.y = Annotation{i}.VarName3;
        temp_SubData{i}.o = wrapToPi(Annotation{i}.VarName5 + Annotation{i}.VarName6);
        temp_SubData{i}.GT = zeros(size(temp_SubData{i}.o));
    end
    for i =1:size(Annotation{1}.VarName1,1)
        temp = find(GT_Annotation(:,1) == Annotation{1}.VarName1(i));
        for j = 1:length(temp)
            tempGroup = GT_Annotation(temp(j),2:end);
            tempGroup(isnan(tempGroup)) = [];
            for k = 1:length(tempGroup)
                temp_SubData{tempGroup(k)}.GT(i) = j;
            end
        end
    end
    for i = 1:length(temp_SubData)
        temp_min(i) = min(temp_SubData{i}.x);
    end
    for i =1:size(FilesName,1)
        for j =1:size(temp_SubData{1}.x,1)
            SubData{j}.x(i)  = temp_SubData{i}.x(j)-min(temp_min);
            SubData{j}.y(i)  = temp_SubData{i}.y(j);
            SubData{j}.o(i)  = temp_SubData{i}.o(j);
            SubData{j}.GT(i)  = temp_SubData{i}.GT(j);
        end
    end
    
    for i = 1:length(SubData)
        temp_max(i) = max(SubData{i}.x);
    end

    for i = 1:length(SubData)
        temp_max2(i) = max(SubData{i}.y);
    end
    
    save(['./ImportedData' filesep DataSet filesep 'SubData.mat'],'SubData');
else
    load(['./ImportedData'  filesep DataSet filesep 'SubData.mat']);
end

if(~exist(['./DataSetAnaLysis/' DataSet]))
    mkdir(['./DataSetAnaLysis/' DataSet]);
end

%% DataSet Analysis [DataSet:'SLASA']
clc;Start_time = datestr(now,31);
counter = 1;DataSet = 'SLASA';
xmax = 15; xmin = 0; ymax = 8;ymin = 0; Dis_step = 0.2;Ori_step = pi/36;

% Testing from every scene and every person
for i = 2:2:size(SubData,2)-1 % extra simple every 6s
    clear Subinform
    for j = 1:size(SubData{i}.x,2)
        Subinform.X(i,j) = nanmean([SubData{i-1}.x(j),SubData{i}.x(j),SubData{i+1}.x(j)]); % mean position. timewindow = 6s
        Subinform.Y(i,j) = nanmean([SubData{i-1}.y(j), SubData{i}.y(j), SubData{i+1}.y(j)]); % mean position. timewindow = 6s
        v = [cos(SubData{i-1}.o(j)),sin(SubData{i-1}.o(j))] + [cos(SubData{i}.o(j)),sin(SubData{i}.o(j))] + [cos(SubData{i+1}.o(j)),sin(SubData{i+1}.o(j))];
        Subinform.O(i,j) = cart2pol(v(:,1),v(:,2)); % mean orientation. timewindow = 6s
        Subinform.GroupSize(i,j) = size(unique(SubData{i}.GT),2);
    end
    for k = 0:max(SubData{i}.GT)
        Groupindex = find(SubData{i}.GT==k);
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
                SaliencyMap{:,counter}= CalculatePotentialFieldContrast([xmin,xmax;ymin,ymax],Dis_step,Ori_step,[VH_x',VH_y',VH_o'],0,1);
                pos_I(:,:,counter) = [Subinform.X(i,I),Subinform.Y(i,I)];
                ori_I(counter) = Subinform.O(i,I);
                counter = counter+1;
            end
        end
    end
    clc;disp([ 'now:   ', num2str(i) '/500' ]);
end

%%%
disp('Analyzing...');
close all
DoP = Contrast_o_Distance(SaliencyMap,pos_I,[xmin,xmax;ymin,ymax],Dis_step,VH_pos);title(DataSet);
saveas(gca,['./DataSetAnaLysis/' DataSet filesep  'DeviationDis.eps']);
saveas(gca,['./DataSetAnaLysis/' DataSet filesep  'DeviationDis.png']);close all
DoO = Contrast_Angle(SaliencyMap,pos_I,ori_I,[xmin,xmax;ymin,ymax],Dis_step,VH_pos);title(DataSet);
saveas(gca,['./DataSetAnaLysis/' DataSet filesep  'DeviationAngle.png']);
saveas(gca,['./DataSetAnaLysis/' DataSet filesep  'DeviationAngle.eps']);close all

save(['./DataSetAnaLysis/' DataSet filesep DataSet]);

csvwrite('./DataSetAnaLysis/S_DoP.csv',[[1,2,3];DoP]);
csvwrite('./DataSetAnaLysis/S_DoO.csv',[[1,2,3];DoO]);


disp(['Time: ' Start_time '   to   ' datestr(now,31)]);
disp('done!');

% DataAnalysis_odis_scatter(SaliencyMap,pos_I,[xmin,xmax;ymin,ymax],Dis_step,[0,3]);title(DataSet);
% DataAnalysis_Angle_scatter(SaliencyMap,pos_I,ori_I,[xmin,xmax;ymin,ymax],Dis_step);

% 
% props = java.lang.System.getProperties;
% props.setProperty('mail.smtp.auth','true');    
% setpref('Internet','SMTP_Username','hanjldx123@gmail.com');
% setpref('Internet','SMTP_Password','');
% % …Ë÷√” œ‰
% setpref('Internet','SMTP_Server','smtp.gmail.com');
% setpref('Internet','E_mail','hanjldx123@gmail.com');
% 
% sendmail('hanjldx123@gmail.com',DataSet, ...
%     DataSet,['./DataSetAnaLysis/' DataSet filesep  'DeviationDis.eps']);
% sendmail('hanjldx123@gmail.com','Hello from MATLAB!', ...
%     'Thanks for using sendmail.',['./DataSetAnaLysis/' DataSet filesep  'DeviationAngle.eps']);

