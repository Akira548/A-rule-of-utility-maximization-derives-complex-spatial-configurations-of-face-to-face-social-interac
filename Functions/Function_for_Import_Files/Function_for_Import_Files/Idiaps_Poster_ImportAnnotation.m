function Annotation = Idiaps_Poster_ImportAnnotation(filename, startRow, endRow)
%IMPORTFILE 将文本文件中的数值数据作为矩阵导入。
%   ANNOTATION1A = IMPORTFILE(FILENAME) 读取文本文件 FILENAME 中默认选定范围的数据。
%
%   ANNOTATION1A = IMPORTFILE(FILENAME, STARTROW, ENDROW) 读取文本文件 FILENAME 的
%   STARTROW 行到 ENDROW 行中的数据。
%
% Example:
%   annotation1A = importfile('annotation_1A.txt', 1, 50);
%
%    另请参阅 TEXTSCAN。

% 由 MATLAB 自动生成于 2020/06/02 11:16:18

%% 初始化变量。
delimiter = {') ('};
if nargin<=2
    startRow = 1;
    endRow = inf;
end

%% 每个文本行的格式:
%   列1: 文本 (%s)
%	列2: 文本 (%s)
%   列3: 文本 (%s)
%	列4: 文本 (%s)
%   列5: 文本 (%s)
%	列6: 文本 (%s)
%   列7: 文本 (%s)
%	列8: 文本 (%s)
%   列9: 文本 (%s)
%	列10: 文本 (%s)
%   列11: 文本 (%s)
%	列12: 文本 (%s)
%   列13: 文本 (%s)
%	列14: 文本 (%s)
%   列15: 文本 (%s)
%	列16: 文本 (%s)
%   列17: 文本 (%s)
%	列18: 文本 (%s)
%   列19: 文本 (%s)
%	列20: 文本 (%s)
%   列21: 文本 (%s)
%	列22: 文本 (%s)
% 有关详细信息，请参阅 TEXTSCAN 文档。
formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';

%% 打开文本文件。
fileID = fopen(filename,'r');

%% 根据格式读取数据列。
% 该调用基于生成此代码所用的文件的结构。如果其他文件出现错误，请尝试通过导入工具重新生成代码。
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% 关闭文本文件。
fclose(fileID);

%% 对无法导入的数据进行的后处理。
% 在导入过程中未应用无法导入的数据的规则，因此不包括后处理代码。要生成适用于无法导入的数据的代码，请在文件中选择无法导入的元胞，然后重新生成脚本。

%% 转化
tempannotation = [dataArray{1:end-1}];
tempannotation = strrep(tempannotation,'(','');
tempannotation = strrep(tempannotation,')','');
for i = 1:size(tempannotation,1)/5
    if(i == 1)
    Annotation.Positions = tempannotation(5*(i-1) +2,:);
    else
    Annotation.Positions = [Annotation.Positions; tempannotation(5*(i-1) +2,:)];
    end
    Annotation.Group{i} = strsplit(tempannotation(5*(i-1) +3,1));
end
