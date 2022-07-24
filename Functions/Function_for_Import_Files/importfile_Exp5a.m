function Sub0 = importfile_Exp5a(filename, startRow, endRow)
%IMPORTFILE1 将文本文件中的数值数据作为矩阵导入。
%   SUB0 = IMPORTFILE1(FILENAME) 读取文本文件 FILENAME 中默认选定范围的数据。
%
%   SUB0 = IMPORTFILE1(FILENAME, STARTROW, ENDROW) 读取文本文件 FILENAME 的
%   STARTROW 行到 ENDROW 行中的数据。
%
% Example:
%   Sub0 = importfile1('Sub_0.txt', 2, 73);
%
%    另请参阅 TEXTSCAN。

% 由 MATLAB 自动生成于 2020/01/06 15:44:55

%% 初始化变量。
delimiter = '\t';
if nargin<=2
    startRow = 2;
    endRow = inf;
end

%% 将数据列作为文本读取:
% 有关详细信息，请参阅 TEXTSCAN 文档。
formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';

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

%% 将包含数值文本的列内容转换为数值。
% 将非数值文本替换为 NaN。
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

for col=[1,2,3,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23]
    % 将输入元胞数组中的文本转换为数值。已将非数值文本替换为 NaN。
    rawData = dataArray{col};
    for row=1:size(rawData, 1)
        % 创建正则表达式以检测并删除非数值前缀和后缀。
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData(row), regexstr, 'names');
            numbers = result.numbers;
            
            % 在非千位位置中检测到逗号。
            invalidThousandsSeparator = false;
            if numbers.contains(',')
                thousandsRegExp = '^[-/+]*\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(numbers, thousandsRegExp, 'once'))
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % 将数值文本转换为数值。
            if ~invalidThousandsSeparator
                numbers = textscan(char(strrep(numbers, ',', '')), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch
            raw{row, col} = rawData{row};
        end
    end
end


%% 将数据拆分为数值和字符串列。
rawNumericColumns = raw(:, [1,2,3,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23]);
rawStringColumns = string(raw(:, 4));


%% 确保包含 <undefined> 的任何文本都已正确转换为 <undefined> 分类值
idx = (rawStringColumns(:, 1) == "<undefined>");
rawStringColumns(idx, 1) = "";

%% 创建输出变量
Sub0 = table;
Sub0.Experiment = cell2mat(rawNumericColumns(:, 1));
Sub0.Sub_Num = cell2mat(rawNumericColumns(:, 2));
Sub0.Age = cell2mat(rawNumericColumns(:, 3));
Sub0.Sex = categorical(rawStringColumns(:, 1));
Sub0.Handness = cell2mat(rawNumericColumns(:, 4));
Sub0.DistanceA = cell2mat(rawNumericColumns(:, 5));
Sub0.ORientaitonA = cell2mat(rawNumericColumns(:, 6));
Sub0.ORientaitonB = cell2mat(rawNumericColumns(:, 7));
Sub0.StartPoint = cell2mat(rawNumericColumns(:, 8));
Sub0.EndPoint = cell2mat(rawNumericColumns(:, 9));
Sub0.DistanceB = cell2mat(rawNumericColumns(:, 10));
Sub0.Origin_x = cell2mat(rawNumericColumns(:, 11));
Sub0.Origin_y = cell2mat(rawNumericColumns(:, 12));
Sub0.Origin_z = cell2mat(rawNumericColumns(:, 13));
Sub0.OriginOrtation_x = cell2mat(rawNumericColumns(:, 14));
Sub0.OriginOrtation_y = cell2mat(rawNumericColumns(:, 15));
Sub0.OriginOrtation_z = cell2mat(rawNumericColumns(:, 16));
Sub0.x = cell2mat(rawNumericColumns(:, 17));
Sub0.y = cell2mat(rawNumericColumns(:, 18));
Sub0.z = cell2mat(rawNumericColumns(:, 19));
Sub0.Ortation_x = cell2mat(rawNumericColumns(:, 20));
Sub0.Ortation_y = cell2mat(rawNumericColumns(:, 21));
Sub0.Ortation_z = cell2mat(rawNumericColumns(:, 22));

