function trackinghp = import_SubData_CocktailParty(filename, startRow, endRow)
%IMPORTFILE 将文本文件中的数值数据作为矩阵导入。
%   TRACKINGHP = IMPORTFILE(FILENAME) 读取文本文件 FILENAME 中默认选定范围的数据。
%
%   TRACKINGHP = IMPORTFILE(FILENAME, STARTROW, ENDROW) 读取文本文件 FILENAME 的
%   STARTROW 行到 ENDROW 行中的数据。
%
% Example:
%   trackinghp = importfile('tracking_hp.log', 1, 24875);
%
%    另请参阅 TEXTSCAN。

% 由 MATLAB 自动生成于 2020/01/10 15:20:11

%% 初始化变量。
if nargin<=2
    startRow = 1;
    endRow = inf;
end

%% 将数据列作为文本读取:
% 有关详细信息，请参阅 TEXTSCAN 文档。
formatSpec = '%17s%7s%5s%5s%6s%7s%5s%5s%6s%7s%5s%5s%6s%7s%5s%5s%6s%7s%5s%5s%6s%7s%5s%5s%6s%s%[^\n\r]';

%% 打开文本文件。
fileID = fopen(filename,'r');

%% 根据格式读取数据列。
% 该调用基于生成此代码所用的文件的结构。如果其他文件出现错误，请尝试通过导入工具重新生成代码。
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
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

for col=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25]
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
rawNumericColumns = raw(:, [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25]);
rawStringColumns = string(raw(:, 26));


%% 将非数值元胞替换为 NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),rawNumericColumns); % 查找非数值元胞
rawNumericColumns(R) = {NaN}; % 替换非数值元胞

%% 确保包含 <undefined> 的任何文本都已正确转换为 <undefined> 分类值
idx = (rawStringColumns(:, 1) == "<undefined>");
rawStringColumns(idx, 1) = "";

%% 创建输出变量
trackinghp1 = table;
trackinghp1.Time = cell2mat(rawNumericColumns(:, 1));
trackinghp1.x = cell2mat(rawNumericColumns(:, 3));
trackinghp1.y = cell2mat(rawNumericColumns(:, 4));
trackinghp1.o = cell2mat(rawNumericColumns(:, 5));

trackinghp2 = table;
trackinghp2.Time = cell2mat(rawNumericColumns(:, 1));
trackinghp2.x = cell2mat(rawNumericColumns(:, 7));
trackinghp2.y = cell2mat(rawNumericColumns(:, 8));
trackinghp2.o = cell2mat(rawNumericColumns(:, 9));

trackinghp3 = table;
trackinghp3.Time = cell2mat(rawNumericColumns(:, 1));
trackinghp3.x = cell2mat(rawNumericColumns(:, 11));
trackinghp3.y = cell2mat(rawNumericColumns(:, 12));
trackinghp3.o = cell2mat(rawNumericColumns(:, 13));

trackinghp4 = table;
trackinghp4.Time = cell2mat(rawNumericColumns(:, 1));
trackinghp4.x = cell2mat(rawNumericColumns(:, 15));
trackinghp4.y = cell2mat(rawNumericColumns(:, 16));
trackinghp4.o = cell2mat(rawNumericColumns(:, 17));

trackinghp5 = table;
trackinghp5.Time = cell2mat(rawNumericColumns(:, 1));
trackinghp5.x = cell2mat(rawNumericColumns(:, 19));
trackinghp5.y = cell2mat(rawNumericColumns(:, 20));
trackinghp5.o = cell2mat(rawNumericColumns(:, 21));

trackinghp6 = table;
trackinghp6.Time = cell2mat(rawNumericColumns(:, 1));
trackinghp6.x = cell2mat(rawNumericColumns(:, 23));
trackinghp6.y = cell2mat(rawNumericColumns(:, 24));
trackinghp6.o = cell2mat(rawNumericColumns(:, 25));

trackinghp{1} = trackinghp1;
trackinghp{2} = trackinghp2;
trackinghp{3} = trackinghp3;
trackinghp{4} = trackinghp4;
trackinghp{5} = trackinghp5;
trackinghp{6} = trackinghp6;



