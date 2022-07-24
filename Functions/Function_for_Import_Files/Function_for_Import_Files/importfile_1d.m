function Sub11 = importfile_1d(filename, startRow, endRow)
%% 初始化变量。
delimiter = '\t';
if nargin<=2
    startRow = 2;
    endRow = inf;
end

%% 将数据列作为文本读取:
% 有关详细信息，请参阅 TEXTSCAN 文档。
formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';

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

for col=[2,3,4,5,8,9,10,11,12,13,14,15]
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

% 使用指定的日期格式将带有日期的列内容转换为 MATLAB 日期时间。
try
    dates{1} = datetime(dataArray{1}, 'Format', 'dd/MM/yyyy', 'InputFormat', 'dd/MM/yyyy');
catch
    try
        % 处理用引号引起来的日期
        dataArray{1} = cellfun(@(x) x(2:end-1), dataArray{1}, 'UniformOutput', false);
        dates{1} = datetime(dataArray{1}, 'Format', 'dd/MM/yyyy', 'InputFormat', 'dd/MM/yyyy');
    catch
        dates{1} = repmat(datetime([NaN NaN NaN]), size(dataArray{1}));
    end
end

dates = dates(:,1);

%% 将数据拆分为数值和字符串列。
rawNumericColumns = raw(:, [2,3,4,5,8,9,10,11,12,13,14,15]);
rawStringColumns = string(raw(:, [6,7]));


%% 确保包含 <undefined> 的任何文本都已正确转换为 <undefined> 分类值
for catIdx = [1,2]
    idx = (rawStringColumns(:, catIdx) == "<undefined>");
    rawStringColumns(idx, catIdx) = "";
end

%% 创建输出变量
Sub11 = table;
Sub11.Date = dates{:, 1};
Sub11.Trial_No = cell2mat(rawNumericColumns(:, 1));
Sub11.Block = cell2mat(rawNumericColumns(:, 2));
Sub11.Sub_Num = cell2mat(rawNumericColumns(:, 3));
Sub11.Age = cell2mat(rawNumericColumns(:, 4));
Sub11.Sex = categorical(rawStringColumns(:, 1));
Sub11.Handedness = categorical(rawStringColumns(:, 2));
Sub11.A_Angle = cell2mat(rawNumericColumns(:, 5));
Sub11.VisualAngle = cell2mat(rawNumericColumns(:, 6));
Sub11.Distance = cell2mat(rawNumericColumns(:, 7));
Sub11.avatar_Sex = cell2mat(rawNumericColumns(:, 8));
Sub11.avatar_No = cell2mat(rawNumericColumns(:, 9));
Sub11.VisualAngleDirection = cell2mat(rawNumericColumns(:, 10));
Sub11.Answer = cell2mat(rawNumericColumns(:, 11));
Sub11.RT = cell2mat(rawNumericColumns(:, 12));

