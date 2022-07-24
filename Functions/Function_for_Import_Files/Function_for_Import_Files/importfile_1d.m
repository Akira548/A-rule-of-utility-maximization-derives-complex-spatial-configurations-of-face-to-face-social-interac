function Sub11 = importfile_1d(filename, startRow, endRow)
%% ��ʼ��������
delimiter = '\t';
if nargin<=2
    startRow = 2;
    endRow = inf;
end

%% ����������Ϊ�ı���ȡ:
% �й���ϸ��Ϣ������� TEXTSCAN �ĵ���
formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';

%% ���ı��ļ���
fileID = fopen(filename,'r');

%% ���ݸ�ʽ��ȡ�����С�
% �õ��û������ɴ˴������õ��ļ��Ľṹ����������ļ����ִ����볢��ͨ�����빤���������ɴ��롣
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% �ر��ı��ļ���
fclose(fileID);

%% ��������ֵ�ı���������ת��Ϊ��ֵ��
% ������ֵ�ı��滻Ϊ NaN��
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

for col=[2,3,4,5,8,9,10,11,12,13,14,15]
    % ������Ԫ�������е��ı�ת��Ϊ��ֵ���ѽ�����ֵ�ı��滻Ϊ NaN��
    rawData = dataArray{col};
    for row=1:size(rawData, 1)
        % ����������ʽ�Լ�Ⲣɾ������ֵǰ׺�ͺ�׺��
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData(row), regexstr, 'names');
            numbers = result.numbers;
            
            % �ڷ�ǧλλ���м�⵽���š�
            invalidThousandsSeparator = false;
            if numbers.contains(',')
                thousandsRegExp = '^[-/+]*\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(numbers, thousandsRegExp, 'once'))
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % ����ֵ�ı�ת��Ϊ��ֵ��
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

% ʹ��ָ�������ڸ�ʽ���������ڵ�������ת��Ϊ MATLAB ����ʱ�䡣
try
    dates{1} = datetime(dataArray{1}, 'Format', 'dd/MM/yyyy', 'InputFormat', 'dd/MM/yyyy');
catch
    try
        % ����������������������
        dataArray{1} = cellfun(@(x) x(2:end-1), dataArray{1}, 'UniformOutput', false);
        dates{1} = datetime(dataArray{1}, 'Format', 'dd/MM/yyyy', 'InputFormat', 'dd/MM/yyyy');
    catch
        dates{1} = repmat(datetime([NaN NaN NaN]), size(dataArray{1}));
    end
end

dates = dates(:,1);

%% �����ݲ��Ϊ��ֵ���ַ����С�
rawNumericColumns = raw(:, [2,3,4,5,8,9,10,11,12,13,14,15]);
rawStringColumns = string(raw(:, [6,7]));


%% ȷ������ <undefined> ���κ��ı�������ȷת��Ϊ <undefined> ����ֵ
for catIdx = [1,2]
    idx = (rawStringColumns(:, catIdx) == "<undefined>");
    rawStringColumns(idx, catIdx) = "";
end

%% �����������
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

