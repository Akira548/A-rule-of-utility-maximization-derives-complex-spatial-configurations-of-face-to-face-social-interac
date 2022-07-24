function Sub0 = importfile_Exp5a(filename, startRow, endRow)
%IMPORTFILE1 ���ı��ļ��е���ֵ������Ϊ�����롣
%   SUB0 = IMPORTFILE1(FILENAME) ��ȡ�ı��ļ� FILENAME ��Ĭ��ѡ����Χ�����ݡ�
%
%   SUB0 = IMPORTFILE1(FILENAME, STARTROW, ENDROW) ��ȡ�ı��ļ� FILENAME ��
%   STARTROW �е� ENDROW ���е����ݡ�
%
% Example:
%   Sub0 = importfile1('Sub_0.txt', 2, 73);
%
%    ������� TEXTSCAN��

% �� MATLAB �Զ������� 2020/01/06 15:44:55

%% ��ʼ��������
delimiter = '\t';
if nargin<=2
    startRow = 2;
    endRow = inf;
end

%% ����������Ϊ�ı���ȡ:
% �й���ϸ��Ϣ������� TEXTSCAN �ĵ���
formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';

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

for col=[1,2,3,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23]
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


%% �����ݲ��Ϊ��ֵ���ַ����С�
rawNumericColumns = raw(:, [1,2,3,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23]);
rawStringColumns = string(raw(:, 4));


%% ȷ������ <undefined> ���κ��ı�������ȷת��Ϊ <undefined> ����ֵ
idx = (rawStringColumns(:, 1) == "<undefined>");
rawStringColumns(idx, 1) = "";

%% �����������
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

