function Sub999 = importfile(filename, startRow, endRow)
%IMPORTFILE ���ı��ļ��е���ֵ������Ϊ�����롣
%   SUB999 = IMPORTFILE(FILENAME) ��ȡ�ı��ļ� FILENAME ��Ĭ��ѡ����Χ�����ݡ�
%
%   SUB999 = IMPORTFILE(FILENAME, STARTROW, ENDROW) ��ȡ�ı��ļ� FILENAME ��
%   STARTROW �е� ENDROW ���е����ݡ�
%
% Example:
%   Sub999 = importfile('Sub_999.txt', 2, 57);
%
%    ������� TEXTSCAN��

% �� MATLAB �Զ������� 2020/08/04 15:47:58

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

for col=[1,2,3,5,6,7,8,9,10,11,12,13,14,15]
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
rawNumericColumns = raw(:, [1,2,3,5,6,7,8,9,10,11,12,13,14,15]);
rawStringColumns = string(raw(:, 4));


%% ȷ������ <undefined> ���κ��ı�������ȷת��Ϊ <undefined> ����ֵ
idx = (rawStringColumns(:, 1) == "<undefined>");
rawStringColumns(idx, 1) = "";

%% �����������
Sub999 = table;
Sub999.Experiment = cell2mat(rawNumericColumns(:, 1));
Sub999.Sub_Num = cell2mat(rawNumericColumns(:, 2));
Sub999.Age = cell2mat(rawNumericColumns(:, 3));
Sub999.Sex = categorical(rawStringColumns(:, 1));
Sub999.Handness = cell2mat(rawNumericColumns(:, 4));
Sub999.Condition = cell2mat(rawNumericColumns(:, 5));
Sub999.Orientation = cell2mat(rawNumericColumns(:, 6));
Sub999.Distance = cell2mat(rawNumericColumns(:, 7));
Sub999.Origin_x = cell2mat(rawNumericColumns(:, 8));
Sub999.Origin_y = cell2mat(rawNumericColumns(:, 9));
Sub999.Origin_z = cell2mat(rawNumericColumns(:, 10));
Sub999.x = cell2mat(rawNumericColumns(:, 11));
Sub999.y = cell2mat(rawNumericColumns(:, 12));
Sub999.z = cell2mat(rawNumericColumns(:, 13));
Sub999.o = cell2mat(rawNumericColumns(:, 14));

