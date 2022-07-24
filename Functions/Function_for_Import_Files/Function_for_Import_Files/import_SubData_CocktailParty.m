function trackinghp = import_SubData_CocktailParty(filename, startRow, endRow)
%IMPORTFILE ���ı��ļ��е���ֵ������Ϊ�����롣
%   TRACKINGHP = IMPORTFILE(FILENAME) ��ȡ�ı��ļ� FILENAME ��Ĭ��ѡ����Χ�����ݡ�
%
%   TRACKINGHP = IMPORTFILE(FILENAME, STARTROW, ENDROW) ��ȡ�ı��ļ� FILENAME ��
%   STARTROW �е� ENDROW ���е����ݡ�
%
% Example:
%   trackinghp = importfile('tracking_hp.log', 1, 24875);
%
%    ������� TEXTSCAN��

% �� MATLAB �Զ������� 2020/01/10 15:20:11

%% ��ʼ��������
if nargin<=2
    startRow = 1;
    endRow = inf;
end

%% ����������Ϊ�ı���ȡ:
% �й���ϸ��Ϣ������� TEXTSCAN �ĵ���
formatSpec = '%17s%7s%5s%5s%6s%7s%5s%5s%6s%7s%5s%5s%6s%7s%5s%5s%6s%7s%5s%5s%6s%7s%5s%5s%6s%s%[^\n\r]';

%% ���ı��ļ���
fileID = fopen(filename,'r');

%% ���ݸ�ʽ��ȡ�����С�
% �õ��û������ɴ˴������õ��ļ��Ľṹ����������ļ����ִ����볢��ͨ�����빤���������ɴ��롣
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
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

for col=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25]
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
rawNumericColumns = raw(:, [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25]);
rawStringColumns = string(raw(:, 26));


%% ������ֵԪ���滻Ϊ NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),rawNumericColumns); % ���ҷ���ֵԪ��
rawNumericColumns(R) = {NaN}; % �滻����ֵԪ��

%% ȷ������ <undefined> ���κ��ı�������ȷת��Ϊ <undefined> ����ֵ
idx = (rawStringColumns(:, 1) == "<undefined>");
rawStringColumns(idx, 1) = "";

%% �����������
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



