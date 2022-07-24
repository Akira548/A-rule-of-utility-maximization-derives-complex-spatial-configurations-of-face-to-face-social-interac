function Annotation = Idiaps_Poster_ImportAnnotation(filename, startRow, endRow)
%IMPORTFILE ���ı��ļ��е���ֵ������Ϊ�����롣
%   ANNOTATION1A = IMPORTFILE(FILENAME) ��ȡ�ı��ļ� FILENAME ��Ĭ��ѡ����Χ�����ݡ�
%
%   ANNOTATION1A = IMPORTFILE(FILENAME, STARTROW, ENDROW) ��ȡ�ı��ļ� FILENAME ��
%   STARTROW �е� ENDROW ���е����ݡ�
%
% Example:
%   annotation1A = importfile('annotation_1A.txt', 1, 50);
%
%    ������� TEXTSCAN��

% �� MATLAB �Զ������� 2020/06/02 11:16:18

%% ��ʼ��������
delimiter = {') ('};
if nargin<=2
    startRow = 1;
    endRow = inf;
end

%% ÿ���ı��еĸ�ʽ:
%   ��1: �ı� (%s)
%	��2: �ı� (%s)
%   ��3: �ı� (%s)
%	��4: �ı� (%s)
%   ��5: �ı� (%s)
%	��6: �ı� (%s)
%   ��7: �ı� (%s)
%	��8: �ı� (%s)
%   ��9: �ı� (%s)
%	��10: �ı� (%s)
%   ��11: �ı� (%s)
%	��12: �ı� (%s)
%   ��13: �ı� (%s)
%	��14: �ı� (%s)
%   ��15: �ı� (%s)
%	��16: �ı� (%s)
%   ��17: �ı� (%s)
%	��18: �ı� (%s)
%   ��19: �ı� (%s)
%	��20: �ı� (%s)
%   ��21: �ı� (%s)
%	��22: �ı� (%s)
% �й���ϸ��Ϣ������� TEXTSCAN �ĵ���
formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';

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

%% ���޷���������ݽ��еĺ���
% �ڵ��������δӦ���޷���������ݵĹ�����˲�����������롣Ҫ�����������޷���������ݵĴ��룬�����ļ���ѡ���޷������Ԫ����Ȼ���������ɽű���

%% ת��
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
