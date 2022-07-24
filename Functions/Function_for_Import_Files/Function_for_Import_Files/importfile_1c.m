function Sub1 = importfile_1c(filename, dataLines)
%IMPORTFILE ���ı��ļ��е�������
%  SUB1 = IMPORTFILE(FILENAME)��ȡ�ı��ļ� FILENAME ��Ĭ��ѡ����Χ�����ݡ�  �Ա���ʽ�������ݡ�
%
%  SUB1 = IMPORTFILE(FILE, DATALINES)��ָ���м����ȡ�ı��ļ� FILENAME
%  �е����ݡ����ڲ��������м�����뽫 DATALINES ָ��Ϊ������������ N��2 �������������顣
%
%  ʾ��:
%  Sub1 = importfile("C:\Users\Hank\OneDrive\Work\MATLAB\Interaction\RawData\Results\Interaction_1c\Sub_1.txt", [2, Inf]);
%
%  ������� READTABLE��
%
% �� MATLAB �� 2020-10-30 18:57:43 �Զ�����

%% ���봦��

% �����ָ�� dataLines���붨��Ĭ�Ϸ�Χ
if nargin < 2
    dataLines = [2, Inf];
end

%% ���õ���ѡ��
opts = delimitedTextImportOptions("NumVariables", 12);

% ָ����Χ�ͷָ���
opts.DataLines = dataLines;
opts.Delimiter = "\t";

% ָ�������ƺ�����
opts.VariableNames = ["Experiment", "Sub_Num", "Age1", "Sex1", "Age2", "Sex2", "Handness", "Theta1", "Theta2", "Distance", "player1", "player2"];
opts.VariableTypes = ["char", "double", "double", "categorical", "double", "categorical", "double", "double", "double", "double", "double", "double"];
opts = setvaropts(opts, 1, "WhitespaceRule", "preserve");
opts = setvaropts(opts, 2, "TrimNonNumeric", true);
opts = setvaropts(opts, 2, "ThousandsSeparator", ",");
opts = setvaropts(opts, [1, 4, 6], "EmptyFieldRule", "auto");
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% ��������
Sub1 = readtable(filename, opts);

end