function Sub1 = importfile_1c(filename, dataLines)
%IMPORTFILE 从文本文件中导入数据
%  SUB1 = IMPORTFILE(FILENAME)读取文本文件 FILENAME 中默认选定范围的数据。  以表形式返回数据。
%
%  SUB1 = IMPORTFILE(FILE, DATALINES)按指定行间隔读取文本文件 FILENAME
%  中的数据。对于不连续的行间隔，请将 DATALINES 指定为正整数标量或 N×2 正整数标量数组。
%
%  示例:
%  Sub1 = importfile("C:\Users\Hank\OneDrive\Work\MATLAB\Interaction\RawData\Results\Interaction_1c\Sub_1.txt", [2, Inf]);
%
%  另请参阅 READTABLE。
%
% 由 MATLAB 于 2020-10-30 18:57:43 自动生成

%% 输入处理

% 如果不指定 dataLines，请定义默认范围
if nargin < 2
    dataLines = [2, Inf];
end

%% 设置导入选项
opts = delimitedTextImportOptions("NumVariables", 12);

% 指定范围和分隔符
opts.DataLines = dataLines;
opts.Delimiter = "\t";

% 指定列名称和类型
opts.VariableNames = ["Experiment", "Sub_Num", "Age1", "Sex1", "Age2", "Sex2", "Handness", "Theta1", "Theta2", "Distance", "player1", "player2"];
opts.VariableTypes = ["char", "double", "double", "categorical", "double", "categorical", "double", "double", "double", "double", "double", "double"];
opts = setvaropts(opts, 1, "WhitespaceRule", "preserve");
opts = setvaropts(opts, 2, "TrimNonNumeric", true);
opts = setvaropts(opts, 2, "ThousandsSeparator", ",");
opts = setvaropts(opts, [1, 4, 6], "EmptyFieldRule", "auto");
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% 导入数据
Sub1 = readtable(filename, opts);

end