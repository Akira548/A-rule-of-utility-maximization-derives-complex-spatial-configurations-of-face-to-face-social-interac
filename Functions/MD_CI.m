function [MD,CI,SD] = MD_CI(Data,isdisp)

    if nargin < 2
        isdisp = 1;
    end
    MD=nanmean(Data(:));
    SD=nanstd(Data(:));
    [~,~,CI]=ttest(Data(:));
% SEM = std(Data(:))./sqrt(Data(:));               % Standard Error
% ts = tinv([0.05  0.95],length(Data(:))-1);      % T-Score
% CI = mean(Data(:)) + ts.*SEM;  
    if isdisp
        disp(['M:' num2str(round(MD,2)) '   95%CI = [' num2str(round(CI(1),2)) ', ' num2str(round(CI(2),2)) ']'  '   SD:' num2str(round(SD,2))]);
    end

end

