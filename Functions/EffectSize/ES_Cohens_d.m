function d = ES_Cohens_d(data1,data2)
%%% Cohen's d = (M2 - M1) ? SDpooled
%%% SDpooled = ¡Ì((SD1^2 + SD2^2) ? 2)
M1 = nanmean(data1);
M2 = nanmean(data2);
SD1 = nanstd(data1);
SD2 = nanstd(data2);
SDpooled = sqrt((SD1.^2 + SD2.^2)/2);
d = abs((M2 - M1)./SDpooled);
end

