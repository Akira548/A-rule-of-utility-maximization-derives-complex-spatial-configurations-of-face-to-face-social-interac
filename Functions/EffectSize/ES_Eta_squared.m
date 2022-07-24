function Eta_squared = ES_Eta_squared(tbl,mode)
if nargin <2
    mode = 1;
end
switch mode 
    case 1
        %%% Partial eta2 squared/Eta_squared = SS_Treatment / SS_Total
        for rowi = 1:size(tbl,2)
            if strcmp(tbl{1,rowi}, 'SS')
                row = rowi;
            end
        end
        for coli = 1:size(tbl,1)
            if strcmp(tbl{coli,1}, 'Columns')
                SS_Treatment = tbl{coli,row};
            end
            if strcmp(tbl{coli,1}, 'Total')
                SS_Total = tbl{coli,row};
            end
        end
        Eta_squared = round(SS_Treatment / SS_Total,3);
        disp(['partial Eta squared:' num2str(Eta_squared)]);
    case 2
        %%% Partial eta2 = SSeffect / SSeffect + SSerror.
        counter = 1;
        for rowi = 3:size(tbl.Row,1)
            if strfind(tbl.Row{rowi},'(Intercept):')
                Eta_squared{counter} = round(tbl.SumSq(rowi)/(tbl.SumSq(rowi) + tbl.SumSq(rowi+1)),3);
                disp(['Eta_squared ' tbl.Row{rowi}(13:end) ' = ' num2str(Eta_squared{counter})  '  ,F = ' num2str(round(tbl.F(rowi),3)) '  ,p = ' num2str(round(tbl.pValue(rowi),3))  '  ,dF = (' num2str(tbl.DF(rowi)) ',' num2str(tbl.DF(rowi+1)) ')']);
                counter = counter+1;
            end
        end
end
end

