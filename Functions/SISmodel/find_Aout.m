function A_out1 = find_Aout(A_voa,A_in1)
starts = [];
ends = [];
A_out1 = [];m=1;
for i = 1:length(A_voa)
    if i == 1
        if(A_voa(i) && ~A_voa(end))
            starts = [starts,i];
        end
        if(~A_voa(i) && A_voa(end))
            ends = [ends,length(A_voa)];
        end
    else
        if(A_voa(i) && ~A_voa(i-1))
            starts = [starts,i];
        end
        if(~A_voa(i) && A_voa(i-1))
            ends = [ends,i-1];
        end
    end
end
if ~isempty(starts)
    if A_voa(end)==1 && (starts(1)>ends(1))
        Aout_index = [starts; [ends(2:end), ends(1)]];
    else
        Aout_index = [starts; ends];
    end
    for i = 1:size(Aout_index,2)
        if Aout_index(1,i)<=Aout_index(2,i)
            A_out1{m}  = A_in1(Aout_index(1,i):Aout_index(2,i),:);
            m = m+1;
        else
            A_out1{m}  = A_in1([Aout_index(1,i):size(A_voa,2),1:Aout_index(2,i)],:);
            m = m+1;
        end
    end
else
    A_out1{1}=A_in1;
end
end

