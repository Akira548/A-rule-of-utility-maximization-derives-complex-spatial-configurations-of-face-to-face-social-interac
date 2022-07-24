function Aout = Estimate_approaching_areas(gk,p_i)
Aout = [];Aout1 = [];
r_g_step = 0.25;
alpha2 = deg2rad(140);
r_p_0 = 0.9;
Tc = 0.14;Tg = 0.14;

if size(p_i,1)==1
    [x,y]=pol2cart(deg2rad(0:.1:359),r_p_0*ones(size(deg2rad(0:.1:359))));
    A_in = [x'+p_i(1), y'+p_i(2)];r = r_p_0;
else
    [x,y]=pol2cart(deg2rad(0:.1:359),gk(3)*ones(size(deg2rad(0:.1:359))));
    A_in = [x'+gk(1), y'+gk(2)];r = gk(3);
end

A_in1 = A_in;
%% fov
while true
    for i = 1:size(A_in1,1)
        for j = 1:size(p_i,1)
            v2 = [A_in1(i,1)-p_i(j,1), A_in1(i,2)-p_i(j,2)];
            v1 = [cos(p_i(j,3)), sin(p_i(j,3))];
            voa(j)=acos(dot(v1,v2)/(norm(v1)*norm(v2)));
        end
        A_voa(i)=all(voa < alpha2/2);
    end
    Aout1 = find_Aout(A_voa, A_in1);
    %%
    Aout_i = 1;
    for m = 1:length(Aout1)
%         if length(Aout1)>1
            temp = Aout1{m};
%         else
%             temp = Aout1{1};
%         end
        tempouti = [];Aout_temp=[];
        for i = 1:size(temp,1)
            if size(p_i,1) == 1
                SII = basic_personal_space(temp(i,:),1,p_i);
                if SII < Tc
                    %                     tempout = [tempout; temp(i,:)];
                    tempouti(i) = 1;
                else
                    tempouti(i) = 0;
                end
            else
                SII = basic_personal_space(temp(i,:),1,p_i);
                SGI = Social_interaction_spce(temp(i,:),1,p_i);
                if (SII < Tc)&&(SGI < Tg)
                    %                     tempout = [tempout; temp(i,:)];
                    tempouti(i) = 1;
                else
                    tempouti(i) = 0;
                end
            end
        end
        if ~isempty(temp)
            Aout_temp = find_Aout(tempouti, temp);
        end
        if ~isempty(Aout_temp)
            if length(Aout_temp)>=1
                for Aout_temp_i = 1:length(Aout_temp)
                    Aout{Aout_i} = Aout_temp{Aout_temp_i};
                    Aout_i = Aout_i +1;
                end
            else
                    Aout{Aout_i} = Aout_temp{1};
                    Aout_i = Aout_i +1;
            end
            
        end
    end
    
    if isempty(Aout)
        r = r + r_g_step;
        if r>10
            Aout = [];
            warning('no result return');
            break
        end
        if size(p_i,1)==2
            [x,y]=pol2cart(deg2rad(0:359),r*ones(size(deg2rad(0:359))));
            A_in = [x', y'];
        else
            [x,y]=pol2cart(deg2rad(0:359),r*ones(size(deg2rad(0:359))));
            A_in = [x', y'];
        end
    else
        break;
    end
    
end

