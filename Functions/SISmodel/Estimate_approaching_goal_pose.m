function q_goal = Estimate_approaching_goal_pose(p_i)
q_goal = [];
if size(p_i,1)==1
    gk = [];
else
    [r,cx,cy]=SIS_circle_fit(p_i);
    gk = [cx,cy,r];
end
A_out = Estimate_approaching_areas(gk,p_i);
if ~isempty(A_out)
    for j = 1:length(A_out)
        Lj=length(A_out{j});
        q_goal(j,1:2) = A_out{j}(round(Lj/2),:);
        if size(p_i,1)==1
            q_goal(j,3) = atan2(p_i(2)-q_goal(j,2), p_i(1)-q_goal(j,1));
        else
            q_goal(j,3) = atan2(gk(2)-q_goal(j,2), gk(1)-q_goal(j,1));
        end
    end
end

