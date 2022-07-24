function [IPD,Ori1,Ori2] = CalculateInterperonalDistance(VH_pos)
% %%% Formula
 
% Theta1_Intensity = [];

%%%
dyadic=nchoosek(1:size(VH_pos,1),2);
for di = 1:size(dyadic,1)
    person1 = dyadic(di,1);person2 = dyadic(di,2);
    Orientation1 = wrapToPi(VH_pos(person1,3));
    Orientation2 = wrapToPi(VH_pos(person2,3));
    v1 = [VH_pos(person2,1)-VH_pos(person1,1),VH_pos(person2,2)-VH_pos(person1,2)];% person2 point to person1
    v2 = -v1;% person2 point to person1
    [Vector_Orientation1,Rho] = cart2pol(v1(:,1),v1(:,2));
    [Vector_Orientation2,~] = cart2pol(v2(:,1),v2(:,2));
    Theta1 = wrapToPi(Orientation1 - Vector_Orientation1);
    Theta2 = wrapToPi(Orientation2 - Vector_Orientation2);
    IPD(di)=Rho;
    Ori1(di)=Theta1;
    Ori2(di)=Theta2;
end
end

