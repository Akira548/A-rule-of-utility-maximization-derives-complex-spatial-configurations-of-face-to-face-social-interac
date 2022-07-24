function [P_sys,P,dyadic] = CalculateSystemPotential(VH_pos)
%%% Formula
para = [1.523, 0.904, 4.173;
        0.042, 1.548, 0.247];
fun = @(beta,x) 1-exp(-( beta(3)./  ( (exp(beta(1).*abs(x(:,1)))).*x(:,2) )).^abs(beta(2)));
 
Theta1_Intensity = [];

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
     P(di) = geomean((((fun(para(1,:),[Theta2,Rho]).*(1-fun(para(2,:),[Theta2,Rho]))).*...
        (fun(para(1,:),[Theta1,Rho]).*(1-fun(para(2,:),[Theta1,Rho])))).^.5));   
end
    P_sys = geomean(P);
%     P_sys = prod(P);
end

