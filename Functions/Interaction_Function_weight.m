function [y] = Interaction_Function_weight(Theta1,Theta2,Rho,w1,w2,mode)

%%%
if nargin < 6
    mode = 'Psocial';
end
Best_ori = [];

para = [-1.523, 0.904, 4.173;
        -0.042, 1.548, 0.247];
fun = @(beta,x) 1-exp(-( beta(3).*(exp(beta(1).*abs(x(:,1))))./x(:,2)).^beta(2));
switch mode
    case 'Patt'
        y = fun([para(1,:) w1],[wrapToPi(Theta1),Rho]);
    case 'Prep'
        y = fun([para(2,:) w2],[wrapToPi(Theta1),Rho]);    
    case 'Pgroup'
%         P1 = (fun([para(1,:) w1],[wrapToPi(Theta1),Rho]).*(1-fun([para(2,:) w2],[wrapToPi(Theta1),Rho])));
%         P2 = (fun([para(2,:) w3],[wrapToPi(Theta2),Rho]).*(1-fun([para(2,:) w4],[wrapToPi(Theta2),Rho])));
%         y = P1.*P2;
        for i = 1:size(w1,1)
            P1= fun(para(1,:),[wrapToPi(Theta1(i)),Rho(i),w1(i,2)]).*(1-fun(para(2,:),[wrapToPi(Theta1(i)),Rho(i),w2(i,2)]));
            P2= fun(para(1,:),[wrapToPi(Theta2(i)),Rho(i),w1(i,1)]).*(1-fun(para(2,:),[wrapToPi(Theta2(i)),Rho(i),w2(i,1)]));
            y(i) = P1.*P2;
        end
%     case 'Pgroup'
%         y = prod(fun([para(1,:) w1],[wrapToPi(Theta1),Rho]).*(1-fun([para(2,:) w2],[wrapToPi(Theta1),Rho])).*...
%             (fun([para(2,:) w3],[wrapToPi(Theta2),Rho]).*(1-fun([para(2,:) w4],[wrapToPi(Theta2),Rho]))));
    otherwise
        error('worng mode')
end
end

