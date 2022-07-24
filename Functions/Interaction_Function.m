function [y] = Interaction_Function(Theta1,Theta2,Rho,mode)

%%%
if nargin < 4
    mode = 'Psocial';
end
Best_ori = [];

para = [-1.523, 0.904, 4.173;
        -0.042, 1.548, 0.247];
fun = @(beta,x) 1-exp(-( beta(3).*(exp(beta(1).*abs(x(:,1))))./x(:,2)).^beta(2));

switch mode
    case 'Patt'
        y = fun(para(1,:),[wrapToPi(Theta1),Rho]);
    case 'Prep'
        y = fun(para(2,:),[wrapToPi(Theta1),Rho]);    
    case 'Psocial'
        P1 = (fun(para(1,:),[wrapToPi(Theta2),Rho]).*(1-fun(para(2,:),[wrapToPi(Theta2),Rho])));
        P2 = (fun(para(1,:),[wrapToPi(Theta1),Rho]).*(1-fun(para(2,:),[wrapToPi(Theta1),Rho])));
        y = P1.*P2;
    case 'Pgroup'
        y = prod(fun(para(1,:),[wrapToPi(Theta2),Rho]).*(1-fun(para(2,:),[wrapToPi(Theta2),Rho])).*...
            (fun(para(1,:),[wrapToPi(Theta1),Rho]).*(1-fun(para(2,:),[wrapToPi(Theta1),Rho]))));
    case 'Pnon'
        y = mean(fun(para(1,:),[wrapToPi(Theta1),Rho]).*(1-fun(para(2,:),[wrapToPi(Theta1),Rho])));
    otherwise
        error('worng mode')
end
end

