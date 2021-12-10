function [L,a] = get_Leslie(Params,F,Conn_scenario,Leslie_factor,Lambda_target)


options = optimoptions('fminunc','Display','none');
switch Conn_scenario
    
    case {'Closed','Porous'} % only add in reproduction if it is closed dynamics
    
if isnan(Leslie_factor) % if there is a target value of Lambda (for the unfished state)
    
    if ~isnan(Lambda_target)
    a0 = 1;
    Params_tmp = Params;
    Params_tmp.F = F;
   % a = fmincon(@(a)Leslie_min(a,Params_tmp,Lambda_target),a0,[],[],[],[],0,10,[],options);
   a = fminunc(@(a)Leslie_min(a,Params_tmp,Lambda_target),a0,options);
    else
    a = 1; % if Lambda_target comes in as NaN; just return the matrix without adjustment
    end
    
else % else if Lambda_target
a = Leslie_factor;
end % end if Lambda_target
    
% Create Leslie Matrix using the returned value of a:
Surv = ones(1,Params.A-1).*exp(-(Params.M + F.*Params.isFish(1:end-1)));
L = diag(Surv,-1);
L(1,:) = a*(Params.BiomassAge.^Params.exp).*Params.isMat;

    
    otherwise % if open pop scenario    
        if isnan(Leslie_factor) % if there is a target value of Lambda (for the unfished state)
    
    if ~isnan(Lambda_target)
    a0 = 1;
    Params_tmp = Params;
    Params_tmp.F = F;
   % a = fmincon(@(a)Leslie_min(a,Params_tmp,Lambda_target),a0,[],[],[],[],0,10,[],options);
   a = fminunc(@(a)Leslie_min(a,Params_tmp,Lambda_target),a0,options);
    else
    a = 1; % if Lambda_target comes in as NaN; just return the matrix without adjustment
    end
    
    else % else if Lambda_target
    a = Leslie_factor;
        end % end if Lambda_target

% Create Leslie Matrix
Surv = ones(1,Params.A-1).*exp(-(Params.M + F.*Params.isFish(1:end-1)));
L = diag(Surv,-1);
L(1,:) = a*(Params.BiomassAge.^Params.exp).*Params.isMat;

end % end switch conn_scenario

% The function to be minimized (an unfished Leslie matrix):
function X = Leslie_min(a,Params,Lambda_target)

% Create Leslie Matrix
%Surv = ones(1,Params.A-1).*exp(-(Params.M));
Surv = ones(1,Params.A-1).*exp(-(Params.M + Params.F.*Params.isFish(1:end-1)));
L = diag(Surv,-1);
L(1,:) = a*(Params.BiomassAge.^Params.exp).*Params.isMat;

X = abs( max(eig(L)) - Lambda_target);



