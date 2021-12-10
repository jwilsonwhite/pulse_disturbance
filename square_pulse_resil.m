function [Noise, N, B, E, EP, Y, Params, Resist,RT]=square_pulse_resil(Species,FLEPs,Reduction,Duration,Delay,...
                                                            T,Ti,Conn_scenario, DD_scenario,...
                                                            RTvar,Thresh,Mfactor,AdultMortality,FReduct)

% N95, B95, Nmin, Bmin, mN, mB, miN, miB
% T is overall length of simulation
% Ti is pre-disturbance period

%Conn_scenario ='Closed';
%DD_scenario ='BH';

if ~exist('Reduction','var')
Reduction = 0.5;
end
if ~exist('Duration','var')
Duration = 5;
end
if ~exist('Delay','var')
Delay = NaN;
end
if ~exist('Mfactor','var')
Mfactor = 1;
end
if ~exist('AdultMortality','var')
AdultMortality=false;
end
if ~exist('FReduct','var')
FReduct=0;
end


Ti2 = Ti+(Duration-1); % subtract 1 so that Duration = 1 leads to a single year of disturbance, etc. in indexing
Ti3 = Ti2+Delay;
Ti4 = Ti3+(Duration-1);
% Ti5 = Ti4+Delay;
% Ti6 = Ti5+Duration;

%Disturbance timeseries
Noise = zeros(T,1);
Noise(Ti:Ti2) = Reduction;
if ~isnan(Delay)
    Noise(Ti3:Ti4) = Reduction;
end

% Set up the population dynamics
% (replace this later with call to a different function)
Params = define_Params(Species);
Params.M = Params.M*(1+(1-Mfactor)); % adjust M up or down

%keyboard

for f = 1:length(FLEPs)
F(f) = get_F(Params,FLEPs(f));
L = get_Leslie(Params,F(f),Conn_scenario,NaN,NaN); % to-be MPA patch
Leslie=L;

% Make Leslie an array so that it can be time-varying
LeslieArray = repmat(Leslie,[1,1,1,T]);
if FReduct ~= 0
    Leslie2 = get_Leslie(Params,F(f)*(1-FReduct),Conn_scenario,NaN,NaN);
    Fdelay = Ti+Params.Af-1; % so fishing is reduced the year a cohort enters the fishery
    LeslieArray(:,:,1,Fdelay:(Fdelay+Duration-1)) = repmat(Leslie2,[1,1,1,length(Fdelay:(Fdelay+Duration-1))]);
end

% Run the model, start at equilibrium
N0 = Params.SAD(:);
B0 = 0;
Y0 = 0; C0 = 0;
% quick loop to get to equilibrium
[N0,B0] = iterate_model(Params,LeslieArray,F(f),N0,Y0,C0,B0,T,Conn_scenario,DD_scenario, 0);

% Now add in pulse disturbance (could add in storing Y, C, etc. later; for
% now they are calculated but not stored)
Fvec = repmat(F(f),[1,T]);
if FReduct ~= 0
Fvec(Fdelay:(Fdelay+Duration-1)) = F(f)*(1-FReduct);
end
[Ntmp, Btmp, Y, E, EP] = iterate_model(Params,LeslieArray,Fvec,N0(:,end),Y0,C0,B0(:,end),T,Conn_scenario,DD_scenario, 0, Noise, AdultMortality);
N(:,:,:,f) = Ntmp;
B(:,:,:,f) = Btmp;




% Calculate magnitudes of initial drop and the return time
if ~exist('RTvar','var') % return-time variable
    RTvar = 'N';
end

if ~exist('Thresh','var') % threshold for return
    Thresh = 0.95;
end

switch RTvar
    case 'N'
        X = sum(squeeze(Ntmp));
    case 'B'
        X = sum(squeeze(Btmp));
end % end switch RTvar

X = X/X(Ti-1); % scale to initial

% Resistance
Resist(f) = min(X);

% Return time
if Resist >= Thresh % if it never dips far enough
    RT(f) = 0;
elseif X(end) < Thresh
    RT(f) = T-Ti;
else
    Xtmp = X(Ti:end);
    RT(f) = find(X(Ti:end) < Thresh,1,'last')+1; % do it from the end, so that we are not deceived by oscillations
end


end % end loop over FLEPs



