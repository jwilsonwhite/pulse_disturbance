function F = get_F(Params,FLEP)

% Translate FLEP into F
FLEPs = Params.FvLEP./Params.FvLEP(1); % Divide by unfished LEP to get FLEP
FLEPmin = abs(FLEPs-FLEP); % absolute difference to target
Index = find(FLEPmin == min(FLEPmin),1); % find the index of the smallest difference
F = Params.Fs(Index); % that's the one!
end


