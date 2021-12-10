function FLEP=get_FLEP(Params,F);
Fs=Params.Fs;
Fmin=abs(Fs-F);
Index=find(Fmin == min(Fmin),1);
FLEPs = Params.FvLEP/Params.FvLEP(1);
FLEP=FLEPs(Index);
end