function Params = define_Params(Species,Mfactor)

% Function to store parameter values for MPA models
% Mfactor argument allows alternative values of M to be applied (only coded
% for blue rockfish as an example for Fig 2 in the text)

switch Species
    case 'Pacific cod'
        % Define the basic parameters
        Params = struct([]);
        Params(1).Linf = 99.46;
        Params.k = 0.188;
        Params.t0 = 0; % Using Female fishbase value
        Params.M = 0.499;
        Params.F = 0.2;
        Params.A = 20; 
        Params.Amat = 5; % 57.3cm
        Params.Af = 4; % 60 cm is 50% for main fleets
        Params.c = 5.631e-6; % biomass-length constant
        %Params.c = 1;
        %Params.d = 3.17; % biomass-length exponent species specific from
        % Kaplan et al. 2019
        Params.d = 3.1306; % biomass-length exponent (general)
        Params.exp = 1; % adjust this to 1.18; in order to perform hyperallometry adjustment from Barneche et al. 2018/Marshall et al. 2019

        % Recruitment rate in Open scenario
        Params.R = 1;
        
       % Params.steepness = 1; % DO NOT USE
        
       % Some derived quantities:
        Params.Ages = 1:Params.A;
        Params.isMat = Params.Ages>=Params.Amat; % logical - is this age mature?
        Params.isFish = Params.Ages>=Params.Af; % logical - is this age fished?
        % Length-at-age
        Params.Len = Params.Linf*(1 - exp(-(Params.k*(Params.Ages-Params.t0))));
        % Biomass-at-age
        Params.BiomassAge = Params.c.*Params.Len.^Params.d;
        % Stable age distribution
        Params.SAD = exp(-Params.M.*(Params.Ages-1));
        % Lifetime Egg Production (LEP, unfished):
        Params.LEP0 = ((Params.BiomassAge(:).^Params.exp).*Params.isMat(:))'*Params.SAD(:);
        Params.EP0 = (Params.BiomassAge(:).^Params.exp).*Params.isMat(:);
        % Generation time
        Params.Tg = ((Params.BiomassAge(:).^Params.exp).*Params.isMat(:).*(1:Params.A)')'*Params.SAD(:)/Params.LEP0;
        
        % F vs. LEP. We will do this in a quick vectorized way
        F = linspace(0,2,1e2); % 
        %F = linspace(0,0.5,1e2); % 
        Fmat = repmat(F(:)',[Params.A,1]); % matrix of all possible Fs (on cols) for each age (on rows)
        Fmat(Params.Ages<Params.Af,:) = 0; % no fishing below age Af
        Fmatbin=Fmat(Params.Ages<Params.Af,:) == 0;% no fishing below age Af
        Agemat = repmat(Params.Ages(:),[1,length(F)]); % matrix of ages
        SADmat = exp(-(Params.M+Fmat).*(Agemat-1)); % stable age distribution given each F
        LEPmat = sum(SADmat.*repmat((Params.BiomassAge(:).^Params.exp).*Params.isMat(:),[1,length(F)])); % SAD * eggs-at-age = LEP 
        Params.Fs = F;
        Params.FvLEP = LEPmat./Params.LEP0; % sum of SAD
        
    case 'Kelp rockfish'
        % Define the basic parameters
        Params = struct([]);
        Params(1).Linf = 37.80;
        Params.k = 0.23;
        Params.t0 = 0; % Using Female fishbase value
        Params.M = 0.20;
        Params.F = 0.17;
        Params.A = 25;
        Params.Amat = 3;
        Params.Af = 3;
        Params.c = 1; % biomass-length constant
        %Params.d = 3.17; % biomass-length exponent species specific from
        % Kaplan et al. 2019
        Params.d = 3; % biomass-length exponent (general)
        Params.exp = 1; % adjust this to 1.18; in order to perform hyperallometry adjustment from Barneche et al. 2018/Marshall et al. 2019

        % Recruitment rate in Open scenario
        Params.R = 1;
        
        % Some derived quantities:
        Params.Ages = 1:Params.A;
        Params.isMat = Params.Ages>=Params.Amat; % logical - is this age mature?
        Params.isFish = Params.Ages>=Params.Af; % logical - is this age fished?
        % Length-at-age
        Params.Len = Params.Linf*(1 - exp(-(Params.k*(Params.Ages-Params.t0))));
        % Biomass-at-age
        Params.BiomassAge = Params.c.*Params.Len.^Params.d;
        % Stable age distribution
        Params.SAD = exp(-Params.M.*(Params.Ages-1));
        % Lifetime Egg Production (LEP, unfished):
        Params.LEP0 = ((Params.BiomassAge(:).^Params.exp).*Params.isMat(:))'*Params.SAD(:);
        Params.EP0 = (Params.BiomassAge(:).^Params.exp).*Params.isMat(:);
        % Generation time
        Params.Tg = ((Params.BiomassAge(:).^Params.exp).*Params.isMat(:).*(1:Params.A)')'*Params.SAD(:)/Params.LEP0;
        
        
        % F vs. LEP. We will do this in a quick vectorized way
        F = linspace(0,2,1e2); % 
        %F = linspace(0,0.5,1e2); % 
        Fmat = repmat(F(:)',[Params.A,1]); % matrix of all possible Fs (on cols) for each age (on rows)
        Fmat(Params.Ages<Params.Af,:) = 0;
        Params.Fmatbin=(Fmat(Params.Ages<Params.Af,:) == 0);% no fishing below age Af
        Agemat = repmat(Params.Ages(:),[1,length(F)]); % matrix of ages
        SADmat = exp(-(Params.M+Fmat).*(Agemat-1)); % stable age distribution given each F
        LEPmat = sum(SADmat.*repmat((Params.BiomassAge(:).^Params.exp).*Params.isMat(:),[1,length(F)])); % SAD * eggs-at-age = LEP 
        Params.Fs = F;
        Params.FvLEP = LEPmat./Params.LEP0; % sum of SAD
    case 'Blue rockfish'
        % Define the basic parameters
        Params = struct([]);
        Params(1).Linf = 38.15;
        Params.k = 0.17; % 0.17
        Params.t0 = -1.34; %-1.34; % Using Female fishbase value
        
        Params.M = 0.14; % 0.14
        if exist('Mfactor','var')
            Params.M = 0.14 * Mfactor;
        end
        Params.A = 44; % should be 44;
        Params.Amat = 6; % should be 6
        Params.Amat_sd = 1; % bc youngest age at maturity is 4
        Params.Af = 8;
        Params.Af_sd = 2.5; % based on length-based selectivities from Dick et al
        Params.c = 1; % biomass-length constant
        Params.d = 3.09; % biomass-length exponent species specific from
        % Kaplan et al. 2019
        %Params.d = 3; % biomass-length exponent (general)
        Params.exp = 1; %1.18; % adjust this to 1.18; in order to perform hyperallometry adjustment from Barneche et al. 2018/Marshall et al. 2019
                
        % Recruitment rate in Open scenario
        Params.R = 1;
        
        % Some derived quantities:
        Params.Ages = 1:Params.A;
        Params.isMat = normcdf(Params.Ages,Params.Amat,Params.Amat_sd); % logical - is this age mature?
        Params.isFish = normcdf(Params.Ages,Params.Af,Params.Af_sd); % logical - is this age fished?
        % Length-at-age
        Params.Len = Params.Linf*(1 - exp(-(Params.k*(Params.Ages-Params.t0))));
        % Biomass-at-age
        Params.BiomassAge = Params.c.*Params.Len.^Params.d;
        % Stable age distribution
        Params.SAD = exp(-Params.M.*(Params.Ages-1));
        % Lifetime Egg Production (LEP, unfished):
        Params.LEP0 = ((Params.BiomassAge(:).^Params.exp).*Params.isMat(:))'*Params.SAD(:);
        Params.EP0 = (Params.BiomassAge(:).^Params.exp).*Params.isMat(:);
        % Generation time
        Params.Tg = ((Params.BiomassAge(:).^Params.exp).*Params.isMat(:).*(1:Params.A)')'*Params.SAD(:)/Params.LEP0;
        
        % F vs. LEP. We will do this in a quick vectorized way
        F = linspace(0,2,1e2); % 
        %F = linspace(0,0.5,1e2); % 
        Fmat = repmat(F(:)',[Params.A,1]); % matrix of all possible Fs (on cols) for each age (on rows)
        Fmat(Params.Ages<Params.Af,:) = 0; % no fishing below age Af
        %Fmatbin=Fmat(Params.Ages<Params.Af,:) == 0;% no fishing below age Af
        Agemat = repmat(Params.Ages(:),[1,length(F)]); % matrix of ages
        SADmat = exp(-(Params.M+Fmat).*(Agemat-1)); % stable age distribution given each F
        LEPmat = sum(SADmat.*repmat((Params.BiomassAge(:).^Params.exp).*Params.isMat(:),[1,length(F)])); % SAD * eggs-at-age = LEP 
        Params.Fs = F;
        Params.FvLEP = LEPmat./Params.LEP0; % sum of SAD
    case 'Black rockfish'
        % Define the basic parameters
        Params = struct([]);
        Params(1).Linf = 45.11;
        Params.k = 0.33;
        Params.t0 = 0; % Using Female fishbase value
        Params.M = 0.18;
        if exist('Mfactor','var')
            Params.M = 0.18 * Mfactor;
        end
        Params.A = 50;
        Params.Amat = 7;
        Params.Af = 4;
        Params.c = 1; % biomass-length constant
        %Params.d = 3.19; % biomass-length exponent species specific from
        % Kaplan et al. 2019
        Params.d = 3; % biomass-length exponent (general)
        %Params.d = 3.19; % biomass-length exponent species specific 
        Params.exp = 1; % adjust this to 1.18; in order to perform hyperallometry adjustment from Barneche et al. 2018/Marshall et al. 2019
        % above and adjustment given hyperallometry adjustment from Barneche et al. 2018/Marshall et al. 2019
        
        % Recruitment rate in Open scenario
        Params.R = 1;
        
        % Some derived quantities:
        Params.Ages = 1:Params.A;
        Params.isMat = Params.Ages>=Params.Amat; % logical - is this age mature?
        Params.isFish = Params.Ages>=Params.Af; % logical - is this age fished?
        % Length-at-age
        Params.Len = Params.Linf*(1 - exp(-(Params.k*(Params.Ages-Params.t0))));
        % Biomass-at-age
        Params.BiomassAge = Params.c.*Params.Len.^Params.d;
        % Stable age distribution
        Params.SAD = exp(-Params.M.*(Params.Ages-1));
        % Lifetime Egg Production (LEP, unfished):
        Params.LEP0 = ((Params.BiomassAge(:).^Params.exp).*Params.isMat(:))'*Params.SAD(:);
        Params.EP0 = (Params.BiomassAge(:).^Params.exp).*Params.isMat(:);
        % Generation time
        Params.Tg = ((Params.BiomassAge(:).^Params.exp).*Params.isMat(:).*(1:Params.A)')'*Params.SAD(:)/Params.LEP0;
        
        % F vs. LEP. We will do this in a quick vectorized way
        F = linspace(0,2,1e2); % 
        %F = linspace(0,0.5,1e2); % 
        Fmat = repmat(F(:)',[Params.A,1]); % matrix of all possible Fs (on cols) for each age (on rows)
        Fmat(Params.Ages<Params.Af,:) = 0; % no fishing below age Af
        Fmatbin=Fmat(Params.Ages<Params.Af,:) == 0;% no fishing below age Af
        Agemat = repmat(Params.Ages(:),[1,length(F)]); % matrix of ages
        SADmat = exp(-(Params.M+Fmat).*(Agemat-1)); % stable age distribution given each F
        LEPmat = sum(SADmat.*repmat((Params.BiomassAge(:).^Params.exp).*Params.isMat(:),[1,length(F)])); % SAD * eggs-at-age = LEP 
        Params.Fs = F;
        Params.FvLEP = LEPmat./Params.LEP0; % sum of SAD
    case 'Gopher rockfish'
        % Define the basic parameters
        Params = struct([]);
        Params(1).Linf = 34.10;
        Params.k = 0.23;
        Params.t0 = 0; % Using Female fishbase value
        Params.M = .2;
        Params.A = 30;
        Params.Amat = 3;
        Params.Af = 6;
        Params.c = 1; % biomass-length constant
        %Params.d = 3.08; % biomass-length exponent species specific from
        % Kaplan et al. 2019
        Params.d = 3; % biomass-length exponent (general)
        %Params.d = 3.08;
         Params.exp = 1; % adjust this to 1.18; in order to perform hyperallometry adjustment from Barneche et al. 2018/Marshall et al. 2019

        % Recruitment rate in Open scenario
        Params.R = 1;
        
        % Some derived quantities:
        Params.Ages = 1:Params.A;
        Params.isMat = Params.Ages>=Params.Amat; % logical - is this age mature?
        Params.isFish = Params.Ages>=Params.Af; % logical - is this age fished?
        % Length-at-age
        Params.Len = Params.Linf*(1 - exp(-(Params.k*(Params.Ages-Params.t0))));
        % Biomass-at-age
        Params.BiomassAge = Params.c.*Params.Len.^Params.d;
        % Stable age distribution
        Params.SAD = exp(-Params.M.*(Params.Ages-1));
        % Lifetime Egg Production (LEP, unfished):
        Params.LEP0 = ((Params.BiomassAge(:).^Params.exp).*Params.isMat(:))'*Params.SAD(:);
        Params.EP0 = (Params.BiomassAge(:).^Params.exp).*Params.isMat(:);
        % Generation time
        Params.Tg = ((Params.BiomassAge(:).^Params.exp).*Params.isMat(:).*(1:Params.A)')'*Params.SAD(:)/Params.LEP0;
        
        % F vs. LEP. We will do this in a quick vectorized way
        F = linspace(0,2,1e2); % 
        %F = linspace(0,0.5,1e2); % 
        Fmat = repmat(F(:)',[Params.A,1]); % matrix of all possible Fs (on cols) for each age (on rows)
        Fmat(Params.Ages<Params.Af,:) = 0; % no fishing below age Af
        Fmatbin=Fmat(Params.Ages<Params.Af,:) == 0;% no fishing below age Af
        Agemat = repmat(Params.Ages(:),[1,length(F)]); % matrix of ages
        SADmat = exp(-(Params.M+Fmat).*(Agemat-1)); % stable age distribution given each F
        LEPmat = sum(SADmat.*repmat((Params.BiomassAge(:).^Params.exp).*Params.isMat(:),[1,length(F)])); % SAD * eggs-at-age = LEP 
        Params.Fs = F;
        Params.FvLEP = LEPmat./Params.LEP0; % sum of SAD   
    case 'Lingcod'
        % Define the basic parameters
        Params = struct([]);
        Params(1).Linf = 96.74;
        Params.k = 0.17;
        Params.t0 = 0; % Using Female fishbase value
        Params.M = .25;
        Params.A = 25;
        Params.Amat = 3;
        Params.Af = 4;
        Params.c = 1; % biomass-length constant
        %Params.d = 3.41; % biomass-length exponent species specific from
        % Kaplan et al. 2019
        Params.d = 3; % biomass-length exponent (general)
        %Params.d = 3.41;
                Params.exp = 1; % adjust this to 1.18; in order to perform hyperallometry adjustment from Barneche et al. 2018/Marshall et al. 2019

                
        % Recruitment rate in Open scenario
        Params.R = 1;
        
        % Some derived quantities:
        Params.Ages = 1:Params.A;
        Params.isMat = Params.Ages>=Params.Amat; % logical - is this age mature?
        Params.isFish = Params.Ages>=Params.Af; % logical - is this age fished?
        % Length-at-age
        Params.Len = Params.Linf*(1 - exp(-(Params.k*(Params.Ages-Params.t0))));
        % Biomass-at-age
        Params.BiomassAge = Params.c.*Params.Len.^Params.d;
        % Stable age distribution
        Params.SAD = exp(-Params.M.*(Params.Ages-1));
        % Lifetime Egg Production (LEP, unfished):
        Params.LEP0 = ((Params.BiomassAge(:).^Params.exp).*Params.isMat(:))'*Params.SAD(:);
        Params.EP0 = (Params.BiomassAge(:).^Params.exp).*Params.isMat(:);
        % Generation time
        Params.Tg = ((Params.BiomassAge(:).^Params.exp).*Params.isMat(:).*(1:Params.A)')'*Params.SAD(:)/Params.LEP0;
        
        % F vs. LEP. We will do this in a quick vectorized way
        F = linspace(0,2,1e2); % 
        %F = linspace(0,0.5,1e2); % 
        Fmat = repmat(F(:)',[Params.A,1]); % matrix of all possible Fs (on cols) for each age (on rows)
        Fmat(Params.Ages<Params.Af,:) = 0; % no fishing below age Af
        Fmatbin=Fmat(Params.Ages<Params.Af,:) == 0;% no fishing below age Af
        Agemat = repmat(Params.Ages(:),[1,length(F)]); % matrix of ages
        SADmat = exp(-(Params.M+Fmat).*(Agemat-1)); % stable age distribution given each F
        LEPmat = sum(SADmat.*repmat((Params.BiomassAge(:).^Params.exp).*Params.isMat(:),[1,length(F)])); % SAD * eggs-at-age = LEP 
        Params.Fs = F;
        Params.FvLEP = LEPmat./Params.LEP0; % sum of SAD
    case 'Copper rockfish'
        % Define the basic parameters
        Params = struct([]);
        Params(1).Linf = 56.50;
        Params.k = 0.14;
        Params.t0 = 0; % Using Female fishbase value
        Params.M = .09;
        Params.A = 100;
        Params.Amat = 5;
        Params.Af = 5;
        Params.c = 1; % biomass-length constant
        %Params.d = 3.13; % biomass-length exponent species specific from
        % Kaplan et al. 2019
        Params.d = 3; % biomass-length exponent (general)
        %Params.d = 3.13;
                Params.exp = 1; % adjust this to 1.18; in order to perform hyperallometry adjustment from Barneche et al. 2018/Marshall et al. 2019

        
        % Recruitment rate in Open scenario
        Params.R = 1;
        
        % Some derived quantities:
        Params.Ages = 1:Params.A;
        Params.isMat = Params.Ages>=Params.Amat; % logical - is this age mature?
        Params.isFish = Params.Ages>=Params.Af; % logical - is this age fished?
        % Length-at-age
        Params.Len = Params.Linf*(1 - exp(-(Params.k*(Params.Ages-Params.t0))));
        % Biomass-at-age
        Params.BiomassAge = Params.c.*Params.Len.^Params.d;
        % Stable age distribution
        Params.SAD = exp(-Params.M.*(Params.Ages-1));
        % Lifetime Egg Production (LEP, unfished):
        Params.LEP0 = ((Params.BiomassAge(:).^Params.exp).*Params.isMat(:))'*Params.SAD(:);
        Params.EP0 = (Params.BiomassAge(:).^Params.exp).*Params.isMat(:);
        % Generation time
        Params.Tg = ((Params.BiomassAge(:).^Params.exp).*Params.isMat(:).*(1:Params.A)')'*Params.SAD(:)/Params.LEP0;
        
        % F vs. LEP. We will do this in a quick vectorized way
        F = linspace(0,2,1e2); % 
        %F = linspace(0,0.5,1e2); % 
        Fmat = repmat(F(:)',[Params.A,1]); % matrix of all possible Fs (on cols) for each age (on rows)
        Fmat(Params.Ages<Params.Af,:) = 0; % no fishing below age Af
        Fmatbin=Fmat(Params.Ages<Params.Af,:) == 0;% no fishing below age Af
        Agemat = repmat(Params.Ages(:),[1,length(F)]); % matrix of ages
        SADmat = exp(-(Params.M+Fmat).*(Agemat-1)); % stable age distribution given each F
        LEPmat = sum(SADmat.*repmat((Params.BiomassAge(:).^Params.exp).*Params.isMat(:),[1,length(F)])); % SAD * eggs-at-age = LEP 
        Params.Fs = F;
        Params.FvLEP = LEPmat./Params.LEP0; % sum of SAD
    case 'California scorpionfish'
        % Define the basic parameters
        Params = struct([]);
        Params(1).Linf = 40.29;
        Params.k = 0.13;
        Params.t0 = -1.90; % Using Female fishbase value
        Params.M = .25;
        Params.A = 21;
        Params.Amat = 2;
        Params.Af = 5;
        Params.c = 1; % biomass-length constant
        %Params.d = 3.17; % biomass-length exponent species specific from
        % Kaplan et al. 2019
        Params.d = 3; % biomass-length exponent (general)
        %Params.d = 3.17;
                Params.exp = 1; % adjust this to 1.18; in order to perform hyperallometry adjustment from Barneche et al. 2018/Marshall et al. 2019

        % Recruitment rate in Open scenario
        Params.R = 1;
        
        % Some derived quantities:
        Params.Ages = 1:Params.A;
        Params.isMat = Params.Ages>=Params.Amat; % logical - is this age mature?
        Params.isFish = Params.Ages>=Params.Af; % logical - is this age fished?
        % Length-at-age
        Params.Len = Params.Linf*(1 - exp(-(Params.k*(Params.Ages-Params.t0))));
        % Biomass-at-age
        Params.BiomassAge = Params.c.*Params.Len.^Params.d;
        % Stable age distribution
        Params.SAD = exp(-Params.M.*(Params.Ages-1));
        % Lifetime Egg Production (LEP, unfished):
        Params.LEP0 = ((Params.BiomassAge(:).^Params.exp).*Params.isMat(:))'*Params.SAD(:);
        Params.EP0 = (Params.BiomassAge(:).^Params.exp).*Params.isMat(:);
        % Generation time
        Params.Tg = ((Params.BiomassAge(:).^Params.exp).*Params.isMat(:).*(1:Params.A)')'*Params.SAD(:)/Params.LEP0;
        
        % F vs. LEP. We will do this in a quick vectorized way
        F = linspace(0,2,1e2); % 
        %F = linspace(0,0.5,1e2); % 
        Fmat = repmat(F(:)',[Params.A,1]); % matrix of all possible Fs (on cols) for each age (on rows)
        Fmat(Params.Ages<Params.Af,:) = 0; % no fishing below age Af
        Fmatbin=Fmat(Params.Ages<Params.Af,:) == 0;% no fishing below age Af
        Agemat = repmat(Params.Ages(:),[1,length(F)]); % matrix of ages
        SADmat = exp(-(Params.M+Fmat).*(Agemat-1)); % stable age distribution given each F
        LEPmat = sum(SADmat.*repmat((Params.BiomassAge(:).^Params.exp).*Params.isMat(:),[1,length(F)])); % SAD * eggs-at-age = LEP 
        Params.Fs = F;
        Params.FvLEP = LEPmat./Params.LEP0; % sum of SAD
    case 'Brown rockfish'
        % Define the basic parameters
        Params = struct([]);
        Params(1).Linf = 51.40;
        Params.k = 0.16;
        Params.t0 = 0; % Using Female fishbase value
        Params.M = .14;
        Params.A = 34;
        Params.Amat = 4;
        Params.Af = 4;
        Params.c = 1; % biomass-length constant
        %Params.d = 3.07; % biomass-length exponent species specific from
        % Kaplan et al. 2019
        Params.d = 3; % biomass-length exponent (general)
        %Params.d = 3.07;
                Params.exp = 1; % adjust this to 1.18; in order to perform hyperallometry adjustment from Barneche et al. 2018/Marshall et al. 2019

        % Recruitment rate in Open scenario
        Params.R = 1;
        
        % Some derived quantities:
        Params.Ages = 1:Params.A;
        Params.isMat = Params.Ages>=Params.Amat; % logical - is this age mature?
        Params.isFish = Params.Ages>=Params.Af; % logical - is this age fished?
        % Length-at-age
        Params.Len = Params.Linf*(1 - exp(-(Params.k*(Params.Ages-Params.t0))));
        % Biomass-at-age
        Params.BiomassAge = Params.c.*Params.Len.^Params.d;
        % Stable age distribution
        Params.SAD = exp(-Params.M.*(Params.Ages-1));
        % Lifetime Egg Production (LEP, unfished):
        Params.LEP0 = ((Params.BiomassAge(:).^Params.exp).*Params.isMat(:))'*Params.SAD(:);
        Params.EP0 = (Params.BiomassAge(:).^Params.exp).*Params.isMat(:);
        % Generation time
        Params.Tg = ((Params.BiomassAge(:).^Params.exp).*Params.isMat(:).*(1:Params.A)')'*Params.SAD(:)/Params.LEP0;
        
        % F vs. LEP. We will do this in a quick vectorized way
        F = linspace(0,2,1e2); % 
        %F = linspace(0,0.5,1e2); % 
        Fmat = repmat(F(:)',[Params.A,1]); % matrix of all possible Fs (on cols) for each age (on rows)
        Fmat(Params.Ages<Params.Af,:) = 0; % no fishing below age Af
        Fmatbin=Fmat(Params.Ages<Params.Af,:) == 0;% no fishing below age Af
        Agemat = repmat(Params.Ages(:),[1,length(F)]); % matrix of ages
        SADmat = exp(-(Params.M+Fmat).*(Agemat-1)); % stable age distribution given each F
        LEPmat = sum(SADmat.*repmat((Params.BiomassAge(:).^Params.exp).*Params.isMat(:),[1,length(F)])); % SAD * eggs-at-age = LEP 
        Params.Fs = F;
        Params.FvLEP = LEPmat./Params.LEP0; % sum of SAD
    case 'Yellowtail rockfish'
        % Define the basic parameters
        Params = struct([]);
        Params(1).Linf = 49.88;
        Params.k = 0.18;
        Params.t0 = -1.95; % Using Female fishbase value
        Params.M = .11;
        Params.A = 64;
        Params.Amat = 6;
        Params.Af = 4;
        Params.c = 1; % biomass-length constant
        %Params.d = 3; % biomass-length exponent species specific from
        % Kaplan et al. 2019
        Params.d = 3; % biomass-length exponent (general)
        %Params.d = 3;
                Params.exp = 1; % adjust this to 1.18; in order to perform hyperallometry adjustment from Barneche et al. 2018/Marshall et al. 2019

        % Recruitment rate in Open scenario
        Params.R = 1;
        
        % Some derived quantities:
        Params.Ages = 1:Params.A;
        Params.isMat = Params.Ages>=Params.Amat; % logical - is this age mature?
        Params.isFish = Params.Ages>=Params.Af; % logical - is this age fished?
        % Length-at-age
        Params.Len = Params.Linf*(1 - exp(-(Params.k*(Params.Ages-Params.t0))));
        % Biomass-at-age
        Params.BiomassAge = Params.c.*Params.Len.^Params.d;
        % Stable age distribution
        Params.SAD = exp(-Params.M.*(Params.Ages-1));
        % Lifetime Egg Production (LEP, unfished):
        Params.LEP0 = ((Params.BiomassAge(:).^Params.exp).*Params.isMat(:))'*Params.SAD(:);
        Params.EP0 = (Params.BiomassAge(:).^Params.exp).*Params.isMat(:);
        % Generation time
        Params.Tg = ((Params.BiomassAge(:).^Params.exp).*Params.isMat(:).*(1:Params.A)')'*Params.SAD(:)/Params.LEP0;
        
        % F vs. LEP. We will do this in a quick vectorized way
        F = linspace(0,2,1e2); % 
        %F = linspace(0,0.5,1e2); % 
        Fmat = repmat(F(:)',[Params.A,1]); % matrix of all possible Fs (on cols) for each age (on rows)
        Fmat(Params.Ages<Params.Af,:) = 0; % no fishing below age Af
        Agemat = repmat(Params.Ages(:),[1,length(F)]); % matrix of ages
        SADmat = exp(-(Params.M+Fmat).*(Agemat-1)); % stable age distribution given each F
        LEPmat = sum(SADmat.*repmat((Params.BiomassAge(:).^Params.exp).*Params.isMat(:),[1,length(F)])); % SAD * eggs-at-age = LEP 
        Params.Fs = F;
        Params.FvLEP = LEPmat./Params.LEP0; % sum of SAD
    case 'Vermilion rockfish'
        % Define the basic parameters
        Params = struct([]);
        Params(1).Linf = 53.92;
        Params.k = 0.16;
        Params.t0 = 0; % Using Female fishbase value
        Params.M = 0.10;
        Params.A = 65;
        Params.Amat = 7;
        Params.Af = 3;
        Params.c = 1; % biomass-length constant
        %Params.d = 3.04; % biomass-length exponent species specific from
        % Kaplan et al. 2019
        Params.d = 3; % biomass-length exponent (general)
        %Params.d = 3.04;
                Params.exp = 1; % adjust this to 1.18; in order to perform hyperallometry adjustment from Barneche et al. 2018/Marshall et al. 2019

        % Recruitment rate in Open scenario
        Params.R = 1;
        
        % Some derived quantities:
        Params.Ages = 1:Params.A;
        Params.isMat = Params.Ages>=Params.Amat; % logical - is this age mature?
        Params.isFish = Params.Ages>=Params.Af; % logical - is this age fished?
        % Length-at-age
        Params.Len = Params.Linf*(1 - exp(-(Params.k*(Params.Ages-Params.t0))));
        % Biomass-at-age
        Params.BiomassAge = Params.c.*Params.Len.^Params.d;
        % Stable age distribution
        Params.SAD = exp(-Params.M.*(Params.Ages-1));
        % Lifetime Egg Production (LEP, unfished):
        Params.LEP0 = ((Params.BiomassAge(:).^Params.exp).*Params.isMat(:))'*Params.SAD(:);
        Params.EP0 = (Params.BiomassAge(:).^Params.exp).*Params.isMat(:);
        % Generation time
        Params.Tg = ((Params.BiomassAge(:).^Params.exp).*Params.isMat(:).*(1:Params.A)')'*Params.SAD(:)/Params.LEP0;
        
        % F vs. LEP. We will do this in a quick vectorized way
        F = linspace(0,2,1e2); % 
        %F = linspace(0,0.5,1e2); % 
        Fmat = repmat(F(:)',[Params.A,1]); % matrix of all possible Fs (on cols) for each age (on rows)
        Fmat(Params.Ages<Params.Af,:) = 0; % no fishing below age Af
        Agemat = repmat(Params.Ages(:),[1,length(F)]); % matrix of ages
        SADmat = exp(-(Params.M+Fmat).*(Agemat-1)); % stable age distribution given each F
        LEPmat = sum(SADmat.*repmat((Params.BiomassAge(:).^Params.exp).*Params.isMat(:),[1,length(F)])); % SAD * eggs-at-age = LEP 
        Params.Fs = F;
        Params.FvLEP = LEPmat./Params.LEP0; % sum of SAD
  case 'Bocaccio'
        % Define the basic parameters
        Params = struct([]);
        Params(1).Linf = 70 ;
        Params.k = 0.22;
        Params.t0 = 0; % Using Female fishbase value
        Params.M = 0.15;
        Params.A = 55;
        Params.Amat = 3;
        Params.Af = 3;
        Params.c = 1; % biomass-length constant
        %Params.d = 3.19; % biomass-length exponent species specific from
        % Kaplan et al. 2019
        Params.d = 3; % biomass-length exponent (general)
        %Params.d = 3.19;
                Params.exp = 1; % adjust this to 1.18; in order to perform hyperallometry adjustment from Barneche et al. 2018/Marshall et al. 2019

        % Recruitment rate in Open scenario
        Params.R = 1;
        
        % Some derived quantities:
        Params.Ages = 1:Params.A;
        Params.isMat = Params.Ages>=Params.Amat; % logical - is this age mature?
        Params.isFish = Params.Ages>=Params.Af; % logical - is this age fished?
        % Length-at-age
        Params.Len = Params.Linf*(1 - exp(-(Params.k*(Params.Ages-Params.t0))));
        % Biomass-at-age
        Params.BiomassAge = Params.c.*Params.Len.^Params.d;
        % Stable age distribution
        Params.SAD = exp(-Params.M.*(Params.Ages-1));
        % Lifetime Egg Production (LEP, unfished):
        Params.LEP0 = ((Params.BiomassAge(:).^Params.exp).*Params.isMat(:))'*Params.SAD(:);
        Params.EP0 = (Params.BiomassAge(:).^Params.exp).*Params.isMat(:);
        % Generation time
        Params.Tg = ((Params.BiomassAge(:).^Params.exp).*Params.isMat(:).*(1:Params.A)')'*Params.SAD(:)/Params.LEP0;
        
        % F vs. LEP. We will do this in a quick vectorized way
        F = linspace(0,2,1e2); % 
        %F = linspace(0,0.5,1e2); % 
        Fmat = repmat(F(:)',[Params.A,1]); % matrix of all possible Fs (on cols) for each age (on rows)
        Fmat(Params.Ages<Params.Af,:) = 0; % no fishing below age Af
        Agemat = repmat(Params.Ages(:),[1,length(F)]); % matrix of ages
        SADmat = exp(-(Params.M+Fmat).*(Agemat-1)); % stable age distribution given each F
        LEPmat = sum(SADmat.*repmat((Params.BiomassAge(:).^Params.exp).*Params.isMat(:),[1,length(F)])); % SAD * eggs-at-age = LEP 
        Params.Fs = F;
        Params.FvLEP = LEPmat./Params.LEP0; % sum of SAD
    case 'China rockfish'
        % Define the basic parameters
        Params = struct([]);
        Params(1).Linf = 33.62;
        Params.k = 0.23;
        Params.t0 = 0; % CHECK ON THIS VALUE IN FISHBASE
        Params.M = 0.06;
        Params.A = 83;
        Params.Amat = 7;
        Params.Af = 10;
        Params.c = 1; % biomass-length constant
        %Params.d = 3.18; % biomass-length exponent species specific from
        % Kaplan et al. 2019
        Params.d = 3; % biomass-length exponent (general)
                Params.exp = 1; % adjust this to 1.18; in order to perform hyperallometry adjustment from Barneche et al. 2018/Marshall et al. 2019

        % Recruitment rate in Open scenario
        Params.R = 1;
        
        % Some derived quantities:
        Params.Ages = 1:Params.A;
        Params.isMat = Params.Ages>=Params.Amat; % logical - is this age mature?
        Params.isFish = Params.Ages>=Params.Af; % logical - is this age fished?
        % Length-at-age
        Params.Len = Params.Linf*(1 - exp(-(Params.k*(Params.Ages-Params.t0))));
        % Biomass-at-age
        Params.BiomassAge = Params.c.*Params.Len.^Params.d;
        % Stable age distribution
        Params.SAD = exp(-Params.M.*(Params.Ages-1));
        % Lifetime Egg Production (LEP, unfished):
        Params.LEP0 = ((Params.BiomassAge(:).^Params.exp).*Params.isMat(:))'*Params.SAD(:);
        Params.EP0 = (Params.BiomassAge(:).^Params.exp).*Params.isMat(:);
        % Generation time
        Params.Tg = ((Params.BiomassAge(:).^Params.exp).*Params.isMat(:).*(1:Params.A)')'*Params.SAD(:)/Params.LEP0;
        
        % F vs. LEP. We will do this in a quick vectorized way
        F = linspace(0,2,1e2); % 
        %F = linspace(0,0.5,1e2); % 
        Fmat = repmat(F(:)',[Params.A,1]); % matrix of all possible Fs (on cols) for each age (on rows)
        Fmat(Params.Ages<Params.Af,:) = 0; % no fishing below age Af
        Agemat = repmat(Params.Ages(:),[1,length(F)]); % matrix of ages
        SADmat = exp(-(Params.M+Fmat).*(Agemat-1)); % stable age distribution given each F
        LEPmat = sum(SADmat.*repmat((Params.BiomassAge(:).^Params.exp).*Params.isMat(:),[1,length(F)])); % SAD * eggs-at-age = LEP 
        Params.Fs = F;
        Params.FvLEP = LEPmat./Params.LEP0; % sum of SAD
     case 'Cabezon'
        % Define the basic parameters
        Params = struct([]);
        Params(1).Linf = 49.90;
        Params.k = 0.28;
        Params.t0 = 0; % Using Female fishbase value
        Params.M = 0.28;
        Params.A = 17;
        Params.Amat = 3;
        Params.Af = 4;
        Params.c = 1; % biomass-length constant
        %Params.d = 3.19; % biomass-length exponent species specific from
        % Kaplan et al. 2019
        Params.d = 3; % biomass-length exponent (general)
        %Params.d = 3.19;
                Params.exp = 1; % adjust this to 1.18; in order to perform hyperallometry adjustment from Barneche et al. 2018/Marshall et al. 2019

        % Recruitment rate in Open scenario
        Params.R = 1;
        
        % Some derived quantities:
        Params.Ages = 1:Params.A;
        Params.isMat = Params.Ages>=Params.Amat; % logical - is this age mature?
        Params.isFish = Params.Ages>=Params.Af; % logical - is this age fished?
        % Length-at-age
        Params.Len = Params.Linf*(1 - exp(-(Params.k*(Params.Ages-Params.t0))));
        % Biomass-at-age
        Params.BiomassAge = Params.c.*Params.Len.^Params.d;
        % Stable age distribution
        Params.SAD = exp(-Params.M.*(Params.Ages-1));
        % Lifetime Egg Production (LEP, unfished):
        Params.LEP0 = ((Params.BiomassAge(:).^Params.exp).*Params.isMat(:))'*Params.SAD(:);
        Params.EP0 = (Params.BiomassAge(:).^Params.exp).*Params.isMat(:);
        % Generation time
        Params.Tg = ((Params.BiomassAge(:).^Params.exp).*Params.isMat(:).*(1:Params.A)')'*Params.SAD(:)/Params.LEP0;
        
        % F vs. LEP. We will do this in a quick vectorized way
        F = linspace(0,2,1e2); % 
        %F = linspace(0,0.5,1e2); % 
        Fmat = repmat(F(:)',[Params.A,1]); % matrix of all possible Fs (on cols) for each age (on rows)
        Fmat(Params.Ages<Params.Af,:) = 0; % no fishing below age Af
        Agemat = repmat(Params.Ages(:),[1,length(F)]); % matrix of ages
        SADmat = exp(-(Params.M+Fmat).*(Agemat-1)); % stable age distribution given each F
        LEPmat = sum(SADmat.*repmat((Params.BiomassAge(:).^Params.exp).*Params.isMat(:),[1,length(F)])); % SAD * eggs-at-age = LEP 
        Params.Fs = F;
        Params.FvLEP = LEPmat./Params.LEP0; % sum of SAD
    case 'Kelp greenling'
        % Define the basic parameters
        Params = struct([]);
        Params(1).Linf = 41.15;
        Params.k = 0.24;
        Params.t0 = 0; % Using Female fishbase value
        Params.M = 0.3;
        Params.A = 25;
        Params.Amat = 3;
        Params.Af = 4;
        Params.c = 1; % biomass-length constant
        Params.d = 3; % biomass-length exponent species specific from
        % Kaplan et al. 2019
        %Params.d = 3; % biomass-length exponent (general)
                Params.exp = 1; % adjust this to 1.18; in order to perform hyperallometry adjustment from Barneche et al. 2018/Marshall et al. 2019

        % Recruitment rate in Open scenario
        Params.R = 1;
        
        % Some derived quantities:
        Params.Ages = 1:Params.A;
        Params.isMat = Params.Ages>=Params.Amat; % logical - is this age mature?
        Params.isFish = Params.Ages>=Params.Af; % logical - is this age fished?
        % Length-at-age
        Params.Len = Params.Linf*(1 - exp(-(Params.k*(Params.Ages-Params.t0))));
        % Biomass-at-age
        Params.BiomassAge = Params.c.*Params.Len.^Params.d;
        % Stable age distribution
        Params.SAD = exp(-Params.M.*(Params.Ages-1));
        % Lifetime Egg Production (LEP, unfished):
        Params.LEP0 = ((Params.BiomassAge(:).^Params.exp).*Params.isMat(:))'*Params.SAD(:);
        Params.EP0 = (Params.BiomassAge(:).^Params.exp).*Params.isMat(:);
        % Generation time
        Params.Tg = ((Params.BiomassAge(:).^Params.exp).*Params.isMat(:).*(1:Params.A)')'*Params.SAD(:)/Params.LEP0;
        
        % F vs. LEP. We will do this in a quick vectorized way
        F = linspace(0,2,1e2); % 
        %F = linspace(0,0.5,1e2); % 
        Fmat = repmat(F(:)',[Params.A,1]); % matrix of all possible Fs (on cols) for each age (on rows)
        Fmat(Params.Ages<Params.Af,:) = 0; % no fishing below age Af
        Agemat = repmat(Params.Ages(:),[1,length(F)]); % matrix of ages
        SADmat = exp(-(Params.M+Fmat).*(Agemat-1)); % stable age distribution given each F
        LEPmat = sum(SADmat.*repmat((Params.BiomassAge(:).^Params.exp).*Params.isMat(:),[1,length(F)])); % SAD * eggs-at-age = LEP 
        Params.Fs = F;
        Params.FvLEP = LEPmat./Params.LEP0; % sum of SAD
    case 'California sheephead'
        % Define the basic parameters
        Params = struct([]);
        Params(1).Linf = 46.7 ;
        Params.k = 0.18;
        Params.t0 = 0; % Using Female fishbase value
        Params.M = 0.25;
        Params.A = 53;
        Params.Amat = 4;
        Params.Af = 6;
        Params.c = 1; % biomass-length constant
        %Params.d = 2.86; % biomass-length exponent species specific from
        % Kaplan et al. 2019
        Params.d = 3; % biomass-length exponent (general)
        %Params.d = 3.17;
                Params.exp = 1; % adjust this to 1.18; in order to perform hyperallometry adjustment from Barneche et al. 2018/Marshall et al. 2019

        % Recruitment rate in Open scenario
        Params.R = 1;
        
        % Some derived quantities:
        Params.Ages = 1:Params.A;
        Params.isMat = Params.Ages>=Params.Amat; % logical - is this age mature?
        Params.isFish = Params.Ages>=Params.Af; % logical - is this age fished?
        % Length-at-age
        Params.Len = Params.Linf*(1 - exp(-(Params.k*(Params.Ages-Params.t0))));
        % Biomass-at-age
        Params.BiomassAge = Params.c.*Params.Len.^Params.d;
        % Stable age distribution
        Params.SAD = exp(-Params.M.*(Params.Ages-1));
        % Lifetime Egg Production (LEP, unfished):
        Params.LEP0 = ((Params.BiomassAge(:).^Params.exp).*Params.isMat(:))'*Params.SAD(:);
        Params.EP0 = (Params.BiomassAge(:).^Params.exp).*Params.isMat(:);
        % Generation time
        Params.Tg = ((Params.BiomassAge(:).^Params.exp).*Params.isMat(:).*(1:Params.A)')'*Params.SAD(:)/Params.LEP0;
        
        % F vs. LEP. We will do this in a quick vectorized way
        F = linspace(0,2,1e2); % 
        %F = linspace(0,0.5,1e2); % 
        Fmat = repmat(F(:)',[Params.A,1]); % matrix of all possible Fs (on cols) for each age (on rows)
        Fmat(Params.Ages<Params.Af,:) = 0; % no fishing below age Af
        Agemat = repmat(Params.Ages(:),[1,length(F)]); % matrix of ages
        SADmat = exp(-(Params.M+Fmat).*(Agemat-1)); % stable age distribution given each F
        LEPmat = sum(SADmat.*repmat((Params.BiomassAge(:).^Params.exp).*Params.isMat(:),[1,length(F)])); % SAD * eggs-at-age = LEP 
        Params.Fs = F;
        Params.FvLEP = LEPmat./Params.LEP0; % sum of SAD
    case 'Red sea urchin'
        % Define the basic parameters
        Params = struct([]);
        Params(1).Linf = 33.15;
        Params.k = 0.23;
        Params.t0 = 0; % CHECK ON THIS VALUE IN FISHBASE
        Params.M = 0.07;
        Params.A = 100;
        Params.Amat = 3;
        Params.Af = 6;
        Params.c = 1; % biomass-length constant
        %Params.d = 2.17; % biomass-length exponent species specific from
        % Kaplan et al. 2019
        Params.d = 3; % biomass-length exponent (general)
        %Params.d = 2.17;
                Params.exp = 1; % adjust this to 1.18; in order to perform hyperallometry adjustment from Barneche et al. 2018/Marshall et al. 2019

        % Recruitment rate in Open scenario
        Params.R = 1;
        
        % Some derived quantities:
        Params.Ages = 1:Params.A;
        Params.isMat = Params.Ages>=Params.Amat; % logical - is this age mature?
        Params.isFish = Params.Ages>=Params.Af; % logical - is this age fished?
        % Length-at-age
        Params.Len = Params.Linf*(1 - exp(-(Params.k*(Params.Ages-Params.t0))));
        % Biomass-at-age
        Params.BiomassAge = Params.c.*Params.Len.^Params.d;
        % Stable age distribution
        Params.SAD = exp(-Params.M.*(Params.Ages-1));
        % Lifetime Egg Production (LEP, unfished):
        Params.LEP0 = ((Params.BiomassAge(:).^Params.exp).*Params.isMat(:))'*Params.SAD(:);
        Params.EP0 = (Params.BiomassAge(:).^Params.exp).*Params.isMat(:);
        % Generation time
        Params.Tg = ((Params.BiomassAge(:).^Params.exp).*Params.isMat(:).*(1:Params.A)')'*Params.SAD(:)/Params.LEP0;
        
        % F vs. LEP. We will do this in a quick vectorized way
        F = linspace(0,2,1e2); % 
        %F = linspace(0,0.5,1e2); % 
        Fmat = repmat(F(:)',[Params.A,1]); % matrix of all possible Fs (on cols) for each age (on rows)
        Fmat(Params.Ages<Params.Af,:) = 0; % no fishing below age Af
        Agemat = repmat(Params.Ages(:),[1,length(F)]); % matrix of ages
        SADmat = exp(-(Params.M+Fmat).*(Agemat-1)); % stable age distribution given each F
        LEPmat = sum(SADmat.*repmat((Params.BiomassAge(:).^Params.exp).*Params.isMat(:),[1,length(F)])); % SAD * eggs-at-age = LEP 
        Params.Fs = F;
        Params.FvLEP = LEPmat./Params.LEP0; % sum of SAD
   case 'Kelp bass'
        % Define the basic parameters
        Params = struct([]);
        Params(1).Linf = 69.8 ;
        Params.k = 0.06;
        Params.t0 = 0; % Using Female fishbase value
        Params.M = 0.18;
        Params.A = 50;
        Params.Amat = 3;
        Params.Af = 6;
        Params.c = 1; % biomass-length constant
        %Params.d = 3.27; % biomass-length exponent species specific from
        % Kaplan et al. 2019
        Params.d = 3; % biomass-length exponent (general)
        %Params.d = 3.27;
                Params.exp = 1; % adjust this to 1.18; in order to perform hyperallometry adjustment from Barneche et al. 2018/Marshall et al. 2019

        % Recruitment rate in Open scenario
        Params.R = 1;
        
        % Some derived quantities:
        Params.Ages = 1:Params.A;
        Params.isMat = Params.Ages>=Params.Amat; % logical - is this age mature?
        Params.isFish = Params.Ages>=Params.Af; % logical - is this age fished?
        % Length-at-age
        Params.Len = Params.Linf*(1 - exp(-(Params.k*(Params.Ages-Params.t0))));
        % Biomass-at-age
        Params.BiomassAge = Params.c.*Params.Len.^Params.d;
        % Stable age distribution
        Params.SAD = exp(-Params.M.*(Params.Ages-1));
        % Lifetime Egg Production (LEP, unfished):
        Params.LEP0 = ((Params.BiomassAge(:).^Params.exp).*Params.isMat(:))'*Params.SAD(:);
        Params.EP0 = (Params.BiomassAge(:).^Params.exp).*Params.isMat(:);
        % Generation time
        Params.Tg = ((Params.BiomassAge(:).^Params.exp).*Params.isMat(:).*(1:Params.A)')'*Params.SAD(:)/Params.LEP0;
        
        % F vs. LEP. We will do this in a quick vectorized way
        F = linspace(0,2,1e2); % 
        %F = linspace(0,0.5,1e2); % 
        Fmat = repmat(F(:)',[Params.A,1]); % matrix of all possible Fs (on cols) for each age (on rows)
        Fmat(Params.Ages<Params.Af,:) = 0; % no fishing below age Af
        Agemat = repmat(Params.Ages(:),[1,length(F)]); % matrix of ages
        SADmat = exp(-(Params.M+Fmat).*(Agemat-1)); % stable age distribution given each F
        LEPmat = sum(SADmat.*repmat((Params.BiomassAge(:).^Params.exp).*Params.isMat(:),[1,length(F)])); % SAD * eggs-at-age = LEP 
        Params.Fs = F;
        Params.FvLEP = LEPmat./Params.LEP0; % sum of SAD
   case 'Olive rockfish'
        % Define the basic parameters
        Params = struct([]);
        Params(1).Linf = 33.62;
        Params.k = 0.23;
        Params.t0 = 0; % Using Female fishbase value
        Params.M = 0.14;
        Params.A = 30;
        Params.Amat = 4;
        Params.Af = 3;
        Params.c = 1; % biomass-length constant
        %Params.d = 2.968; % biomass-length exponent species specific from
        % Kaplan et al. 2019
        Params.d = 3; % biomass-length exponent (general)
        %Params.d = 2.968;
                Params.exp = 1; % adjust this to 1.18; in order to perform hyperallometry adjustment from Barneche et al. 2018/Marshall et al. 2019

        % Recruitment rate in Open scenario
        Params.R = 1;
        
        % Some derived quantities:
        Params.Ages = 1:Params.A;
        Params.isMat = Params.Ages>=Params.Amat; % logical - is this age mature?
        Params.isFish = Params.Ages>=Params.Af; % logical - is this age fished?
        % Length-at-age
        Params.Len = Params.Linf*(1 - exp(-(Params.k*(Params.Ages-Params.t0))));
        % Biomass-at-age
        Params.BiomassAge = Params.c.*Params.Len.^Params.d;
        % Stable age distribution
        Params.SAD = exp(-Params.M.*(Params.Ages-1));
        % Lifetime Egg Production (LEP, unfished):
        Params.LEP0 = ((Params.BiomassAge(:).^Params.exp).*Params.isMat(:))'*Params.SAD(:);
        Params.EP0 = (Params.BiomassAge(:).^Params.exp).*Params.isMat(:);
        % Generation time
        Params.Tg = ((Params.BiomassAge(:).^Params.exp).*Params.isMat(:).*(1:Params.A)')'*Params.SAD(:)/Params.LEP0;
        
        % F vs. LEP. We will do this in a quick vectorized way
        F = linspace(0,2,1e2); % 
        %F = linspace(0,0.5,1e2); % 
        Fmat = repmat(F(:)',[Params.A,1]); % matrix of all possible Fs (on cols) for each age (on rows)
        Fmat(Params.Ages<Params.Af,:) = 0; % no fishing below age Af
        Agemat = repmat(Params.Ages(:),[1,length(F)]); % matrix of ages
        SADmat = exp(-(Params.M+Fmat).*(Agemat-1)); % stable age distribution given each F
        LEPmat = sum(SADmat.*repmat((Params.BiomassAge(:).^Params.exp).*Params.isMat(:),[1,length(F)])); % SAD * eggs-at-age = LEP 
        Params.Fs = F;
        Params.FvLEP = LEPmat./Params.LEP0; % sum of SAD
    case 'Black and yellow rockfish'
        % Define the basic parameters
        Params = struct([]);
        Params(1).Linf = 24.95 ;
        Params.k = 0.23;
        Params.t0 = 0; % Using Female fishbase value
        Params.M = 0.2;
        Params.A = 30;
        Params.Amat = 5;
        Params.Af = 14;
        Params.c = 1; % biomass-length constant
        %Params.d = 3.114; % biomass-length exponent species specific from
        % Kaplan et al. 2019
        Params.d = 3; % biomass-length exponent (general)
        %Params.d = 3.114*1.18; % biomass-length exponent species specific as 
        % above and adjustment given hyperallometry adjustment from Barneche et al. 2018/Marshall et al. 2019
                Params.exp = 1; % adjust this to 1.18; in order to perform hyperallometry adjustment from Barneche et al. 2018/Marshall et al. 2019

        % Recruitment rate in Open scenario
        Params.R = 1;
        
        % Some derived quantities:
        Params.Ages = 1:Params.A;
        Params.isMat = Params.Ages>=Params.Amat; % logical - is this age mature?
        Params.isFish = Params.Ages>=Params.Af; % logical - is this age fished?
        % Length-at-age
        Params.Len = Params.Linf*(1 - exp(-(Params.k*(Params.Ages-Params.t0))));
        % Biomass-at-age
        Params.BiomassAge = Params.c.*Params.Len.^Params.d;
        % Stable age distribution
        Params.SAD = exp(-Params.M.*(Params.Ages-1));
        % Lifetime Egg Production (LEP, unfished):
        Params.LEP0 = ((Params.BiomassAge(:).^Params.exp).*Params.isMat(:))'*Params.SAD(:);
        Params.EP0 = (Params.BiomassAge(:).^Params.exp).*Params.isMat(:);
        % Generation time
        Params.Tg = ((Params.BiomassAge(:).^Params.exp).*Params.isMat(:).*(1:Params.A)')'*Params.SAD(:)/Params.LEP0;
        
        % F vs. LEP. We will do this in a quick vectorized way
        F = linspace(0,2,1e2); % 
        %F = linspace(0,0.5,1e2); % 
        Fmat = repmat(F(:)',[Params.A,1]); % matrix of all possible Fs (on cols) for each age (on rows)
        Fmat(Params.Ages<Params.Af,:) = 0; % no fishing below age Af
        Agemat = repmat(Params.Ages(:),[1,length(F)]); % matrix of ages
        SADmat = exp(-(Params.M+Fmat).*(Agemat-1)); % stable age distribution given each F
        LEPmat = sum(SADmat.*repmat((Params.BiomassAge(:).^Params.exp).*Params.isMat(:),[1,length(F)])); % SAD * eggs-at-age = LEP 
        Params.Fs = F;
        Params.FvLEP = LEPmat./Params.LEP0; % sum of SAD
    otherwise
        error('Unknown species requested')
end 