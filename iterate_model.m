function [N, B, Y, E, EP, C, Yextra] = iterate_model(Params,L,F,N0,Y0,C0,B0,T,Conn_scenario,DD_scenario, MPA_frac,Noise,doAdultNoise,AdultNoise)

% Iterate a Leslie-matrix-based model
% Note that F is a vector with a value for each patch
% Updated 25 June 2019 by WW to correct yield calculations

%Add in noise to recruitment;
if ~exist('Noise','var')
    Noise = ones(T,1);
else

    if length(Noise) < T % if sending parameters to create colored noise
         Noise = exp(normrnd(Noise(1),Noise(2),T,1)); % white noise (ok for now...add other features later)
        % Noise = [mean, sd] of white noise process
        % Can add other parameters to indicate what type of noise &
        % parameters (red, 1/f, etc)
    else
        %In this case, Noise is the vector of disturbance at each t
         Noise=exp(Noise); % This + uncomment section above, and hard coding the noise change in the
        %closed patch seems to work...
    end
end

if ~exist('doAdultNoise','var')
    doAdultNoise = false;
end
 if ~exist('AdultNoise','var')
    AdultNoise = ones(T,1);
 end

if size(F,2)==1
    F = repmat(F(:),[1,T]); % may need to be adjusted for multi-patch...
end


% Initialize
N = zeros(size(N0,1),size(N0,2),T);
B = N;
%E = N;
E=zeros(size(N0,2),T);
EP = zeros(size(N0,1),T);
D_BM = N(1:end-1,:,:); % has one fewer dimension
Y = zeros(size(N0,2),T);
Y(:,1) = Y0;
C = zeros(size(N0,2),T);
C(:,1) = C0;
N(:,:,1)= N0; % initialize density
%B(:,:,1) = N0.*repmat(Params.BiomassAge(:),[1,size(N0,2)]); % initialize biomass
B(:,:,1) = B0;

% Add in extra state variable for 'extra' recruits in fished patch
Nextra = zeros(size(N));
Bextra = zeros(size(B));
Yextra = zeros(size(Y));

% storage for 100 iterations
% random number generation on e^-random number in iterate_model-
% randomizing recruitment given environmental variability.
% Nt0 = zeros(100,Params.A,2); % Abundance with 100 iterations x maximum age x number of patches
% Bt0 = zeros(100,Params.A, 2); % Biomass ditto
% Yt0 = zeros(100,2,T); % Yield ditto
% Ct0 = zeros(100,2,T); % CPUE ditto

switch Conn_scenario
    case {'Open'}
            for t = 2:T
                
                R = Params.R * Noise(t); % noise in recruitment
                 
                for i = 1:size(N0,2) % number of patches
                N(:,i,t) = L(:,:,i,t)*N(:,i,t-1); 
                if doAdultNoise
                   N(:,i,t) = N(:,i,t)*AdultNoise(t);
                end
                N(1,i,t)=R; %N(1,i,t)=R*MPA_frac; % !!!! CHANGE THIS IF MPA_FRAC>0 %% open recruitment
                B(:,i,t)= N(:,i,t).* Params.BiomassAge(:); %(N(1:end,t-1)-N(2:end)).*(Ba(1:end-1)+Ba(2:end))/2;
                
                
                
                if F(i) > 0
                Y(i,t)=sum((B(Params.isFish(1:end-1),i,t).*F(i,t))/(Params.M+F(i,t))); %
                C(i,t)=Y(i,t)./F(i,t); % CPUE
                end
%                 % Now 100 iterations to get at natural variability in
%                 % recruitment.
%                     for j=1:100-1;
%                             Nt0(j,:,t) = L(:,:,i)*N(:,i,t-1); % add in recruits to open population - randomizing from a standard normal distribution one point to simulate environmental stochasticity in recruitment.
%                             Nt0(1,:,t) = Params.R*exp(-randn(1)); % add in recruits to open population - randomizing from a standard normal distribution one point to simulate environmental stochasticity in recruitment.
%                             Bt0(j,:,t) = Nt0(j,i,t).* Params.BiomassAge(:);
%                             %D_BMt0(j,:,t) = Bt0(1:end-1,i,t-1)-Bt0(2:end,i,t);
%                             Yt0(j,:,t) = sum((Bt0(Params.isFish(1:end-1),i,t).*F(i))/(Params.M+F(i)));
%                             %Ct0(j,:,t) = Y(j,i,t)./F(i);
%                     end                    
                % I don't think we want to do the next step
                %CPUE(:,i,t)=(Params.F.*CPUE(:,i,2:t))./sum(CPUE(:,2:t)); % standardized CPUE
                end % end loop over patches
                
                if Noise(t) < 1
               % keyboard
                end
                
            end% end loop over time    
    
    case {'Porous'}
        
                        

        for t = 2:T % time
            
            R = Params.R * Noise(t); % noise in recruitment
            
                for i = 1:size(N0,2) % number of patches
                    N(:,i,t) = L(:,:,i,t)*N(:,i,t-1);
                    if doAdultNoise
                   N(:,i,t) = N(:,i,t)*AdultNoise(t);
                    end
                    Nextra(:,i,t) = L(:,:,i)*Nextra(:,i,t-1);
                    if i==1 % MPA
                    E(i,t)=N(1,i,t); %Eggs produced by mpa
                    
                    % NEED TO SCALE THIS BY HOW MANY ARE PRODUCED AT TIME 0
                    
                    N(1,i,t)=(R*MPA_frac)./N(1,1,1); % what is this for???????
                    N(1,i,t)=(R*MPA_frac); %./N(1,1,1); % why is it changed here?????

                    %N(1,:,t) = (Params.R*exp(-randn(1)))*MPA_frac; % add
                    %in recruits to open population - randomizing from a standard normal distribution one point to simulate environmental stochasticity in recruitment
                    B(:,i,t)= N(:,i,t).* Params.BiomassAge(:); %(N(1:end,t-1)-N(2:end)).*(Ba(1:end-1)+Ba(2:end))/2;
                    %Y(i,t)=sum((B(Params.isFish(1:end-1),i,t).*F(i))/(Params.M+F(i))); %
                    %C(i,t)=Y(i,t)./F(i); % CPUE
                    end
                
                    if i==2 % Fished patch
                        if F(1) == 0 % if patch 1 is a reserve (i.e., F=0 there) only do this if that is true)
                            Etmp = E(1,t)-E(1,1); % scale by initial egg production so we are only taking the overage
                    E_obase=R*(1-MPA_frac)+Etmp*(1-MPA_frac); 
                        else
                            Etmp = 0;
                    E_obase=R*(1-MPA_frac);
                        end
                    N(1,i,t) = E_obase; % add in recruits to protected population
                    
                    Nextra(1,i,t) = Etmp*(1-MPA_frac);
%                                 for j=1:100-1;
%                             Nt0(j,:,t) = L(:,:,i)*N(:,i,t-1); % add in recruits to open population - randomizing from a standard normal distribution one point to simulate environmental stochasticity in recruitment.
%                             Nt0(1,:,t) = E_obase;
%                             Nt0(1,:,t) = (Params.R*exp(-randn(1)))*MPA_frac; % add in recruits to open population - randomizing from a standard normal distribution one point to simulate environmental stochasticity in recruitment.
%                             Bt0(j,:,t) = Nt0(j,i,t).* Params.BiomassAge(:);
%                             %D_BMt0(j,:,t) = Bt0(1:end-1,i,t-1)-Bt0(2:end,i,t);
%                             Yt0(j,:,t) = sum((Bt0(Params.isFish(1:end-1),i,t).*F(i))/(Params.M+F(i)));
%                             Ct0(j,:,t) = Y(j,i,t)./F(i);
%                                 end
                
                    end
                    B(:,i,t)= N(:,i,t).* Params.BiomassAge(:); %(N(1:end,t-1)-N(2:end)).*(Ba(1:end-1)+Ba(2:end))/2; 
                    Bextra(:,i,t) = Nextra(:,i,t).* Params.BiomassAge(:); % 'extra' biomass due to MPA
                    
                    if F(i) > 0
                    Y(i,t)=sum((B(Params.isFish(1:end-1),i,t).*F(i,t))/(Params.M+F(i,t))); %
                    C(i,t)=Y(i,t)./F(i,t); % CPUE
                    
                    Yextra(i,t) = sum((Bextra(Params.isFish(1:end-1),i,t).*F(i,t))/(Params.M+F(i,t))); 
                    
                    end
                end
        end

    case {'Closed'}
                    
        if size(N0,2) == 1
            dispmat = 1;
        elseif size(N0,2) == 2
        dispmat=[MPA_frac,MPA_frac;(1-MPA_frac),(1-MPA_frac)];
        else
            error('Can only have 1 or 2 patches')
        end

        for t = 2:T %this was 2:T
            for i = 1:size(N0,2) % number of patches
                N(:,i,t)=L(:,:,i,t)*N(:,i,t-1);%closed pop has no external recruits
                if doAdultNoise
                   N(:,i,t) = N(:,i,t)*AdultNoise(t);
                end
                %N(1,:,t) = (Params.R*exp(-randn(1)))*MPA_frac; % add noise
                %to recruitment
                %E(i,t)=N(:,i,t).*Noise(t);% %Eggs produced in that patch
%               subjected to a full age structure decrease due to noise
                E(i,t)=N(1,i,t)*Noise(t);% %Eggs produced in that patch
                %E(i,t) = N(1,i,t).*exp(-normrnd(1,0.3)); % white noise
                %E(i,t) = N(1,i,t).*exp(diffENSOnorm(t));
                B(:,i,t)= N(:,i,t).* Params.BiomassAge(:);%(N(1:end,t-1)-N(2:end)).*(Ba(1:end-1)+Ba(2:end))/2;
                EP(:,t) = Params.EP0(:)*Noise(t);
                if F(i) > 0
                Y(i,t)=sum((B(Params.isFish(1:end-1),i,t).*F(i,t))/(Params.M+F(i,t))); %
                C(i,t)=Y(i,t)./F(i,t); % CPUE
                end
             end % end loop over patches
           
            % Larval dispersal & density-dependence 
            Rdis=dispmat*E(:,t);
                 
            
     switch DD_scenario
         case {'BH'}
             a = 0.20; 
            Rdis1=(1/(a*Params.LEP0)).*Rdis./(1+((1/(a*Params.LEP0))/1).*Rdis);
       %  if t == 150; keyboard; end
         case {'Steepness'}
            h = Params.steepness;
            a = (h-0.2)/(0.8*h); % convert steepness to slope
           % a = 1/a; % convert to CRT
            %keyboard
            Rdis1=(1/(a*Params.LEP0)).*Rdis./(1+((1/(a*Params.LEP0))/1).*Rdis);

         otherwise
             Rdis1 = Rdis;
     end % end switch
     
            N(1,:,t)=Rdis1;
            
            for i = 1:size(N0,2)
            B(:,i,t)= N(:,i,t) .* Params.BiomassAge(:); %(N(1:end,t-1)-N(2:end)).*(Ba(1:end-1)+Ba(2:end))/2;
            end 
            %  if F(i) > 0
          %  Y(i,t)=sum((B(Params.isFish(1:end-1),i,t).*F(i))/(Params.M+F(i))); %
          %  C(i,t)=Y(i,t)./F(i); % CPUE
          
                  % % Now 100 iterations to get at natural variability in
%                 % recruitment.
%                     for j=1:100-1;
%                             Nt0(j,:,t) = L(:,:,i)*N(:,i,t-1); % add in recruits to open population - randomizing from a standard normal distribution one point to simulate environmental stochasticity in recruitment.
%                             Nt0(1,:,t) = Params.R*exp(-randn(1)); % add in recruits to open population - randomizing from a standard normal distribution one point to simulate environmental stochasticity in recruitment.
%                             Bt0(j,:,t) = Nt0(j,i,t).* Params.BiomassAge(:);
%                             %D_BMt0(j,:,t) = Bt0(1:end-1,i,t-1)-Bt0(2:end,i,t);
%                             Yt0(j,:,t) = sum((Bt0(Params.isFish(1:end-1),i,t).*F(i))/(Params.M+F(i)));
%                             %Ct0(j,:,t) = Y(j,i,t)./F(i);
%                     end 
    

        end % end loop over time
        
% % Loop a few times to simulate infinite (or semi-infinite) coastline:
% switch type
%     case 'inf'
% for i = -50:50
% D = D + normcdf(Dist+0.5+x*i,Mu,Sig) - normcdf(Dist-0.5+x*i,Mu,Sig);    
% end
%     case 'semi'
% 
% D = normcdf(Dist+0.5,Mu,Sig) - normcdf(Dist-0.5,Mu,Sig);    

            %end % end loop over patches
        %end% end loop over time 
         end %end switch over connectivity

         
