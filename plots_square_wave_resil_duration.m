function [Species, N95, Nmin, intN, Mvarab,red]=plots_square_wave_resil_duration(Conn_scenario,DD_scenario)

%All 16 species
Species = {'China rockfish','Copper rockfish','Vermilion rockfish',...
 'Yellowtail rockfish','Black rockfish','Blue rockfish','Brown rockfish',...
 'Olive rockfish','Kelp bass','Black and yellow rockfish','Gopher rockfish',...
 'Kelp rockfish','California scorpionfish','Lingcod','Cabezon','Kelp greenling'};

red = 0.5; % pulse severity
%red = [0.2:0.1:0.7]; % varying severity of the pulse

%dur = 5; %duration of the pulse
dur = [1:1:5]; % varying duration of pulse

flep = 1; % fishing pressure in patch
%flep = [0:0.1:1]; % varying FLEP

FLEPs=flep;
Reduction=log(red); %0.5
Duration=dur; 
Delay=NaN; %for more than one pulse
T = 500; % total duration of simulation
Ti = 100; % initial pre-disturbance period
Thresh = 0.95; % proportion of initial abundance

%Initiate parameters
N=zeros(length(Species),length(dur),T);
B=zeros(length(Species),length(dur),T);
Y=zeros(length(Species),length(dur),T);

N95=zeros(length(Species),length(dur));%,length(FLEPs))
B95=zeros(length(Species),length(dur));%length(FLEPs)

Nmin=zeros(length(Species),length(dur));
Bmin=zeros(length(Species),length(dur));
Ymin=zeros(length(Species),length(dur));

subNtmp = zeros(length(Species),length(dur),length((Ti-1):T));%length((Ti-1):T));
invsubNtmp = zeros(length(Species),length(dur),length((Ti-1):T));
subBtmp = zeros(length(Species),length(dur),length((Ti-1):T));
invsubBtmp = zeros(length(Species),length(dur),length((Ti-1):T));
subYtmp = zeros(length(Species),length(dur),length((Ti-1):T));
invsubYtmp = zeros(length(Species),length(dur),length((Ti-1):T));

intN = zeros(length(Species),length(dur));
intB = zeros(length(Species),length(dur));
Bpost = zeros(length(Species),length(dur),T-Ti-1);
%Npost = zeros(length(Species),length(dur),400); %Npost = zeros(length(Species),length(dur),(T-Ti)-(max(Duration)-1));

Ypost = zeros(length(Species),length(dur),T-Ti-1);
intY = zeros(length(Species),length(dur));
one_Ntmp=zeros(length(Species),length(dur),491);
log_Ntmp=zeros(length(Species),length(dur),491);
Tg =zeros(length(Species),length(dur));
Mvarab =zeros(length(Species),length(dur));
minus1 = zeros(length(Species),length(dur),500);
Nrate=zeros(length(Species),length(dur));

Ntmp=zeros(length(Species),length(dur),T);
Btmp=zeros(length(Species),length(dur),T);
Ytmp=zeros(length(Species),length(dur),T);

for i=1:length(Species);
    for j=1:length(Duration);  
        %for k=1:length(Duration);

    [Noise, N, B, E, EP, Y, Params]=square_wave_resil(Species{i},FLEPs,Reduction,Duration(j),Delay,T,Ti,Conn_scenario,DD_scenario);

Tg(i)=Params.Tg;
Mvarab(i)=Params.M;

    Ntmp(i,j,:) = squeeze(sum(N(:,1,:)));    
    Btmp(i,j,:) = squeeze(sum(B(:,1,:)));
    Ntmp(i,j,:) = Ntmp(i,1,:)./Ntmp(i,1,Ti-1); % standardized abundance trajectory
    Btmp(i,j,:) = Btmp(i,1,:)./Btmp(i,1,Ti-1); % standardized biomass trajectory
    Ytmp(i,j,:) = Y(:)./Y(end); % standardized yield trajectory

    %Calculate impact 
    Nmin(i,j,:) = min(Ntmp(i,j,:)); %IMPACT
    Bmin(i,j,:) = min(Btmp(i,j,(Ti+Duration(j)):end)); %IMPACT
    Ymin(i,j,:) = min(Ytmp(i,j,(Ti+Duration(j)):end)); %IMPACT

    %miN = (Ntmp(:,103)-Ntmp(:,101))/(103-101);
    %Calculate area over/under the curve due to impact
    subNtmp(i,j,:) = Ntmp(i,j,(Ti-1):end); % for abundance 
    subBtmp(i,j,:) = Btmp(i,j,(Ti-1):end); % for biomass
    subYtmp(i,j,:) = Ytmp(i,j,(Ti-1):end); % for yield
    xLast = find(Ntmp(i,j,:) < 1, 1, 'last');
    
    % Calculate integral (area between curve and 1)
    int1(i,j) = max(T);
    intN(i,j) = int1(i,j)-sum(abs(Ntmp(i,j,1:xLast))); 

    %intN(i)=(area(Ntmp(i,98:end))<=1);%-300)); %perturbation area
    %intB(i)=trapz(invsubBtmp(i,:));%-300; %perturbation area
    %intY(i)=trapz(invsubYtmp(i,:));%-300; %perturbation area

%% % Now calculate trajectory recovery time
Npost(:,:,:) = zeros(length(Species),length(Duration(j)),T-(Ti-Duration(j))); %Npost = zeros(length(Species),length(dur),(T-Ti)-(max(Duration)-1));
%keyboard
% For abundance
Npost(i,j,:) = sum(N(:,1,(Ti+Duration(j)):end)); % just the trajectory after disturbance
Ninit = squeeze(sum(N(:,1,Ti))); % baseline pre-disturbance
% for f = 1:length(FLEPs)
keyboard
     %Ndiff(i,j,:) = abs(sum(N(:,1,(Ti+Duration(j)+j):end))-Ninit)/Ninit;
     Ndiff(i,j,:) = abs(Npost(i,j,:) - Ninit)/Ninit; %was divided by Ninit...not sure why
     %N95(i) = find(Ndiff(i,:)>=(1-Thresh),1,'last');
     if all(Ndiff(i,j,:) <= 1- Thresh) % if it never gets below the threshold
     N95(i) = 0;
     elseif all(Ndiff(i,j,:) > 1- Thresh) % if it never recovers enough
         N95(i,j) = T-Ti;
     else % otherwise, use the standard procedure
        N95(i,j) = find(Ndiff(i,j,:)>=(1-Thresh),1,'last');
     end
      
     Nrate(i,j) = max(gradient(log(1-(sum(squeeze(N(:,1,(Ti+Duration(j)):end)))))));

% %text(Nmin_10,N95_10,Species','FontSize',15)
% box on
% hold on
        end
    end
%end

Mvar =repmat(Mvarab(:,1),1,6)
%% Ingrish and Bahn biplots


%col={'ro', 'r*', 'rd', 'bo', 'b*', 'bd', 'ko', 'k*', 'kd','go', 'g*', 'gd', 'mo', 'm*', 'md', 'co'};
%subplot(1,2,2)
%for i=1:length(Species);
    %for j=1:length(Reduction);
    %surfcolors = flipud(jet(length(Mvarab)));
    %hold on
    %contourf(Nmin, N95, Mvarab)
    surf(Nmin,Mvar,N95)
    shading interp;
    colormap(cool(20));
%plot(Nmin(i,1:6),N95(i,1:6),'-o','color',surfcolors(i))
    %scatter(Nmin(i,j),N95(i,j),i*100,Mvarab(i),col{i})%'d','filled')
    xlim([0,1])
    ylim([0.05, 0.35])
    zlim([0, 30])
    
    ylabel('Natural Mortality Rate')
    xlabel('N impact')
    zlabel('Recovery time (t)')

    %end
%end
%hold on
%legend(col,Species)
%set(gca,'Colormap',surfcolors)
%h=colorbar('FontSize',12)
%h.Label.String = 'Natural Mortality (M)';


% %Ym_10= abs(Ym(:,1));
% %Nrt_10= Nrt(:,1);
% maxNgrad_10= abs(maxNgrad(:,1));
% gmeanNgrad_10 = gmeanNgrad(:,1);

% X = [N95_10, intN_10];
% R=corrcoef(X)
%rsq = R(2,1)^2
%keyboard

%for generation time
% % hold on
% %subplot(1,2,1)
% figure(2)
% scatter(Nmin_10,Tg,100,intN_10,'filled')%'d','filled'
%ylim([3,12.5])
%xlim([0.75, 1])
% xlabel('N impact (%)')
% ylabel('Generation time (t)')
% h=colorbar
% h.Label.String = 'Perturbation N';
%text(Nmin_10,Tg,Species')
% box on
%scatterplot


  
