function runme_pulse(Fig1,Fig2,Fig3,Fig4a,Fig4b,Fig5,Fig6,FigS1)

% Master runme file for analysis of response of age-structured populations
% to pulse disturbances

% Inputs: arguments are logical scalars indicating whether to create
% particular figures

% Outputs: figures

if ~exist('Fig1','var'); Fig1 = false; end % example pulse-response
if ~exist('Fig2','var'); Fig2 = false; end % example influence functions
if ~exist('Fig3','var'); Fig3 = false; end % example pulse-response
if ~exist('Fig3x','var'); Fig3x = false; end % example pulse-response as a convolution
if ~exist('Fig4a','var'); Fig4a = false; end % bivariate plot, multiple species (open or closed, closed is 4c)
if ~exist('Fig4b','var'); Fig4b = false; end % bivariate plot, multiple species (open or closed)
if ~exist('Fig5','var'); Fig5 = false; end % harvest, one species
if ~exist('Fig6','var'); Fig6 = false; end % harvest, Pacific cod, harvest stops
if ~exist('Fig6b','var'); Fig6b = true; end % P cod with & without recruitment effects
if ~exist('Fig7','var'); Fig7 = false; end % Pacific cod, no harvest intervention (deprecated)
if ~exist('FigS1','var'); FigS1 = false; end % illustration of convolution approach, multipanel

%% Figure 1: Example pulse-response figure 
if Fig1
    
Species = 'Black rockfish';
FLEPs = 0.3; % heavy fishing
Reduction = log(0.1); % 50% reduction in juvenile survival
Duration = 5; % 5 year disturbance
Delay = NaN; % only used if multiple disturbances (it gives the interval)
T = 200; % total time
Ti = 100; % time prior to disturbance
Conn_scenario = 'Closed'; % Connectivity: open (linear subsidy) or closed
DD_scenario = 'BH'; % density dependence: Bev-Holt, or none?


[~, N, B, ~, ~, ~, ~,Resist,RT] = square_pulse_resil(Species,FLEPs,Reduction,Duration,Delay,T,Ti,...
                                                     Conn_scenario, DD_scenario,'B',0.95);

Nplot = sum(squeeze(N));
Nplot = Nplot/Nplot(Ti-1);

Bplot = sum(squeeze(B));
Bplot = Bplot/Bplot(Ti-1);

% Create plot
resil_trajec_plot([Nplot(:),Bplot(:)],Duration,T,Ti,'Relative biomass')

end %% end if Fig1

%% Figure 2: Influence functions
if Fig2
Species = 'Black rockfish';
FLEP = 0.33; % from Dick et al assessment
Mfactor = [0.5, 1, 1.5];

infl_fxn_plot(Species,FLEP,Mfactor)
end
%% End Fig 2



%% Figure 3: Example pulse-response figures with multiple durations or intensities, multipanel
if Fig3
    
Species = 'Black rockfish';
FLEPs = 1; % no fishing
Reductions = log([0.9, 0.7, 0.5]); % 50% reduction in juvenile survival
Duration = 10; % 10 year disturbance
Mfactors = [1.5 1 0.5]; % multipliers on M
Durations = [5 10 15];
Delay = NaN; % only used if multiple disturbances (it gives the interval)
T = 200; % total time
Ti = 100; % time prior to disturbance
Conn_scenario = 'Open'; % Connectivity: open (linear subsidy) or closed
DD_scenario = 'none'; % density dependence: Bev-Holt, or none?
AdultMortality = false; % does disturbance affect all age classes or just age 0?

Col = [0 0 0; 0.2 0.2 0.9; 0.9, 0.2 0.2]; % K B R

figure(2)
clf
set(gcf,'units','cent','position',[10,10,21,12])...

for i = 1:length(Reductions)
    [~, N, B, ~, ~, ~, ~,Resist,RT] = square_pulse_resil(Species,FLEPs,Reductions(i),Duration,Delay,T,Ti,...
                                                     Conn_scenario, DD_scenario,'B',0.95,1,AdultMortality);
       
Nplot_tmp = sum(squeeze(N));
Nplot_tmp = Nplot_tmp/Nplot_tmp(Ti-1);
Nplot(:,i) = Nplot_tmp(:);

Bplot_tmp = sum(squeeze(B));
Bplot_tmp = Bplot_tmp/Bplot_tmp(Ti-1);
Bplot(:,i) = Bplot_tmp(:);

end

sh(1) = subplot(2,3,1);
hold on
resil_trajec_plot(Nplot,Duration,T,Ti,'Relative abundance',Col)
xl = get(gca,'xlim');
set(gca,'xlim',[xl(1),50])

sh(4) = subplot(2,3,4);
hold on
resil_trajec_plot(Bplot,Duration,T,Ti,'Relative biomass',Col)
xl = get(gca,'xlim');
set(gca,'xlim',[xl(1),50])

for i = 1:length(Durations)
    [~, N, B, ~, ~, ~, ~,Resist,RT] = square_pulse_resil(Species,FLEPs,Reductions(2),Durations(i),Delay,T,Ti,...
                                                     Conn_scenario, DD_scenario,'B',0.95,1);

                                                 
Nplot_tmp = sum(squeeze(N));
Nplot_tmp = Nplot_tmp/Nplot_tmp(Ti-1);
Nplot(:,i) = Nplot_tmp(:);

Bplot_tmp = sum(squeeze(B));
Bplot_tmp = Bplot_tmp/Bplot_tmp(Ti-1);
Bplot(:,i) = Bplot_tmp(:);                                               
                                                 
end

sh(2) = subplot(2,3,2);
hold on
% Consider specifying plot window size, etc.
resil_trajec_plot(Nplot,Duration,T,Ti,'Relative abundance',Col)
xl = get(gca,'xlim');
set(gca,'xlim',[xl(1),50])

sh(5) = subplot(2,3,5);
hold on
% Consider specifying plot window size, etc.
resil_trajec_plot(Bplot,Duration,T,Ti,'Relative biomass',Col)
xl = get(gca,'xlim');
set(gca,'xlim',[xl(1),50])


for i = 1:length(Mfactors)
    [~, N, B, ~, ~, ~, ~,Resist,RT] = square_pulse_resil(Species,FLEPs,Reductions(2),Duration,Delay,T,Ti,...
                                                     Conn_scenario, DD_scenario,'B',0.95,Mfactors(i));

                                                 
Nplot_tmp = sum(squeeze(N));
Nplot_tmp = Nplot_tmp/Nplot_tmp(Ti-1);
Nplot(:,i) = Nplot_tmp(:);

Bplot_tmp = sum(squeeze(B));
Bplot_tmp = Bplot_tmp/Bplot_tmp(Ti-1);
Bplot(:,i) = Bplot_tmp(:);                                               
                                                 
end

sh(3) = subplot(2,3,3);
hold on
% Consider specifying plot window size, etc.
resil_trajec_plot(Nplot,Duration,T,Ti,'Relative abundance',Col)
xl = get(gca,'xlim');
set(gca,'xlim',[xl(1),50])

sh(6) = subplot(2,3,6);
hold on
% Consider specifying plot window size, etc.
resil_trajec_plot(Bplot,Duration,T,Ti,'Relative biomass',Col)
xl = get(gca,'xlim');
set(gca,'xlim',[xl(1),50])

set(sh,'xcolor','k','ycolor','k','ylim',[0.6 1])

end %% end if Fig2

%% Fig 3X
% Examples of convolution
if Fig3x
    
end % end if Fig3x




%% Fig 4a
% one FLEP, range of species
if Fig4a

Species = {'China rockfish','Copper rockfish','Vermilion rockfish',...
 'Yellowtail rockfish','Black rockfish','Blue rockfish','Brown rockfish',...
 'Olive rockfish','Kelp bass','Black and yellow rockfish','Gopher rockfish',...
 'Kelp rockfish','California scorpionfish','Lingcod','Cabezon','Kelp greenling'};
FLEPs = 1; % no fishing
Reductions = log([0.4, 0.3, 0.2]); % 50% reduction in juvenile survival
%Reductions = log(0.0001);
Duration = 5; % 5 year disturbance
Delay = NaN; % only used if multiple disturbances (it gives the interval)
T = 200; % total time
Ti = 100; % time prior to disturbance

Conn_scenario = 'Closed'; % Connectivity: open (linear subsidy) or closed
DD_scenario = 'BH'; % density dependence: Bev-Holt, or none?
AdultMortality = false; % does disturbance affect all age classes or just age 0?


for r = 1:length(Reductions)
for s = 1:length(Species)
[~, N, ~, ~, ~, ~, Params_tmp, Resist(s,r),RT(s,r)] = square_pulse_resil(Species{s},FLEPs,...
                                                 Reductions(r),Duration,Delay,T,Ti,Conn_scenario,...
                                                 DD_scenario,'B',0.95,1,AdultMortality);
Ms(s) = Params_tmp.M;
hold on
end
end

% IB plot
% create color scale for M
Coltmp = flipud(winter(100));
Map = linspace(min(Ms),max(Ms),100);
for m = 1:length(Ms)
    MapM = abs(Map-Ms(m));
    IndM(m) = find(MapM==min(MapM),1);
end

Col = Coltmp(IndM,:);
Symb = {'o','^','d'};


IB_plot(Resist,RT,T-Ti,Symb,Col,Coltmp,Map,'Natural mortality (\itM\rm, y^-^1)')
xlim([0 0.5])
ylim([0 40])
text(0.05,35,'Low M','color',Col(3,:),'fontsize',12)
text(0.4,15,'High M','color',Col(end,:),'fontsize',12)
%text(Resist(1,2)*0.9,RT(1,2),'Closed','color','r','fontsize',12)




%Conn_scenario = 'Closed'; % Connectivity: open (linear subsidy) or closed
%DD_scenario = 'BH'; % density dependence: Bev-Holt, or none?
    
%[~, ~, ~, ~, ~, ~, ~, Resistc,RTc] = square_pulse_resil(Species,FLEPs,Reduction,Duration,Delay,T,Ti,Conn_scenario, DD_scenario);

end % end if Fig4a
%%

%% Fig4b: Disturbance affects more than first age class
if Fig4b
Species = 'Black rockfish';
FLEPs = 1.0; % heavy fishing
Reduction = log(0.5); % 50% reduction in juvenile survival
Duration = 5; % 5 year disturbance
Delay = NaN; % only used if multiple disturbances (it gives the interval)
T = 200; % total time
Ti = 100; % time prior to disturbance
Conn_scenario = 'Closed'; % Connectivity: open (linear subsidy) or closed
DD_scenario = 'BH'; % density dependence: Bev-Holt, or none?


[~, N, B, ~, ~, ~, ~,Resist,RT] = square_pulse_resil(Species,FLEPs,Reduction,Duration,Delay,T,Ti,...
                                                     Conn_scenario, DD_scenario,'B',0.95,true);

                                                 Nplot = sum(squeeze(N));
Nplot = Nplot/Nplot(Ti-1);

Bplot = sum(squeeze(B));
Bplot = Bplot/Bplot(Ti-1);

% Consider specifying plot window size, etc.
resil_trajec_plot([Nplot(:),Bplot(:)],Duration,T,Ti,'Relative biomass')
end % end if Fig4b
%%

%% Fig 5
% one species, range of harvest rates
if Fig5

Species = 'Black rockfish';
FLEPs = 0.3:0.1:1; % no fishing
Reduction = log(0.1); % 50% reduction in juvenile survival
Duration = 5; % 5 year disturbance
Delay = NaN; % only used if multiple disturbances (it gives the interval)
T = 200; % total time
Ti = 100; % time prior to disturbance

Conn_scenario = 'Open'; % Connectivity: open (linear subsidy) or closed
DD_scenario = 'None'; % density dependence: Bev-Holt, or none?

[~, ~, ~, ~, ~, ~, ~, Resisto,RTo] = square_pulse_resil(Species,FLEPs,Reduction,Duration,Delay,T,Ti,Conn_scenario, DD_scenario,'B');

Conn_scenario = 'Closed'; % Connectivity: open (linear subsidy) or closed
DD_scenario = 'BH'; % density dependence: Bev-Holt, or none?
    
[~, ~, ~, ~, ~, ~, ~, Resistc,RTc] = square_pulse_resil(Species,FLEPs,Reduction,Duration,Delay,T,Ti,Conn_scenario, DD_scenario,'B');


% IB plot
% organize into arrays for plotting
Resist = nan(length(FLEPs),2);
RT = Resist;
Resist(:,1) = Resisto(:);
Resist(:,2) = Resistc(:);
RT(:,1) = RTo(:);
RT(:,2) = RTc(:);
Symb = {'o','d'};
Col = autumn(length(FLEPs));
Cmap = flipud(Col);

IB_plot(Resist,RT,T-Ti,Symb,Col,Cmap,(1-FLEPs),'Fishing depletion')
text(Resist(1,1)*0.9,RT(1,1),'Open','color','r','fontsize',12)
text(Resist(1,2)*0.9,RT(1,2),'Closed','color','r','fontsize',12)
end % end if Fig5
%%




%% Fig 6 - Pacific cod, fishing stops
if Fig6
    Species = 'Black rockfish';
    % Mean F of 0.5803 in 10 years prior to 2014, that corresponds to this
    % FLEP:(Table 2.24 in 2018 assessment)
    FLEPs = 0.2:0.05:0.95; %0.28; % actual estimated value for Pacific cod in the Barbeaux assessement: 0.2864; 
    %Reduction: 58% of normal for 3 years (Laurel & Rogers 2000)
    % age-0 recruitment in past 10 years
    %Age0 = [1.002, 1.703,0.989, 0.678, 0.49, 0.942, 0.761, 0.869, 0.734, 0.372];
    Reduction = log(1-0.58);
    Duration = 3;
    Delay = NaN; % only used if multiple disturbances (it gives the interval)
    T = 200; % total time
    Ti = 100; % time prior to disturbance
    Conn_scenario = 'Closed'; % Connectivity: open (linear subsidy) or closed
    DD_scenario = 'BH';
    F_reduct = 0:0.1:1;
    
    figure(1)
    clf
    set(gcf,'units','cent','position',[10 10 19 20])
    
    
    Cols = flipud(winter(length(F_reduct)));
   
    for ff = 1:length(FLEPs)
    for f = 1:length(F_reduct)
    
    [~, N, B, ~, ~, Y, ~,Resist,RT] = square_pulse_resil(Species,FLEPs(ff),Reduction,Duration,Delay,T,Ti,...
                                                     Conn_scenario, DD_scenario,'B',0.95,1,false,F_reduct(f));
                        
  %keyboard                                               
    Bplot_tmp = sum(squeeze(B));
    Bplot_tmp = Bplot_tmp/Bplot_tmp(Ti-1);
    Bplot(:,f,ff) = Bplot_tmp(:); 
    
    Yplot_tmp = Y/Y(Ti-1);
    Yplot(:,f,ff) = Yplot_tmp(:);
    
        
        Ym(f,ff) = mean(Y(Ti:(Ti+9)))/Y(Ti-1);
        
        
        Ymax(ff) = max(Y);
    end
    end
    
    % standardize yield to be relative
   % Ym = Ym/Ym(1);
    Whichone = 2;
    Best = find(Ym(:,Whichone)==max(Ym(:,Whichone)),1);
    
    
    Col = winter(length(F_reduct));
    
    subplot(3,1,1)
    resil_trajec_plot(Bplot(:,:,Whichone),Duration,T,Ti,'Relative biomass',Col)
   
    plot(-5:100,Bplot((Ti-5):end,Best,Whichone),'k-','linewidth',2)
    xl = get(gca,'xlim');
    set(gca,'xlim',[-5,20])
    
    colormap(flipud(Cols))
    cb = colorbar;
    set(cb,'tickdir','out','ticklength',0.015)
    caxis([F_reduct(1),F_reduct(end)])
    ylabel(cb,'Prop. reduction in F','rotation',270,'verticalalignment','bottom','fontsize',12)
    
    
    subplot(3,1,2)
    resil_trajec_plot(Yplot(:,:,Whichone),Duration,T,Ti,'Relative yield',Col)
    plot(-5:100,Yplot((Ti-5):end,Best,Whichone),'k-','linewidth',2)

    xl = get(gca,'xlim');
    set(gca,'xlim',[-5,20])
    
   % xlim([99 115])
    colormap(flipud(Cols))
    cb = colorbar;
    set(cb,'tickdir','out','ticklength',0.015)
    caxis([F_reduct(1),F_reduct(end)])
    ylabel(cb,'Prop. reduction in F','rotation',270,'verticalalignment','bottom','fontsize',12)
    
    
    subplot(3,1,3)
    pcolor(1-FLEPs,F_reduct,Ym); shading interp
    ylabel('Prop. reduction in F','fontsize',14)
    xlabel('Depletion','fontsize',14)
    set(gca,'xcolor','k','ycolor','k','tickdir','out','ticklength',[0.02 0.02])
    %ylabel('Mean yield')
    xlim([0.15 0.8])
    cb = colorbar;
    ylabel(cb,'Relative yield','rotation',270,'verticalalignment','bottom','fontsize',12)
    set(cb,'ycolor','k','tickdir','out','ticklength',[0.015])
                             
   
                                              
end % end if Fig6
%%

%% Fig 6b - Pacific cod w & w/o recruitment
if Fig6b
    Species = 'Pacific cod';
    % Mean F of 0.5803 in 10 years prior to 2014, that corresponds to this
    % FLEP:(Table 2.24 in 2018 assessment)
    FLEPs = 0.26; %0.28; % actual estimated value for Pacific cod in the Barbeaux assessement: 0.2864; 
    %Reduction: 58% of normal for 3 years (Laurel & Rogers 2000)
    % age-0 recruitment in past 10 years
    %Age0 = [1.002, 1.703,0.989, 0.678, 0.49, 0.942, 0.761, 0.869, 0.734, 0.372];
    Reduction = [log(0.99999999), log(1-0.58), log(1-0.58)];
    Duration = 3;
    Delay = NaN; % only used if multiple disturbances (it gives the interval)
    T = 200; % total time
    Ti = 100; % time prior to disturbance
    Conn_scenario = 'Closed'; % Connectivity: open (linear subsidy) or closed
    DD_scenario = 'BH';
    F_reduct = 0; %:0.1:1;
    AdultMort = [exp(0.499 - 0.8448);... % difference factor between two natural mortality rates
                 exp(0.499 - 0.8448);...
                 1];
    figure(1)
    clf
    set(gcf,'units','cent','position',[10 10 19 20])
    
    
    %Cols = flipud(winter(length(F_reduct)));
    Cols = [0.1 0.1 0.9; 0 0 0;  0.9, 0.1 0.1];
   
    for ff = 1:length(Reduction)
    
    [~, N, B, ~, ~, Y, ~,Resist,RT] = square_pulse_resil(Species,FLEPs,Reduction(ff),Duration,Delay,T,Ti,...
                                                     Conn_scenario, DD_scenario,'B',0.95,1,true,F_reduct,AdultMort(ff));
                        
  %keyboard                                               
    Bplot_tmp = sum(squeeze(B));
    Bplot_tmp = Bplot_tmp/Bplot_tmp(Ti-1);

    Bplot(:,ff) = Bplot_tmp(:); 
    
    Yplot_tmp = Y/Y(Ti-1);
    Yplot(:,ff) = Yplot_tmp(:);
    
        
        Ym(ff) = mean(Y(Ti:(Ti+9)))/Y(Ti-1);
        
        
        Ymax(ff) = max(Y);
    end
   % keyboard
    
    % standardize yield to be relative
   % Ym = Ym/Ym(1);
  %  Whichone = 2;
  %  Best = find(Ym(:,Whichone)==max(Ym(:,Whichone)),1);
    
    
 %   Col = winter(length(F_reduct));
    
    subplot(2,1,1)
    hold on
    resil_trajec_plot(Bplot,Duration,T,Ti,'Relative biomass',Cols)
   
  %  plot(-5:100,Bplot((Ti-5):end,Best,Whichone),'k-','linewidth',2)
    xl = get(gca,'xlim');
    set(gca,'xlim',[-5,20])
    
  %  colormap(flipud(Cols))
  %  cb = colorbar;
  %  set(cb,'tickdir','out','ticklength',0.015)
  %  caxis([F_reduct(1),F_reduct(end)])
  %  ylabel(cb,'Prop. reduction in F','rotation',270,'verticalalignment','bottom','fontsize',12)
    
    
    subplot(2,1,2)
    resil_trajec_plot(Yplot,Duration,T,Ti,'Relative yield',Cols)
  %  plot(-5:100,Yplot((Ti-5):end,Best,Whichone),'k-','linewidth',2)

    xl = get(gca,'xlim');
    set(gca,'xlim',[-5,20])
    
   % xlim([99 115])
 %   colormap(flipud(Cols))
 %   cb = colorbar;
 %   set(cb,'tickdir','out','ticklength',0.015)
   % caxis([F_reduct(1),F_reduct(end)])
   % ylabel(cb,'Prop. reduction in F','rotation',270,'verticalalignment','bottom','fontsize',12)
    
                                 
end % end if Fig6b
%%


%% Loo's Convolution Figure
if FigS1
Species = 'Black rockfish';
FLEP = 1.0; % from Dick et al assessment
% Check that maturity ogive is specified correctly...

convo_plot(Species,FLEP)
end %
%% End FigS1

%% Fig 7 - Pacific cod (deprecated)
if Fig7
    Species = 'Black rockfish';
    % Mean F of 0.5803 in 10 years prior to 2014, that corresponds to this
    % FLEP:(Table 2.24 in 2018 assessment)
    FLEPs = 0.2864;
    %Reduction: 58% of normal for 3 years (Laurel & Rogers 2000)
    % age-0 recruitment in past 10 years
    %Age0 = [1.002, 1.703,0.989, 0.678, 0.49, 0.942, 0.761, 0.869, 0.734, 0.372];
    Reduction = log(1-0.58);
    Duration = 3;
    Delay = NaN; % only used if multiple disturbances (it gives the interval)
    T = 200; % total time
    Ti = 100; % time prior to disturbance
    Conn_scenario = 'Open'; % Connectivity: open (linear subsidy) or closed
    DD_scenario = 'None';
    
    
    [~, N, B, ~, ~, ~, ~,Resist,RT] = square_pulse_resil(Species,FLEPs,Reduction,Duration,Delay,T,Ti,...
                                                     Conn_scenario, DD_scenario,'B',0.95,true);

                                                 Nplot = sum(squeeze(N));
Nplot = Nplot/Nplot(Ti-1);

Bplot = sum(squeeze(B));
Bplot = Bplot/Bplot(Ti-1);

% Consider specifying plot window size, etc.
resil_trajec_plot([Nplot(:),Bplot(:)],Duration,T,Ti,'Relative biomass')
% Add actual values:
Time = 2013:2019;
Time = Time-2015;
B = [0.9958270 1.0000000 0.7029892 0.5465771 0.4048187 0.3282937 0.3095126]; % relative biomass - from Suryan et al & Barbeaux assessment
plot(Time,B,'ko')

% use cod scenario, reduce F by various % for different lengths of time


end % end if Fig7
%%


