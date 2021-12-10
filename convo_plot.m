function convo_plot(Species,FLEP)

% Plot the four influence functions for a species

Params = define_Params(Species);
F = get_F(Params,FLEP);
T = 0:200;
Ti = 50;
Duration = 21;
Dist = ones(size(T));
Dist(Ti:(Ti+Duration-1)) = 0.1;
Indices = [Ti-8, Ti+10, Ti+Duration-1+2,Ti+Duration-1+10];
Col = winter(length(Indices));

% Simulate the disturbance
[~, N, B]=square_pulse_resil(Species,FLEP,log(0.1),Duration,NaN,T(end),Ti,'Open', 'None');

% SAD + effect of disturbance
for i = 1:length(Indices)
   AgeX = (Indices(i)-(Duration-1)-Ti+1):(Indices(i)-Ti); % +1 accounts for us starting ages at 1 not 0
   
   AgeX = AgeX(AgeX > 0 & AgeX <= Params.A);
   
   SADx(:,i) = Params.SAD;
   SADx(AgeX,i) = SADx(AgeX,i)*Dist(Ti);
end


figure(1)
clf
set(gcf,'units','cent')
%set(gcf,'position',[10 10 20 18])
set(gcf,'position',[10 10 20 24])
LW = 2;


for i = 1:length(Indices)
    %A(1) = subplot(3,1,1);
 A(1) = subplot(4,1,1);
hold on
%plot(T+0.5,Dist,'k-','linewidth',LW)
patch([0 0 Duration-1 Duration-1]+Ti,[0 1 1 0],[-1 -1 -1 -1],...
    'facecolor',[0 0 0],'facealpha',0.1,'edgecolor','none')
    % vector of x-indices
    Xtmp = Indices(i)-(Params.Ages-1); % -1 so starts at 0
plot(Xtmp,Params.SAD,'k:','linewidth',LW,'color',Col(i,:))
%plot(Xtmp,SADx(:,i),'k-','linewidth',LW,'color',Col(i,:))
area(Xtmp,SADx(:,i),'linewidth',LW,'facecolor',Col(i,:),'facealpha',0.5)
plot(Xtmp,SADx(:,i),'ko','markersize',4)
end

ylabel(A(1),'Relative abundance','fontsize',14);


%A(2) = subplot(3,1,2);
A(2) = subplot(4,1,2);
%A(5) = subplot(2,4,5:8);
hold on

% Set time interval to be plotted
Tgap = 20;
Tv = (Ti-Tgap):(Ti+length(T)/2);
Tv0 = Tv - Ti;
Indices2 = Indices-Ti;
X = sum(squeeze(N));
X = X/X(1);
%keyboard
plot(Tv0,X(Tv),'k','linewidth',2,'linestyle','-')
plot(Tv0,X(Tv),'ko','markersize',4)
patch([0 0 Duration-1 Duration-1],[0 1 1 0],[-1 -1 -1 -1],...
    'facecolor',[0 0 0],'facealpha',0.25,'edgecolor','none')
for i = 1:length(Indices)
plot(Indices2(i),X(Indices(i)),'ko','markerfacecolor',Col(i,:),'markersize',8)
end
ylabel('Total abundance','fontsize',14)
%%%%%%%%%%%%%%%%%%%%%

%A(3) = subplot(3,1,3);
A(3) = subplot(4,1,3);
hold on

% Set time interval to be plotted
Tgap = 20;
Tv = (Ti-Tgap):(Ti+length(T)/2);
Tv0 = Tv - Ti;
Indices2 = Indices-Ti;
X = sum(squeeze(B));
X = X/X(1);
%keyboard
plot(Tv0,X(Tv),'k','linewidth',2,'linestyle','-')
plot(Tv0,X(Tv),'ko','markersize',4)
patch([0 0 Duration-1 Duration-1],[0 1 1 0],[-1 -1 -1 -1],...
    'facecolor',[0 0 0],'facealpha',0.25,'edgecolor','none')
for i = 1:length(Indices)
plot(Indices2(i),X(Indices(i)),'ko','markerfacecolor',Col(i,:),'markersize',8)
end
ylabel('Total biomass','fontsize',14)
%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mean ages
Mean_age_infl = sum(Params.Ages.*Params.SAD)/sum(Params.SAD);
Mean_biomass_infl = sum(Params.Ages.*Params.BiomassAge)/sum(Params.BiomassAge);

% get response time
Thresh = 0.95;
B0 = squeeze(B);
B0 = sum(B0)/sum(B0(:,Ti-1));
N0 = squeeze(N);
N0 = sum(N0)/sum(N0(:,Ti-1));
BRtime = find(B0<Thresh,1,'last')-Ti-Duration;
NRtime = find(N0<Thresh,1,'last')-Ti-Duration;

disp(strcat('Mean of age influence',num2str(Mean_age_infl)))
disp(strcat('Abundance recovery time',num2str(NRtime)))
disp(strcat('Mean of biomass influence',num2str(Mean_biomass_infl)))
disp(strcat('Biomass recovery time',num2str(BRtime)))


%%%%%%%%%%%%%%%%%%%%
%%%%%%%
% Bubble plot to show cohorts
A(4) = subplot(4,1,4);
hold on
N = squeeze(N);
% convert from area to radius
N = N.^2
Ages = 1:(size(N,1));
Years = 1:(size(N,2));

% expand to full dimensions for plotting
Ages2 = repmat(Ages(:),[size(N,2),1]);
Years2 = repmat(Years(:)',[size(N,1),1]);
Years2 = Years2(:);


bubblechart(Years2,Ages2,N(:),'markeredgecolor','k','markerfacecolor','none','linewidth',1)
bubblesize([1e-3,10])
patch([50 50 50+Duration-1 50+Duration-1],[0 15 15 0],[-1 -1 -1 -1],...
    'facecolor',[0 0 0],'facealpha',0.25,'edgecolor','none')
for i = 1:length(Indices)
bubblechart(repmat(Indices(i),[length(Ages),1]),Ages,N(:,Indices(i)),'markeredgecolor','k','markerfacecolor',Col(i,:))
end
bubblesize([1e-3,10])

%%%%%%%%%%%%%%%%
%%%% Formatting
set(A(:),'xcolor','k','ycolor','k','tickdir','out','ticklength',[0.015 0.015])
set(A(1),'xlim',[30 85],'xtick',10:10:100,'xticklabels',-40:10:100,'ytick',[],'ylim',[-0.03 1.03])
set(A(2:3),'xlim',[-20 35],'xtick',-20:10:100,'ytick',0:0.2:1)
xlabel(A(4),'Time (y)','fontsize',14)
set(A(4),'xlim',[30 85],'ylim',[0 15],'xtick',30:10:80,'xticklabels',-20:10:30)
ylabel(A(4),'Age (y)','fontsize',14)










