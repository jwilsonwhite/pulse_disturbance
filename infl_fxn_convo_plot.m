function infl_fxn_convo_plot(Species,FLEP,Mfactor)

% Plot the pulse response as a convolution of influence function &
% disturbacne



T = 80;
Duration = 10;
T1 = zeros(T,1);
T1(20:(20+Duration-1)) = 0.75;
T2 = zeros(T,1);
T2a = normpdf(19:29,29,Duration/4)*6 + normrnd(0,0.2,[11,1])';
T2b = normpdf(30:40,29,Duration/8)*6+ normrnd(0,0.2,[11,1])';
T2(19:29) = T2a/max(T2a);
T2(30:40) = T2b/max(T2b);
T2 = max(T2,0);


Cols = {'k',[0.1 0.7 0.1],[0.7 0 0.8]};

figure(1)
clf
set(gcf,'units','cent')
set(gcf,'position',[10 10 18 16])
LW = 2;

SPs = [1 2; 3 4; 5 6];

for s = 1:2
    
Params_master = define_Params(Species{s});
F = get_F(Params_master,FLEP);
if ~exist('Mfactor','var')
    Mfactor = 1;
end
    
    
for m = 1:length(Mfactor)
% Calculate fishery yield pattern
Params = define_Params(Species{s},Mfactor(m));
%L = get_Leslie(Params,F,'Open',NaN,NaN);
%Lf = diag(L,-1);
%Lf = cumprod(Lf); % cumulative survival
%Nf = Params.BiomassAge(:).*Params.isFish(:).*[Lf(:);0].*(F/(F+Params.M));


A(1)=subplot(3,2,SPs(1,s));
hold on
plot(Params.Ages,Params.SAD/max(Params.SAD),'k-','linewidth',LW,'color',Cols{m})
plot(Params.Ages,Params.SAD.*Params.BiomassAge/max(Params.SAD.*Params.BiomassAge),'k--','linewidth',LW,'color',Cols{m})
ylabel({'Relative abundance or', 'relative biomass'},'fontsize',14);
xlabel('Age','fontsize',14);
set(gca,'ytick',[])

A(2)=subplot(3,2,SPs(2,s));
hold on
plot(1-T1,'k-','linewidth',LW,'color',Cols{2})
plot(1-T2,'k-','linewidth',LW,'color',Cols{3})
ylabel('Recruit survival','fontsize',14);
xlabel('Time','fontsize',14);
ylim([-0.25 1]);
xlim([15 70])
set(gca,'xtick',20:10:80,'xticklabels',0:10:100,'xgrid','on','ytick',[])


A(3)=subplot(3,2,SPs(3,s));
hold on
convA1 = conv(T1,Params.SAD);
convA1 = convA1/max(convA1);
convB1 = conv(T1,Params.SAD.*Params.BiomassAge);
convB1 = convB1/max(convB1);
plot(1:(length(convA1)),1-convA1,'k-','color',Cols{2},'linewidth',LW);
plot(1:(length(convB1)),1-convB1,'k--','color',Cols{2},'linewidth',LW);
convA2 = conv(T2,Params.SAD);
convA2 = convA2/max(convA2);
convB2 = conv(T2,Params.SAD.*Params.BiomassAge);
convB2 = convB2/max(convB2);
plot(1:(length(convA2)),1-convA2,'k-','color',Cols{3},'linewidth',LW);
plot(1:(length(convB2)),1-convB2,'k--','color',Cols{3},'linewidth',LW);
xlim([15 70])
ylabel({'Total abundance or', 'Biomass'},'fontsize',14);
xlabel('Time','fontsize',14);
set(gca,'xtick',20:10:80,'xticklabels',0:10:100,'xgrid','on','ytick',[])




end % end loop over Mfactor
end % end loop over species


set(A(:),'xcolor','k','ycolor','k','tickdir','out','ticklength',[0.015 0.015])





