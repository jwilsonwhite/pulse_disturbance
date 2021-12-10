function infl_fxn_plot(Species,FLEP,Mfactor)

% Plot the four influence functions for a species, also permitting
% variation in M

Params_master = define_Params(Species);
F = get_F(Params_master,FLEP);
if ~exist('Mfactor','var')
    Mfactor = 1;
end

Cols = {'b','k','r'};

figure
set(gcf,'units','cent')
set(gcf,'position',[10 10 9 18])
LW = 2;


for m = 1:length(Mfactor)
% Calculate fishery yield pattern
Params = define_Params(Species,Mfactor(m));
L = get_Leslie(Params,F,'Open',NaN,NaN);
Lf = diag(L,-1);
Lf = cumprod(Lf); % cumulative survival
Nf = Params.BiomassAge(:).*Params.isFish(:).*[Lf(:);0].*(F/(F+Params.M));


A(1)=subplot(4,1,1);
hold on
plot(Params.Ages,Params.SAD,'k-','linewidth',LW,'color',Cols{m})
ylabel('Abundance','fontsize',14);

A(2)=subplot(4,1,2);
hold on
plot(Params.Ages,Params.SAD.*Params.BiomassAge,'k-','linewidth',LW,'color',Cols{m})
ylabel('Biomass','fontsize',14);

A(3)=subplot(4,1,3);
hold on
plot(Params.Ages,Params.SAD(:).*Params.EP0(:),'k-','linewidth',LW,'color',Cols{m})
ylabel('Fecundity','fontsize',14);

A(4)=subplot(4,1,4);
hold on
plot(Params.Ages,Nf,'k-','linewidth',LW,'color',Cols{m})
ylabel('Fishery yield','fontsize',14);


end % end loop over Mfactor

set(A(:),'xcolor','k','ycolor','k','tickdir','out','ticklength',[0.015 0.015],...
    'xlim',[1 Params.A],'xtick',[1,5:5:100],'ytick',[])
xlabel(A(4),'Age (y)','fontsize',14)




