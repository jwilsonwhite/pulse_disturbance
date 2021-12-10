function IB_plot(Resist,RT,T,Symb,Col,Cmap,Cval,Clab)

% Resist = resistance (vector)
% RT = return time (vector)
% T = max time (scalar)

if ~exist('Symb','var')
    Symb = cell(size(Resist,2));
    Symb{:} = 'o';
end
if ~exist('Col','var')
    Col = zeros(size(Resist,1),3);
end


figure
hold on

%keyboard
for i = 1:size(Resist,2)
    for j = 1:size(Resist,1)
plot(1-Resist(j,i),RT(j,i),'marker',Symb{i},...
    'markerfacecolor',Col(j,:),'markeredgecolor',Col(j,:),...
    'markersize',8)
    end
end


% Plot formatting
set(gca,'ylim',[0 T],'xlim',[0 1],'xtick',0:0.1:1) %,'xdir','reverse')
set(gca,'fontsize',10,'color',repmat(0.8,1,3))
set(gca,'xcolor','k','ycolor','k','tickdir','out','ticklength',[0.015 0.015])
xlabel('Impact (proportional decline)','fontsize',14)
ylabel('Return time (y)','fontsize',14)
%set(gca,'ygrid','on','xgrid','on')


colormap(Cmap)
cb = colorbar;
set(cb,'tickdir','out','ticklength',0.015)
caxis([min(Cval),max(Cval)])
ylabel(cb,Clab,'fontsize',14,'rotation',270,'verticalalignment','bottom')
keyboard
