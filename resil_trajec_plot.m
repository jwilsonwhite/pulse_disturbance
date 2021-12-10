function resil_trajec_plot(X,Duration,T,Ti,YLabel,Col,fname)

% Plots the trajectory of a disturbance and recovery
if exist('fname','var')
   figure
end

% include option for multipanel

%LS = {'-','--',':','-.'}; % possible line styles
LS = cell(1,size(X,2)); 
LS(:) = {'-'};

if ~exist('Col','var')
    Col = repmat([0 0 0],[size(X,2),1]);
end

% Set time interval to be plotted
Tv = (Ti-Duration):(Ti+100);
Tv0 = Tv - Ti;

hold on

% If there are multiple vectors for plotting
if any(size(X)==1)
plot(Tv0,X(Tv),'b','linewidth',2)
else
    
    % ensure correct dimensionality
    if size(X,1)<size(X,2); X = X'; end
    
    for i = 1:size(X,2)
        plot(Tv0,X(Tv,i),'b','linewidth',2,'linestyle',LS{i},'color',Col(i,:))
        plot(Tv0,X(Tv,i),'ko','markeredgecolor',Col(i,:),'markersize',4)
        %keyboard
    end 
end % end if column vector
    

% Plot the 'disturbance' interval behind
patch([0 0 Duration-1 Duration-1],[0 1 1 0],[-1 -1 -1 -1],...
    'facecolor',[0 0 0],'facealpha',0.25,'edgecolor','none')
%plot(Tv0,Noise(Tv),'k')


% Plot formatting
set(gca,'ylim',[0 1.1],'xlim',[-5 100],'xtick',0:10:100)
set(gca,'fontsize',10)
set(gca,'xcolor','k','ycolor','k','tickdir','out','ticklength',[0.015 0.015])
xlabel('Time','fontsize',14)
ylabel(YLabel,'fontsize',14)
set(gca,'ygrid','on','xgrid','on')

% Save
if exist('fname','var')
   print -depsc2 fname 
end

