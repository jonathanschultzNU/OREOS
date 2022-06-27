function beatcontours(fig,freqinds,i,limcor, params, xdata, ydata, zdata,w2val,cmap)
%CREATEFIGURE(xdata1, ydata1, zdata1, XData2)
%  XDATA1:  contour x
%  YDATA1:  contour y
%  ZDATA1:  contour z

cmap = cmap(size(cmap,1)/2:end,:);

zdata = normdim(zdata);
ax = axes();
hold(ax);
subplot2 = subplot(1,length(freqinds),i,'Parent',fig);
hold(subplot2,'on');
p = pcolor(xdata,ydata,zdata);
p.FaceAlpha = 0.8;

colormap(cmap)
shading interp
colorbar;
clist = 0.1:0.1:1; 
c = contour(xdata,ydata,zdata,'k','LineWidth',0.25,'LevelList',clist);

hold on

% Create lines
line('XData',[min(ydata) max(ydata)],'YData',[min(ydata) max(ydata)],'Linewidth',1.5);
hold on
yline(params.e1/1000, 'k--', 'LineWidth', 0.5); hold on
xline(params.e1/1000, 'k--', 'LineWidth', 0.5); hold on

% Define limits
xlim([min(xdata)+limcor max(xdata)-limcor])
ylim([min(ydata)+limcor max(ydata)-limcor])

ylabel('\omega_{3} (10^{3} cm^{-1})');
xlabel('\omega_{1} (10^{3} cm^{-1})');

% Create title
title(strcat('\omega_{2} =', num2str(w2val), ' cm^{-1}'));

box on
set(gca,'BoxStyle','full','CLim',[0 1],'DataAspectRatio',[1 1 1],...
    'FontSize',18,'Layer','top');
grid minor

