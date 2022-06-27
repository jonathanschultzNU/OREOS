function fig1 = fig_open()
 set(0,'defaultfigurecolor',[1 1 1]) % white background
 set(0,'defaultaxesfontname','cambria math') % beautify the axes a bit
 scrn_size = get(0,'ScreenSize'); % get size of screen
 shrink_pct = 0.15; % shrink the figure by 10%
 %
 fig1 = figure('Visible','on','DefaultAxesFontSize',20,'Position',...
     [scrn_size(1)+(scrn_size(3)*shrink_pct) scrn_size(2)+(scrn_size(4)*shrink_pct)...
     scrn_size(3)-(scrn_size(3)*2*shrink_pct) scrn_size(4)-(scrn_size(4)*2*shrink_pct)]); % shrinking the figure
 %
 set(gca,'FontSize',20,'YTickLabel',[]);

 grid on % grid for more contrast in the axes
 end